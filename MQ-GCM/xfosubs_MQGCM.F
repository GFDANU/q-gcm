c***********************************************************************
c     Q-GCM Version 1.5.0 : last modified 10/08/2020, by S. Kravtsov,
c     to include additions associated with the moist version of Q-GCM
c     (this option is invoked by the flag "moist"). In this model, the
c     variables hmixa,  hmixam denote the specific humidity in the atmospheric
c     mixed layer,  whose thickness (previously stored in these variables)
c     is now assumed constant. The code has an additional function esat
c     which computes saturation specific humidity as a function of temperature.
c     The latent heat loss by the ocean is computed here, but the latent heat
c     forcing of the atmosphere - later on in amlsubs, after the precipiation
c     has been computed from moisture equations. Additional shared variables 
c     evap and evapa have also been added to the intrfac module
c============================================================================
c
c     Modified by Sergey Kravtsov to open the S/N boundaries in the a.m.l.
c
c     This version goes with the modified radiation scheme and a.m.l.
c     emissivity (zm) < 1
c
c     These changes were made in 07-08/2020.
c
c     Before this, the following  changes were also made:
c 
c     12/28/2017: Previously added dependence of wind-stress on 
c     mixed-layer temperature difference, as per  Hogg et al. (2009)
c     (in particular, an 
c     additional interpolation of atmospheric
c     temperature to ocean resolution)
c
c=======================================================================
c     Also previousoy modified by S. Kravtsov to include temperature-driven
c     wind (via hydrostatic adjustment) in the atmospheric mixed layer. 24-07-2018
c
c     ***Both of the temperature dependent options above are invoked, 
c        simultaneously, by the flag temp_fdbck*** (09/09/2020)
c
c     The individual execution of the two options is controlled below by
c     the parameters alp and alp1. 
c
c=======================================================================
c
c=========================================================================
c    Also  included  here is the ability of the ocean to "see", instead
c    of the  full  wind stress, its temperature independent
c    version (Kravtsov, 02/08/2019) -- a partially coupled case --
c    with a possibility to run the code in the fully coupled mode as well
c
c    ***This option is invoked by the flag partial_coup*** (09/09/2020)
c
c=========================================================================
c
c***********************************************************************
c
c     Copyright 2013 Jeff Blundell, Andy Hogg and Bill Dewar.
c     This file is part of Q-GCM.
c
c     Q-GCM is free software: you can redistribute it and/or modify
c     it under the terms of the GNU General Public License as
c     published by the Free Software Foundation, either version 3
c     of the License, or (at your option) any later version.
c
c     Q-GCM is distributed in the hope that it will be useful,
c     but WITHOUT ANY WARRANTY; without even the implied warranty
c     of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
c     See the GNU General Public License for more details.
c
c     You should have received a copy of the GNU General Public License
c     along with Q-GCM.  If not, see <http://www.gnu.org/licenses/>.
c
c***********************************************************************
c
      MODULE xfosubs

*     Contains subroutines xforc, bilint, auvbcu, bcuini and wts2bb,
*     and the atmospheric radiative forcing shape function fsprim, of
*     which only xforc and fsprim are called from outside this module.
*     These routines compute the forcings: stresses (7.1-7.5), Ekman
*     velocities on both T and p grids (7.6-7.7), and net diabatic
*     forcings of oceanic and atmos. mixed layers (7.8-7.10).

*     Modules
      USE parameters, ONLY : ndxr

      IMPLICIT NONE

      PRIVATE

      PUBLIC  :: xforc, fsprim, esat
      PRIVATE :: bilint, auvbcu, bcuini, wts2bb

*     Storage for bicubic interpolation weights on high-resolution grid
      double precision, PRIVATE, SAVE ::
     &                 stbbb(16,0:ndxr,0:ndxr),
     &                 stbus(16,0:ndxr,0:ndxr),stbun(16,0:ndxr,0:ndxr),
     &                 stbvs(16,0:ndxr,0:ndxr),stbvn(16,0:ndxr,0:ndxr)

      CONTAINS

c***********************************************************************
c
      SUBROUTINE xforc
*
*     Compute forcings: stresses (7.1-7.5), Ekman velocities
*     on both T and p grids (7.6-7.7), and net diabatic
*     forcings of oceanic and atmos. mixed layers (7.8-7.10).
*
*     Modified at v1.4.0 to compute tau at oceanic resolution, to
*     enable tau optionally to depend on ocean mixed layer velocity
*     (depending on setting of preprocessor option tau_udiff)

*     Modules
      USE parameters
      USE atconst
      USE occonst
#ifndef ocean_only
      USE athomog, ONLY : txisat, txinat
#endif
#if defined tau_udiff || !defined ocean_only
      USE atstate, ONLY : pam, wekpa, wekta
#endif
#ifndef atmos_only
      USE ochomog
      USE ocstate, ONLY : wekpo, wekto
#endif
      USE intrfac, ONLY : sstm, astm, tauxo, tauyo, tauxa, tauya,
     &                    uekat, vekat, uekat1, vekat1, 
     &                    fnetoc, fnetat, hmixam,
     &                    hmadmp, hmoc, hmat, cdat, raoro, ssta,
     &                    evapa,evap,Ch,Ce
      USE radiate, ONLY : xlamda, Adown, Bmup, B1down,
     &                    Cmup, C1down, D0up, Dmup, Dmdown,zm,
     &                    gamma,Flambar
      USE monitor, ONLY : arlaav, slhfav, oradav, arocav
*     Optionally required for checking
c      USE intsubs, ONLY : xintt,xintp

      IMPLICIT NONE

*     Subroutine arguments
*
*     Local parameters

      integer iocoff,jocoff,ipobeg,ipoend,jpobeg,jpoend,jsou,jnor
      parameter ( iocoff = (nx1-1)*ndxr, jocoff = (ny1-1)*ndxr,
     &            ipobeg = iocoff + 1, ipoend = ipobeg + nxto,
     &            jpobeg = jocoff + 1, jpoend = jpobeg + nyto,
     &            jsou = 1 + ndxr/2, jnor = nypaor - ndxr/2 )
*     These are the offsets of ocean p points in the ocean-resolution
*     atmospheric arrays, and their start and end indices therein; also
*     the boundaries of the cells surrounding the timestepped q points
*

#  ifdef temp_fdbck
c---------------------------------------------------------------------
c     The parameter  alp below governs the strength of the
c     mixed-layer temperature difference feedback, with
c     alp=0 corresponding to the default case of no feedback;
c     The suggested value for the simulation with the feedback
c     is 0.15 (cf. Hogg et al. 2009)
c
c     The parameter alp1 should be set to either 0.d0 or 1.d0.
c     alp1=0 corresponds to the original formulation, without
c     hydrostatic adjustment of mixed-layer wind to the temperature
c     field. alp1=1 is the case with hydrostatic  adjustment.
c
c
c     Sergey Kravtsov, 09/06/2020
c
      double precision alp, alp1
      parameter ( alp = 0.15d0, alp1 = 1.d0 )

#  endif
*---------------------------------------------------------------------
*
*     Local variables
      integer i,j,ioff,joff,io,jo
      double precision hxafac,hxofac,txsums,txsumn
      double precision radeta(nxta,nyta),radeto(nxto,nyto)
*     Optional variables for extra checking
c      double precision wekint,wpsumi

#if defined tau_udiff || !defined ocean_only
      integer nijwid,ia,ja,ibeg,jbeg,jlo,jhi,it
      double precision cdhfaa,cdhfab,cdrfaa,cdrfab,qu2faa,qu2fab,
     &                 cdhfaa1, qu2faa1,
     &                 zbfcat,zbfcoc,u1at(nxpa,nypa),v1at(nxpa,nypa),
     &                 u1ator(nxpaor,nypaor),v1ator(nxpaor,nypaor),
     &                 cdrfac(nxpaor),qu2fac(nxpaor),delu1,delv1,delt,
     &                 scasqd,scashr,cdochi,tausum,uvekfc,hmrdxa,
     &                wt(0:ndxr),cdhfac(nxpaor),wspdor(nxpaor,nypaor),
     &                 wspdot(nxto,nyto),
     &                 tauxaor(nxpaor,nypaor),tauyaor(nxpaor,nypaor),
     &                 wektaor(nxtaor,nytaor),wtasum,asto(nxto,nyto),
     &              wsum,wtj,fsp,ocfrac,fmafac,fmafac1,fmatop,hmafac,
     &                 ocnrad,slhf,atmrad,arlasm,arocsm,slhfsm,oradsm

#  ifdef temp_fdbck
      double precision astfac,astp(nxpa,nypa),astpor(nxpaor,nypaor),
     &                 sstp(nxpo,nypo),astpor1(nxpaor,nypaor)
#   ifdef partial_coup
      double precision u1at1(nxpa,nypa),v1at1(nxpa,nypa),
     &                 u1ator1(nxpaor,nypaor),v1ator1(nxpaor,nypaor),
     &               tauxaor1(nxpaor,nypaor),tauyaor1(nxpaor,nypaor)
#   endif
#  endif

#  ifdef moist
      double precision hmixao(nxto,nyto),sstr,hsat,sstr1,hsat1,hsat2,
     &                 fact,evapm
      double precision Gc,g,RR,RV,eeps,Gtilde
      parameter (
     &           Gc=6.5d-3,g=9.82d0,RR=287.0d0,RV=461.d0,eeps=RR/RV,
     &           Gtilde=g/(RR*Gc))

#  endif

*     Optional variables for extra checking
!!    double precision uvint
      logical ndxodd
#  ifdef tau_udiff
      double precision u1oc,v1oc
#  endif
*     Extra variables for MONITOR
      integer natlan,natocn
#endif

      hxafac = 0.5d0*rdxaf0
      hxofac = 0.5d0*rdxof0
#ifdef temp_fdbck
      astfac= alp1*0.5d0*9.82d0*hmat/300.d0
#endif

*======================================================================
*     9.82 is the gravity acceleration, 300.d0 is the
*      reference potential temperature for the atmospheric mixed layer.
*     Setting astfac to zero would remove temperature dependent flow
*     in the atmospheric mixed layer, returning back to the original
*     model formulation. This should be done by setting the parameter alp1=0.
*     S. Kravtsov, 24-07-2018
*======================================================================


#if defined tau_udiff || !defined ocean_only

      uvekfc = 1.0d0/( hmat*fnot*dble(ndxr) )
      hmrdxa = hmat/dxa

*     Work out factors in the windstress equations for both
*     atmos. vel. only (aa) and velocity shear (ab) cases
*     -----------------------------------------------------
*     cdhfab = "(a+b)" in the Userguide, cdhfaa = "a"
*     Define both; decide which to use later in code
c
c-----------------------------------------------------------------------
c     Added weaker drag coefficient over ocean (2/3 of that  over land)
c     this modifies cdhfaa1, cdhfab (and qu2fab), as well as qu2faa1
c
c     Sergey  Kravtsov (09/04/2020)
c----------------------------------------------------------------------
      cdhfaa = (cdat/fnot)/hmat
      cdhfaa1 = (2.d0/3.d0)*cdhfaa
      cdhfab = (2.d0/3.d0)*(cdat/fnot)*( 1.0d0/hmat + raoro/hmoc )
      cdrfaa = cdat/( abs(cdhfaa) )
      cdrfab = (2.d0/3.d0)*cdat/( abs(cdhfab) )
      qu2faa = 4.0d0*cdhfaa*cdhfaa
      qu2faa1 = 4.0d0*cdhfaa1*cdhfaa1
      qu2fab = 4.0d0*cdhfab*cdhfab

*     Work out quantities needed for averaging wektaor to wekpa, and
*     for computing taux integrals needed for the momentum constraints
*     ----------------------------------------------------------------
      ndxodd = mod(ndxr,2) .eq. 1
*     Compute width of averaging box
      nijwid = ndxr + mod(ndxr,2)
*     Set default weights
      do i=0,ndxr
        wt(i) = 1.0d0
      enddo
*     Modify as appropriate
      if ( ndxodd ) then
        wt(  0 ) = 0.5d0
        wt(ndxr) = 0.5d0
       else
        wt(ndxr) = 0.0d0
      endif

#ifdef temp_fdbck
*============================================================================
*     Added 12/28/2017: Sergey Kravtsov
*      Block of changes for  temperature dependence of the  wind stress
*      the drag coefficient being a function of ocean-atmosphere
*      temperature difference, following Hogg et al. (2009)
*
*     Compute atmospheric mixed-layer temperature at p points

*     Zonal boundaries
*     ----------------
      astp(1, 1)=    0.5d0*( astm(1, 1) + astm(nxta, 1) )
      astp(1, nypa) = 0.5d0*( astm(1, nyta) + astm(nxta, nyta) )
      do i=2,nxpa-1
*       Southern boundary (j = 1)
        astp(i,  1 ) = 0.5d0*( astm(i, 1) + astm(i-1, 1) )
*       Northern boundary (j = nypa)
        astp(i,nypa) = 0.5d0*( astm(i, nyta) + astm(i-1, nyta) )
      enddo
      astp(nxpa, 1) = astp(1, 1)
      astp(nxpa, nypa) = astp(1, nypa)

!$OMP PARALLEL DO DEFAULT (NONE)
!$OMP&         PRIVATE (i,j)
!$OMP&         SHARED  (astm,astp)
!$OMP&         SCHEDULE (STATIC)
*     Inner points + meridional boundaries
*     ------------------------------------
      do j=2,nypa-1
*       Western boundary (periodic)
        astp( 1, j) = 0.25d0*( astm(1, j) + astm(1, j-1)
     &                 + astm(nxta, j) + astm(nxta, j-1))
*       Inner points (standard case)
        do i=2,nxpa-1
           astp( i, j) = 0.25d0*( astm(i, j) + astm(i, j-1)
     &                 + astm(i-1, j) + astm(i-1, j-1))
        enddo
*       Eastern boundary (periodic)
        astp(nxpa, j) = astp(1, j)
      enddo
!$OMP END PARALLEL DO

*     Interpolate atmospheric mixed-layer temperature to ocean resolution
*
*     Using the subroutine for u,v velocities, at the expense of
*     some additional storage
*     ================================================================
      call auvbcu (astp, astpor, astp, astpor1)

*     Compute ocean mixed layer temperature  at p points, and
*     amend atmospheric temperature to O-A temperature diff.
*     zero out the temperature difference field over land
*     ======================================================
*     Uses exactly the same scheme as for the atmosphere,
*     except that it must also handle the box ocean case.
*     Use scalars for ocean geostrophic velocity components
*     as they aren't needed anywhere else in this routine.

*     Zonal boundaries
*     ----------------
#      ifdef cyclic_ocean
      sstp(1, 1)=    0.5d0*( sstm(1, 1) + sstm(nxto, 1) )
      sstp(1, nypo) = 0.5d0*( sstm(1, nyto) + sstm(nxto, nyto) )
#      else
      sstp(1, 1)=     sstm(1, 1)
      sstp(1, nypo) = sstm(1, nyto)
#      endif   
      do i=2,nxpo-1
*       Southern boundary (j = 1)
        sstp(i,  1 ) = 0.5d0*( sstm(i, 1) + sstm(i-1, 1) )
*       Northern boundary (j = nypa)
        sstp(i,nypo) = 0.5d0*( sstm(i, nyto) + sstm(i-1, nyto) )
      enddo
#      ifdef cyclic_ocean
      sstp(nxpo, 1) = sstp(1, 1)
      sstp(nxpo, nypo) = sstp(1, nypo)
#      else
      sstp(nxpo, 1) = sstm(nxto, 1)
      sstp(nxpo, nypo) = sstm(nxto, nyto)
#      endif

!$OMP PARALLEL DO DEFAULT (NONE)
!$OMP&         PRIVATE (i,j)
!$OMP&         SHARED  (sstm,sstp)
!$OMP&         SCHEDULE (STATIC)
*     Inner points + meridional boundaries
*     ------------------------------------
      do j=2,nypo-1
#      ifdef cyclic_ocean
*       Western boundary (periodic)
        sstp( 1, j) = 0.25d0*( sstm(1, j) + sstm(1, j-1)
     &                 + sstm(nxto, j) + sstm(nxto, j-1))
#      else
        sstp( 1, j) = 0.5d0*( sstm(1, j) + sstm(1, j-1))
#      endif
*       Inner points (standard case)
        do i=2,nxpo-1
           sstp( i, j) = 0.25d0*( sstm(i, j) + sstm(i, j-1)
     &                 + sstm(i-1, j) + sstm(i-1, j-1))
        enddo
#      ifdef cyclic_ocean
*       Eastern boundary (periodic)
        sstp(nxpo, j) = sstp(1, j)
#      else
        sstp(nxpo, j) = 0.5d0*( sstm(nxto, j) + sstm(nxto, j-1))
#      endif
      enddo
!$OMP END PARALLEL DO

!$OMP PARALLEL DO DEFAULT (NONE)
!$OMP&         PRIVATE (i,j)
!$OMP&         SHARED  (astpor1)
!$OMP&         SCHEDULE (STATIC)
*     Inner points + meridional boundaries
*     ------------------------------------
      do j=1,nypaor
        do i=1,nxpaor
          astpor1(i,j) = 0.d0
        enddo
      enddo
!$OMP END PARALLEL DO

!$OMP PARALLEL DO DEFAULT (NONE)
!$OMP&         PRIVATE (i,j)
!$OMP&         SHARED  (astpor,astpor1,sstp)
!$OMP&         SCHEDULE (STATIC)
*     Inner points + meridional boundaries
*     ------------------------------------
      do j=1,nypo
        do i=1,nxpo
          astpor1(iocoff+i,jocoff+j) = sstp(i,j)-astpor(iocoff+i,jocoff+j)
        enddo
      enddo
!$OMP END PARALLEL DO

*
*    END OF CHANGES FOR TEMPERATURE-DEPENDENT WIND STRESS FORMULATION
*       (Kravtsov, 12/28/2017)   
*==============================================================================

#endif


*     Compute atmospheric geostrophic velocity at p points
*     ====================================================
*     Bilinearly interpolated fields aren't suitable for further
*     differentiation, so we must derive au1 and av1 from ap1 at
*     the normal atmospheric resolution, and then interpolate
*     them onto the higher (oceanic) resolution grid.

*     Zonal boundaries
*     ----------------
*     Version with nondimensional bccoat
      zbfcat = rdxaf0/( 0.5d0*bccoat + 1.0d0 )
      do i=1,nxpa
*       Southern boundary (j = 1)
        u1at(i,  1 ) = -zbfcat*( pam(i,  2 ,1) - pam(i,   1  ,1) )
        v1at(i,  1 ) = 0.0d0
*       Northern boundary (j = nypa)
        u1at(i,nypa) = -zbfcat*( pam(i,nypa,1) - pam(i,nypa-1,1) )
        v1at(i,nypa) = 0.0d0
      enddo
*
*       corners
*-------------------------------------------
        v1at(nxpa, 1)=v1at(1, 1)
        v1at(nxpa, nypa)=v1at(1, nypa)

!$OMP PARALLEL DO DEFAULT (NONE)
!$OMP&         PRIVATE (i,j)
!$OMP&         SHARED  (hxafac,pam,astm,u1at,v1at)
!$OMP&         SCHEDULE (STATIC)
*     Inner points + meridional boundaries
*     ------------------------------------
      do j=2,nypa-1
*       Western boundary (periodic)
        u1at( 1 ,j) = -hxafac*( pam( 1 ,j+1,1) - pam(   1  ,j-1,1))
        v1at( 1 ,j) =  hxafac*( pam( 2 , j ,1) - pam(nxpa-1, j ,1))
*       Inner points (standard case)
        do i=2,nxpa-1
          u1at(i,j) = -hxafac*( pam( i ,j+1,1) - pam( i ,j-1,1))
          v1at(i,j) =  hxafac*( pam(i+1, j ,1) - pam(i-1, j ,1))
        enddo
*       Eastern boundary (periodic)
        u1at(nxpa,j) = u1at( 1 ,j)
        v1at(nxpa,j) = v1at( 1 ,j)
      enddo
!$OMP END PARALLEL DO
c------------------------------------------------------------------------
c
c New "temperature-dependent" code segment starts here
c  (09/09/2020)
c

#ifdef temp_fdbck
# ifdef partial_coup
*
*      For partial coupling, store temperature-independent velocities
*      in u1at1, v1at1
*
!$OMP PARALLEL DO DEFAULT (NONE)
!$OMP&         PRIVATE (i,j)
!$OMP&         SHARED  (u1at,v1at,u1at1,v1at1)
!$OMP&         SCHEDULE (STATIC)
      do j=1,nypa
        do i=1,nxpa
          u1at1(i,j) = u1at(i,j)
          v1at1(i,j) = v1at(i,j)
        enddo
      enddo
!$OMP END PARALLEL DO
# endif
*
*   Now amend u1at and v1at to include temperature-dependent component
*
**********************************************************************
*
* 02/08/2019: added partially coupled version of temperature dependent
* wind stress, in which the ocean "sees" the wind stress computed by
* neglecting temperature-dependent flow component in the atmospheric
* mixed layer. This is needed to avoid eddy activity suppression in the
* ocean by temperature-dependent wind stress
*
*  In this version, u1at1, v1at1 are the geostrophic component of the atmospheric
*  mixed layer velocities, w/o addition of temperature dependent currents.
*  u1ator1, v1ator1 are these velocities interpolated to ocean resolution.
**********************************************************************

*     Zonal boundaries
*     ----------------
*========================================================================
* dT/dy=0 at N/S boundary, and does not contribute to u component of vel.
* No flow would mean dT/dx=0 there, which can be a problem (?)
* The loop below adds temperature-dependent velocity component
* to geostrophic component in the atmospheric mixed layer.
* Note that this component is (i) non-divergent (so it does not explicitly
* add to Ekman velocity - but does so implicitly via changes to wind stress);
*  and (ii) does not advect atmospheric mixed layer temperature.
* Sergey Kravtsov (24-07-2018)
*========================================================================
*

      do i=2,nxpa-1
*       Southern boundary (j = 1)
        v1at(i,  1 ) = v1at(i,  1 )-hxafac*astfac*
     &  (astm(i,1)-astm(i-1,1)+astm(i,1)-astm(i-1,1))
*       Northern boundary (j = nypa)
        v1at(i,nypa) = v1at(i,nypa)-hxafac*astfac*
     &  (astm(i,nyta)-astm(i-1,nyta)
     &  +astm(i,nyta)-astm(i-1,nyta))
      enddo
*
*       corners
*-------------------------------------------
        v1at(1,  1 ) = v1at(1,  1 )-hxafac*astfac*
     &  (astm(1,1)-astm(nxta,1)+astm(1,1)-astm(nxta,1))
        v1at(1,nypa) = v1at(1,nypa)-hxafac*astfac*
     &  (astm(1,nyta)-astm(nxta,nyta)
     &  +astm(1,nyta)-astm(nxta,nyta))
        v1at(nxpa, 1)=v1at(1, 1)
        v1at(nxpa, nypa)=v1at(1, nypa)

!$OMP PARALLEL DO DEFAULT (NONE)
!$OMP&         PRIVATE (i,j)
!$OMP&         SHARED  (hxafac,astfac,astm,u1at,v1at)
!$OMP&         SCHEDULE (STATIC)
*     Inner points + meridional boundaries
*     ------------------------------------
      do j=2,nypa-1
*       Western boundary (periodic)
        u1at( 1 ,j) = u1at( 1 ,j)+hxafac*
     &  astfac*(astm(1,j)-astm(1,j-1)+astm(nxta,j)-astm(nxta,j-1))
        v1at( 1 ,j) =  v1at(1 ,j)-hxafac*
     &  astfac*(astm(1,j)-astm(nxta,j)+astm(1,j-1)-astm(nxta,j-1))
*       Inner points (standard case)
        do i=2,nxpa-1
          u1at(i,j) = u1at(i,j)+hxafac*
     &  astfac*(astm(i,j)-astm(i,j-1)+astm(i-1,j)-astm(i-1,j-1))
          v1at(i,j) = v1at(i,j)- hxafac*
     &  astfac*(astm(i,j)-astm(i-1,j)+astm(i,j-1)-astm(i-1,j-1))
        enddo
*       Eastern boundary (periodic)
        u1at(nxpa,j) = u1at( 1 ,j)
        v1at(nxpa,j) = v1at( 1 ,j)
      enddo
!$OMP END PARALLEL DO

#endif
c
c End of  "temperature-dependent" segment (09/09/2020)
c--------------------------------------------------------------------------

*     Interpolate atmospheric geostrophic velocity to ocean resolution
*     ================================================================
      call auvbcu (u1at, u1ator, v1at, v1ator)

#ifdef partial_coup
      call auvbcu (u1at1, u1ator1, v1at1, v1ator1)
#endif

*     Check conservation of momentum
*     (area integrals of velocity components)
*     This is only really a useful test for the u-component in the
*     cyclic case; in all other cases the integrals should be zero,
*     so this would check numerical precision rather than the code.
#  ifdef cyclic_ocean
**    call xintp (uvint, u1at, nxpa, nypa)
**    write (*,299) '  xforc: Integ ( u1at ) daA   = ',dxa*dya*uvint
**    call xintp (uvint, u1ator, nxpaor, nypaor)
**    write (*,299) '  xforc: Integ (u1ator) daA   = ',dxo*dyo*uvint
#  endif

#  ifndef atmos_only
*     Logically enough, can't compute the dependence on the
*     oceanic velocity in an atmosphere only configuration!
#    ifdef tau_udiff

*     Compute geostrophic ocean velocity at p points, and
*     amend atmospheric velocity to Q-G layer velocity diff.
*     ======================================================
*     Uses exactly the same scheme as for the atmosphere,
*     except that it must also handle the box ocean case.
*     Use scalars for ocean geostrophic velocity components
*     as they aren't needed anywhere else in this routine.

*     Zonal boundaries
*     ----------------
*     Note that v1ator is unchanged on zonal boundaries because the
*     geostrophic ocean vel. v1oc vanishes identically since p = const
*     Version with nondimensional bccooc
      zbfcoc = rdxof0/( 0.5d0*bccooc + 1.0d0 )
      do i=1,nxpo
*       Southern boundary (j = 1)
        u1oc = -zbfcoc*( pom(i,  2 ,1) - pom(i,   1  ,1) )
        u1ator(iocoff+i,jocoff+1) = u1ator(iocoff+i,jocoff+1) - u1oc
*       Northern boundary (j = nypo)
        u1oc = -zbfcoc*( pom(i,nypo,1) - pom(i,nypo-1,1) )
        u1ator(iocoff+i,jocoff+nypo) = u1ator(iocoff+i,jocoff+nypo)
     &                                - u1oc
      enddo

!$OMP PARALLEL DO DEFAULT (NONE)
!$OMP&         PRIVATE (i,j,u1oc,v1oc)
!$OMP&         SHARED  (hxofac,pom,zbfcoc,u1ator,v1ator)
!$OMP&         SCHEDULE (STATIC)
*     Inner points + meridional boundaries
*     ------------------------------------
      do j=2,nypo-1
*       Western boundary (i = 1)
#      ifdef cyclic_ocean
*       Zonally cyclic ocean
        u1oc = -hxofac*( pom( 1 ,j+1,1) - pom(   1  ,j-1,1) )
        v1oc =  hxofac*( pom( 2 , j ,1) - pom(nxpo-1, j ,1) )
#      else
*       Finite box ocean
        u1oc = 0.0d0
        v1oc = zbfcoc*( pom(  2 ,j,1) - pom(   1  ,j,1) )
#      endif
        u1ator(iocoff+1,jocoff+j) = u1ator(iocoff+1,jocoff+j) - u1oc
        v1ator(iocoff+1,jocoff+j) = v1ator(iocoff+1,jocoff+j) - v1oc

*       Inner points (standard case)
        do i=2,nxpo-1
          u1oc = -hxofac*( pom( i ,j+1,1) - pom( i ,j-1,1) )
          v1oc =  hxofac*( pom(i+1, j ,1) - pom(i-1, j ,1) )
          u1ator(iocoff+i,jocoff+j) = u1ator(iocoff+i,jocoff+j) - u1oc
          v1ator(iocoff+i,jocoff+j) = v1ator(iocoff+i,jocoff+j) - v1oc
        enddo
*       Eastern boundary (i = nxpo)
#      ifdef cyclic_ocean
*       Zonally cyclic ocean
        u1oc = -hxofac*( pom(nxpo,j+1,1) - pom( nxpo ,j-1,1) )
        v1oc =  hxofac*( pom(  2 , j ,1) - pom(nxpo-1, j ,1) )
#      else
*       Finite box ocean
        u1oc = 0.0d0
        v1oc = zbfcoc*( pom(nxpo,j,1) - pom(nxpo-1,j,1) )
#      endif
        u1ator(iocoff+nxpo,jocoff+j) = u1ator(iocoff+nxpo,jocoff+j)
     &                                - u1oc
        v1ator(iocoff+nxpo,jocoff+j) = v1ator(iocoff+nxpo,jocoff+j)
     &                                - v1oc

      enddo
!$OMP END PARALLEL DO

c---------------------------------------------------------------------
c 
c Now the same for the partially coupled case  (09/09/2020)
c

#ifdef partial_coup

*     Zonal boundaries
*     ----------------
*     Note that v1ator is unchanged on zonal boundaries because the
*     geostrophic ocean vel. v1oc vanishes identically since p = const
*     Version with nondimensional bccooc
      zbfcoc = rdxof0/( 0.5d0*bccooc + 1.0d0 )
      do i=1,nxpo
*       Southern boundary (j = 1)
        u1oc = -zbfcoc*( pom(i,  2 ,1) - pom(i,   1  ,1) )
        u1ator1(iocoff+i,jocoff+1) = u1ator1(iocoff+i,jocoff+1) - u1oc
*       Northern boundary (j = nypo)
        u1oc = -zbfcoc*( pom(i,nypo,1) - pom(i,nypo-1,1) )
        u1ator1(iocoff+i,jocoff+nypo) = u1ator1(iocoff+i,jocoff+nypo)
     &                                - u1oc
      enddo

!$OMP PARALLEL DO DEFAULT (NONE)
!$OMP&         PRIVATE (i,j,u1oc,v1oc)
!$OMP&         SHARED  (hxofac,pom,zbfcoc,u1ator1,v1ator1)
!$OMP&         SCHEDULE (STATIC)
*     Inner points + meridional boundaries
*     ------------------------------------
      do j=2,nypo-1
*       Western boundary (i = 1)
#      ifdef cyclic_ocean
*       Zonally cyclic ocean
        u1oc = -hxofac*( pom( 1 ,j+1,1) - pom(   1  ,j-1,1) )
        v1oc =  hxofac*( pom( 2 , j ,1) - pom(nxpo-1, j ,1) )
#      else
*       Finite box ocean
        u1oc = 0.0d0
        v1oc = zbfcoc*( pom(  2 ,j,1) - pom(   1  ,j,1) )
#      endif
        u1ator1(iocoff+1,jocoff+j) = u1ator1(iocoff+1,jocoff+j) - u1oc
        v1ator1(iocoff+1,jocoff+j) = v1ator1(iocoff+1,jocoff+j) - v1oc

*       Inner points (standard case)
        do i=2,nxpo-1
          u1oc = -hxofac*( pom( i ,j+1,1) - pom( i ,j-1,1) )
          v1oc =  hxofac*( pom(i+1, j ,1) - pom(i-1, j ,1) )
          u1ator1(iocoff+i,jocoff+j) = u1ator1(iocoff+i,jocoff+j) - u1oc
          v1ator1(iocoff+i,jocoff+j) = v1ator1(iocoff+i,jocoff+j) - v1oc
        enddo
*       Eastern boundary (i = nxpo)
#      ifdef cyclic_ocean
*       Zonally cyclic ocean
        u1oc = -hxofac*( pom(nxpo,j+1,1) - pom( nxpo ,j-1,1) )
        v1oc =  hxofac*( pom(  2 , j ,1) - pom(nxpo-1, j ,1) )
#      else
*       Finite box ocean
        u1oc = 0.0d0
        v1oc = zbfcoc*( pom(nxpo,j,1) - pom(nxpo-1,j,1) )
#      endif
        u1ator1(iocoff+nxpo,jocoff+j) = u1ator1(iocoff+nxpo,jocoff+j)
     &                                - u1oc
        v1ator1(iocoff+nxpo,jocoff+j) = v1ator1(iocoff+nxpo,jocoff+j)
     &                                - v1oc

      enddo
!$OMP END PARALLEL DO

#endif

#    endif
#  endif

!$OMP PARALLEL DEFAULT (NONE)
!$OMP&         PRIVATE (i,j)
!$OMP&         SHARED  (cdrfaa,cdrfab,qu2faa,qu2faa1,qu2fab,cdat)
!$OMP&         SHARED  (cdhfaa,cdhfaa1,cdhfab)
!$OMP&         SHARED  (u1ator,v1ator,tauxaor,tauyaor,tauxa,tauya,
!$OMP&                 uvekfc,hmrdxa,uekat1,vekat1,wektaor,hxofac,
!$OMP&                 wekta,nijwid,wt,wekpa,wspdor,wspdot)
#ifdef temp_fdbck
!$OMP&        SHARED  (astpor1)
# ifdef partial_coup
!$OMP&        SHARED  (u1ator1,v1ator1,tauxaor1,tauyaor1)
# endif
#endif

*     Compute dynamic windstress from atmospheric geostrophic velocity
*     or Q-G velocity difference as appropriate, using quadratic drag law
*     ===================================================================
*     Atmospheric stress from equations (7.1-7.4):
!$OMP DO SCHEDULE (STATIC)
!$OMP&  PRIVATE (cdrfac,qu2fac,cdhfac,delu1,delv1,
!$OMP&           delt,scasqd,scashr,cdochi)
      do j=1,nypaor
*       Set a & b coefficients as appropriate for (7.1-7.4)
*       Default case (uses just atmospheric velocity,
*       not velocity shear; always applies over land)
        do i=1,nxpaor
          cdrfac(i) = cdrfaa
          qu2fac(i) = qu2faa
          cdhfac(i) = cdhfaa
        enddo
c-----------------------------------------------------------------------       
c       Lower drag over ocean (Sergey Kravtsov, 09/04/2020)
c       
        if ( j.ge.jpobeg .and. j.le.jpoend ) then
          do i=ipobeg,ipoend
            qu2fac(i) = qu2faa1
            cdhfac(i) = cdhfaa1
          enddo
        endif
c-----------------------------------------------------------------------

#  ifdef tau_udiff
*       If using the velocity-shear formulation, modify the
*       coefficients when over the ocean; land areas always
*       use the original dependence on atmospheric velocity only
        if ( j.ge.jpobeg .and. j.le.jpoend ) then
          do i=ipobeg,ipoend
            cdrfac(i) = cdrfab
            qu2fac(i) = qu2fab
            cdhfac(i) = cdhfab
          enddo
        endif
#  endif
*       The i-loops setting the coefficients are separated to
*       ensure that the main computational loop will vectorise
        do i=1,nxpaor
*         Compute magnitude of velocity or velocity difference
*         by solving the quadratic for the square of
*         the shear magnitude "M" (as defined by 7.4)
*         This version solves for scashr = |a+b|*M
          delu1 = u1ator(i,j)
          delv1 = v1ator(i,j)
          delt=1.d0
#ifdef temp_fdbck
          delt = max(1.d0+alp*astpor1(i,j),0.1d0)
#endif
          scasqd = -0.5d0 + 0.5d0*sqrt( 1.0d0 + delt*delt*qu2fac(i)*
     &                            ( delu1*delu1 + delv1*delv1 ) )
          scashr = sqrt( scasqd )
          cdochi = cdrfac(i)*scashr/( 1.0d0 + scasqd )
          tauxaor(i,j) = cdochi*( delu1 - scashr*delv1 )
          tauyaor(i,j) = cdochi*( delv1 + scashr*delu1 )
c----------------------------------------------------------------
c      (Relative) wind speed M, at ocean resolution
c         (Sergey Kravtsov; 09/04/2020)
c
c          wspdor(i,j)=scashr/(cdhfac(i)*delt)
c
c      In fact, it makes sense to include the temperature-difference
c      dependence in the heat-exchange coefficients too! If so,
c      no division on delt is necessary:
c
          wspdor(i,j)=scashr/cdhfac(i)
c----------------------------------------------------------------

*******************************************************************
*
* Now use the same formulas, but assume there are no temperature
* dependent currents, and no drag-coefficient dependence on the  ocean-
* atmosphere temperature difference, in the atmospheric mixed layer.
* The  ocean will "see" the wind stress based on this assumption.
*
* added by S. Kravtsov 02/08/2019
*
******************************************************************

#ifdef partial_coup

          delu1 = u1ator1(i,j)
          delv1 = v1ator1(i,j)
          scasqd = -0.5d0 + 0.5d0*sqrt( 1.0d0 + qu2fac(i)*
     &                            ( delu1*delu1 + delv1*delv1 ) )
          scashr = sqrt( scasqd )
          cdochi = cdrfac(i)*scashr/( 1.0d0 + scasqd )
          tauxaor1(i,j) = cdochi*( delu1 - scashr*delv1 )
          tauyaor1(i,j) = cdochi*( delv1 + scashr*delu1 )
#endif
        enddo
      enddo
!$OMP END DO
c-----------------------------------------------------------------
c
c     Compute relative wind speed at ocean T points, for the use in
c     the  heat-flux bulk formulas later on
c    
c     (Sergey Kravtsov,  09/04/2020)
c

!$OMP DO SCHEDULE (STATIC)
!$OMP&   PRIVATE  (io,jo)
      do jo=jpobeg,jpoend-1
        do io=ipobeg,ipoend-1
          wspdot(io-ipobeg+1,jo-jpobeg+1) = 0.25d0*
     &             ( wspdor(io+1,jo)+wspdor(io+1,jo+1)
     &              + wspdor(io,jo )+wspdor(io,jo+1) )
        enddo
      enddo
!$OMP END DO
c-----------------------------------------------------------------

*
*     Sample tau onto standard resolution atmos. grid
*     ===============================================
*     Copy across values at those points common to the two grids
!$OMP DO SCHEDULE (STATIC)
!$OMP&   PRIVATE  (ia,ja)
      do ja=1,nypa
        do ia=1,nxpa
          tauxa(ia,ja) = tauxaor(1+(ia-1)*ndxr,1+(ja-1)*ndxr)
          tauya(ia,ja) = tauyaor(1+(ia-1)*ndxr,1+(ja-1)*ndxr)
        enddo
      enddo
!$OMP END DO NOWAIT
*
*     Compute horizontal Ekman velocity components at stan-
*     dard resolution, for use in advection terms in amladf
*     =====================================================
*     These use averages of the high-resolution tau compo-
*     nents along the sides of the atmospheric T cells
!$OMP DO SCHEDULE (STATIC)
!$OMP&   PRIVATE  (ia,ja,ioff,joff,tausum)
      do ja=1,nypa
        joff = 1 + (ja-1)*ndxr
        do ia=1,nxta
          ioff = 1 + (ia-1)*ndxr
*         Integrate taux along zonal side of cell
          tausum = 0.5d0*tauxaor(ioff,joff)
          do i=1,ndxr-1
            tausum = tausum + tauxaor(ioff+i,joff)
          enddo
          tausum = tausum + 0.5d0*tauxaor(ioff+ndxr,joff)
          vekat1(ia,ja) = uvekfc*tausum
        enddo
      enddo
!$OMP END DO

*     Also compute wekta from the divergence of uekat, vekat
*     ------------------------------------------------------
*     The array wekta is used in subroutine aml to timestep ast,
*     and is otherwise only used for diagnostics and checking.
!$OMP DO SCHEDULE (STATIC)
!$OMP&   PRIVATE  (ia,ja,ioff,joff,tausum)
      do ja=1,nyta
        joff = 1 + (ja-1)*ndxr
        do ia=1,nxta
          ioff = 1 + (ia-1)*ndxr
*         Integrate tauy along meridional side of cell
          tausum = 0.5d0*tauyaor(ioff,joff)
          do j=1,ndxr-1
            tausum = tausum + tauyaor(ioff,joff+j)
          enddo
          tausum = tausum + 0.5d0*tauyaor(ioff,joff+ndxr)
          uekat1(ia,ja) = -uvekfc*tausum
        enddo
*       Ensure cyclicity
        uekat1(nxpa,ja) = uekat1(1,ja)
*       Compute wekta = -hmat*( d/dx(uekat) + d/dy(vekat) )
        do ia=1,nxta
          wekta(ia,ja) = -hmrdxa*(  uekat1(ia+1,ja) - uekat1(ia,ja)
     &                            + vekat1(ia,ja+1) - vekat1(ia,ja) )
        enddo
      enddo
!$OMP END DO NOWAIT

*     Compute atmospheric Ekman velocity at T points - equation (7.6)
*     ===============================================================
*     Include values of tau-x on the S & N boundaries.
*     If these are to be zero, set them explicitly before use.
*     This removes the need for separate S & N boundary code.
!$OMP DO SCHEDULE (STATIC)
      do j=1,nytaor
        do i=1,nxtaor
          wektaor(i,j) = hxofac*( tauyaor(i+1, j )+tauyaor(i+1,j+1)
     &                         -( tauyaor( i , j )+tauyaor( i ,j+1) )
     &                          + tauxaor( i , j )+tauxaor(i+1, j )
     &                         -( tauxaor( i ,j+1)+tauxaor(i+1,j+1) ) )
        enddo
      enddo
!$OMP END DO

*     Derive atmospheric Ekman velocity at p points
*     ---------------------------------------------
*     The array wekpa is used in subroutine qgastep to timestep qa,
*     and is otherwise only used for diagnostics and checking.
*     Get wekpa by averaging those values of wektaor nearest to
*     the required p-point. This is simple when ndxr is even, but
*     requires some weighting when ndxr is odd, as T-points
*     can then be equidistant from 2 (or even 4) p-points.
*     This approach gives the simplest Stokes' thm equivalent.
!$OMP DO SCHEDULE (STATIC)
!$OMP&   PRIVATE  (ia,ja,jbeg,jlo,jhi,ibeg,it,wsum,wtasum,wtj)
      do ja=1,nypa
        jbeg = (ja-1)*ndxr - (ndxr-1)/2
        jlo = max( 1, jbeg )
        jhi = min( jbeg+nijwid-1, nytaor )
**      write (*,'(a,4i6)') '  xforc: wekpa: ja, jbeg, jlo, jhi = ',
**   &                      ja,jbeg,jlo,jhi
        do ia=1,nxpa
          ibeg = (ia-1)*ndxr - (ndxr-1)/2
          wsum = 0.0d0
          wtasum = 0.0d0
          do j=jlo,jhi
            wtj = wt(j-jbeg)
            do i=ibeg,ibeg+nijwid-1
              it = 1 + mod( i-1+nxtaor, nxtaor )
              wsum = wsum + wt(i-ibeg)*wtj
              wtasum = wtasum + wt(i-ibeg)*wtj*wektaor(it,j)
            enddo
          enddo
          wekpa(ia,ja) = wtasum/wsum
        enddo
*       Check cyclicity
**      if ( wekpa(1,ja).ne.wekpa(nxpa,ja) ) then
**        write (*,'(a,i6,1p,2d20.12)') '  xforc: ja, wekpa = ',
**   &                                  ja,wekpa(1,ja),wekpa(nxpa,ja)
**      endif
      enddo
!$OMP END DO NOWAIT

!$OMP END PARALLEL


#ifndef temp_fdbck

!$OMP PARALLEL DO DEFAULT (NONE)
!$OMP&         PRIVATE (i,j)
!$OMP&         SHARED  (uekat,uekat1,vekat,vekat1)
!$OMP&         SCHEDULE (STATIC)
      do j=1,nypa
        do i=1,nxpa
          uekat(i,j) = uekat1(i,j)
          vekat(i,j) = vekat1(i,j)
        enddo
      enddo
!$OMP END PARALLEL DO

#else

*=======================================================================
*=======================================================================
*     Add temperature-driven component to atmospheric mixed-layer
*     velocities (Kravtsov, 08/28/2018)
*     ====================================================
*     The advection contribution from these velocity components are
*     zero in the temperature advection equation (so uekat1,  vekat1
*     are used to advect temperature there in amlsubs.F, but it is
*     non-trivial in the moisture equation. Their contribution to Ekman
*     pumping (computed just above) is also identically zero.
*

*     Zonal boundaries (no corners)
*     ------------------------------------
*
      do i=2,nxta
*       Southern boundary (j = 1)
        uekat(i,  1 ) = uekat1(i,  1 )+0.5d0*astfac*hxafac*
     &  (astm(i, 2)-astm(i, 1)+astm(i-1, 2)-astm(i-1, 1))
*       Northern boundary (j = nyta)
        uekat(i,nyta) = uekat1(i,nyta)+0.5d0*astfac*hxafac*
     & (astm(i, nyta)-astm(i, nyta-1)
     & +astm(i-1, nyta)-astm(i-1, nyta-1))
      enddo
      do i=2,nxta-1
*       Southern boundary (j = 1)
        vekat(i,  1 ) = vekat1(i,  1 )-astfac*hxafac*
     &  (astm(i+1, 1)-astm(i-1, 1))
*       Northern boundary (j = nypa)
        vekat(i,nypa) = vekat1(i,nypa)-astfac*hxafac*
     & (astm(i+1, nyta)-astm(i-1, nyta))
      enddo
*
*     Corners
*     ------------------------------------
*
        uekat(1,  1 ) = uekat1(1,  1 )+0.5d0*astfac*hxafac*
     &  (astm(1, 2)-astm(1, 1)+astm(nxta, 2)-astm(nxta, 1))
        uekat(nxpa,  1 ) = uekat(1,  1 )
        uekat(1,nyta) = uekat1(1,nyta)+0.5d0*astfac*hxafac*
     & (astm(1, nyta)-astm(1, nyta-1)
     & +astm(nxta, nyta)-astm(nxta, nyta-1))
        uekat(nxpa,  nyta ) = uekat(1,  nyta )
*
        vekat(1,  1 ) = vekat1(1,  1 )-astfac*hxafac*
     &  (astm(2, 1)-astm(nxta, 1))
        vekat(1,nypa) = vekat1(1,nypa)-astfac*hxafac*
     & (astm(2, nyta)-astm(nxta, nyta))
        vekat(nxta,  1 ) = vekat1(nxta,  1 )-astfac*hxafac*
     &  (astm(1, 1)-astm(nxta-1, 1))
        vekat(nxta,nypa) = vekat1(nxta,nypa)-astfac*hxafac*
     & (astm(1, nyta)-astm(nxta-1, nyta))
*
*========================================================================
* v~=0 on zonal boundaies. No flow would mean dT/dx=0 there, is it a
* problem (?)
* Used above are no-flux condition on zonal boundaries and periodic
* conditions in x (for temperature).
*========================================================================
*
*
!$OMP PARALLEL DO DEFAULT (NONE)
!$OMP&         PRIVATE (i,j)
!$OMP&         SHARED  (hxafac,astfac,astm,vekat,vekat1)
!$OMP&         SCHEDULE (STATIC)
*     Inner points + meridional boundaries
*     ------------------------------------
*
      do j=2,nypa-1
*       Western boundary (periodic)
          vekat(1,j) = vekat1(1,j)-0.5d0*astfac*hxafac*
     &   (astm(2,j)-astm(nxta,j)+astm(2,j-1)-astm(nxta,j-1))
*       Inner points (standard case)
        do i=2,nxta-1
          vekat(i,j) = vekat1(i,j)-0.5d0*astfac*hxafac*
     &   (astm(i+1,j)-astm(i-1,j)+astm(i+1,j-1)-astm(i-1,j-1))
        enddo
*       Eastern boundary (periodic)
          vekat(nxta,j) = vekat1(nxta,j)-0.5d0*astfac*hxafac*
     &   (astm(1,j)-astm(nxta-1,j)+astm(1,j-1)-astm(nxta-1,j-1))
      enddo
!$OMP END PARALLEL DO

!$OMP PARALLEL DO DEFAULT (NONE)
!$OMP&         PRIVATE (i,j)
!$OMP&         SHARED  (hxafac,astfac,astm,uekat,uekat1)
!$OMP&         SCHEDULE (STATIC)
*     Inner points + meridional boundaries
*     ------------------------------------
*
      do j=2,nyta-1
*       Western boundary (periodic)
          uekat(1,j) = uekat1(1,j)+0.5d0*astfac*hxafac*
     &   (astm(1,j+1)-astm(1,j-1)+astm(nxta,j+1)-astm(nxta,j-1))
*       Inner points (standard case)
        do i=2,nxpa-1
          uekat(i,j) = uekat1(i,j)+0.5d0*astfac*hxafac*
     &   (astm(i,j+1)-astm(i,j-1)+astm(i-1,j+1)-astm(i-1,j-1))
        enddo
*       Eastern boundary (periodic)
          uekat(nxpa,j) = uekat(1,j)
      enddo
!$OMP END PARALLEL DO

*     End of addition of temperature-driven component to atmospheric
*     mixed-layer
*     velocities (Kravtsov, 08/28/2018)
*
*=================================================================================
*==================================================================================

#endif


c*     Check conservation of atmos. Ekman velocity
c      call xintt (wekint, wektaor, nxtaor, nytaor)
c      write (*,299) '  xforc: Integ (wektaor) daA  = ',dxo*dyo*wekint
c      call xintt (wekint, wekta, nxta, nyta)
c      write (*,299) '  xforc: Integ ( wekta ) daA  = ',dxa*dya*wekint
c      call xintp (wekint, wekpa, nxpa, nypa)
c      write (*,299) '  xforc: Integ ( wekpa ) daA  = ',dxa*dya*wekint
c 

#  ifndef ocean_only
*     Calculate contributions of integrated Ekman velocity at
*     internal p points to the atmos. momentum constraints.
*     -------------------------------------------------------
*     Compute by applying Stokes' theorem 1/2 (atmos. res)
*     gridlength in from the Southern & Northern boundaries,
*     i.e. surrounding those p-cells which are timestepped.
*     This is simple when ndxr is even, but
*     requires some averaging when ndxr is odd.
      if ( ndxodd ) then
        txsums = 0.5d0*( tauxaor(1,jsou) + tauxaor(1,jsou+1) )
        txsumn = 0.5d0*( tauxaor(1,jnor) + tauxaor(1,jnor-1) )
        do i=2,nxpaor-1
          txsums = txsums + ( tauxaor(i,jsou) + tauxaor(i,jsou+1) )
          txsumn = txsumn + ( tauxaor(i,jnor) + tauxaor(i,jnor-1) )
        enddo
        txsums = txsums + 0.5d0*(  tauxaor(nxpaor,jsou)
     &                           + tauxaor(nxpaor,jsou+1) )
        txsumn = txsumn + 0.5d0*(  tauxaor(nxpaor,jnor)
     &                           + tauxaor(nxpaor,jnor-1) )
        txisat = 0.5d0*dxo*txsums
        txinat = 0.5d0*dxo*txsumn
       else
        txsums = 0.5d0*tauxaor(1,jsou)
        txsumn = 0.5d0*tauxaor(1,jnor)
        do i=2,nxpaor-1
          txsums = txsums + tauxaor(i,jsou)
          txsumn = txsumn + tauxaor(i,jnor)
        enddo
        txsums = txsums + 0.5d0*tauxaor(nxpaor,jsou)
        txsumn = txsumn + 0.5d0*tauxaor(nxpaor,jnor)
        txisat = dxo*txsums
        txinat = dxo*txsumn
      endif

!!*     Compute Integ (wekpa) daAI for checking
!!*     ---------------------------------------
!!*     Should agree with Stokes' thm estimate computed above
!!      wekint = 0.0d0
!!!$OMP PARALLEL DEFAULT (NONE)
!!!$OMP&         PRIVATE (i,j,wpsumi)
!!!$OMP&         SHARED  (wekpa,wekint)
!!!$OMP DO SCHEDULE (STATIC)
!!!$OMP&   REDUCTION(+:wekint)
!!*     Inner points + 0.5d0*( W & E boundaries)
!!      do j=2,nypa-1
!!        wpsumi = 0.5d0*wekpa(1,j)
!!        do i=2,nxpa-1
!!          wpsumi = wpsumi + wekpa(i,j)
!!        enddo
!!        wpsumi = wpsumi + 0.5d0*wekpa(nxpa,j)
!!        wekint = wekint + wpsumi
!!      enddo
!!!$OMP END DO NOWAIT
!!!$OMP END PARALLEL
!!      write (*,299) '  xforc: Stokes theorem wekpa = ',
!!     &              (txisat-txinat)/fnot
!!      write (*,299) '  xforc: Integ ( wekpa ) daAI = ',dxa*dya*wekint
#  endif

#  ifndef atmos_only
*     Rescale atmospheric stresses tauxaor, tauyaor
*     to derive oceanic stresses tauxo and tauyo
*     =============================================
*     taux, tauy are actually dynamic stresses, so need to multiply
*     by density ratio across interface - equation (7.5)
*
*
*****************************************************************
* in this version of the code, the ocean can "see" temperature independent
* wind stress (Kravtsov, 02/08/2019);  
* This is governed by the partial_coup  switch 
******************************************************************
*
*
#ifndef partial_coup

!$OMP PARALLEL DO DEFAULT (NONE)
!$OMP&         PRIVATE (io,jo)
!$OMP&         SHARED  (tauxo,tauyo,raoro,tauxaor,tauyaor)
!$OMP&         SCHEDULE (STATIC)
      do jo=1,nypo
        do io=1,nxpo
          tauxo(io,jo) = raoro*tauxaor(iocoff+io,jocoff+jo)
          tauyo(io,jo) = raoro*tauyaor(iocoff+io,jocoff+jo)
        enddo
      enddo
!$OMP END PARALLEL DO

#else

!$OMP PARALLEL DO DEFAULT (NONE)
!$OMP&         PRIVATE (io,jo)
!$OMP&         SHARED  (tauxo,tauyo,raoro,tauxaor1,tauyaor1)
!$OMP&         SCHEDULE (STATIC)
      do jo=1,nypo
        do io=1,nxpo
          tauxo(io,jo) = raoro*tauxaor1(iocoff+io,jocoff+jo)
          tauyo(io,jo) = raoro*tauyaor1(iocoff+io,jocoff+jo)
        enddo
      enddo
!$OMP END PARALLEL DO

#endif



#  endif

#endif

*     In the ocean only case, the code resumes here with
*     tauxo, tauyo having been specified in the main program

#ifndef atmos_only

!$OMP PARALLEL DEFAULT (NONE)
!$OMP&         SHARED  (wekto,hxofac,tauxo,tauyo,wekpo)

*     Compute oceanic Ekman velocity at T points - equation (7.7)
*     ===========================================================
*     Internal latitudes AND S & N boundaries
!$OMP DO SCHEDULE (STATIC)
!$OMP&   PRIVATE  (i,j)
      do j=1,nyto
*       Internal points AND W & E boundaries
        do i=1,nxto
          wekto(i,j) = hxofac*(   tauyo(i+1,j+1)+tauyo(i+1,j)
     &                         -( tauyo( i ,j+1)+tauyo( i ,j) )
     &                          + tauxo(i+1, j )+tauxo(i, j )
     &                         -( tauxo(i+1,j+1)+tauxo(i,j+1) ) )
        enddo
      enddo
!$OMP END DO


c!$OMP PARALLEL DEFAULT (NONE)
c!$OMP&         SHARED  (wekto,hxofac,tauxo,tauyo,wekpo)
c
c*     Compute oceanic Ekman velocity at T points - equation (7.7)
c*     ===========================================================
c*     Internal latitudes AND S & N boundaries
c!$OMP DO SCHEDULE (STATIC)
c!$OMP&   PRIVATE  (i,j)
c
c------------------------------------------------------------------
c
c     ADDED: zero out boundary taus to conserve mass in o.m.l.
c      (Sergey Kravtsov, 07/27/2020)
c
c
c      do j=2,nyto-1
c
c*         W boundary
c          i=1
c#  ifndef cyclic_ocean
c          wekto(i,j) = hxofac*(   tauyo(i+1,j+1)+tauyo(i+1,j)
c     &                          + tauxo(i+1, j )+tauxo(i, j )
c     &                         -( tauxo(i+1,j+1)+tauxo(i,j+1) ) )
c#  else
c          wekto(i,j) = hxofac*(   tauyo(i+1,j+1)+tauyo(i+1,j)
c     &                         -( tauyo( i ,j+1)+tauyo( i ,j) )
c     &                          + tauxo(i+1, j )+tauxo(i, j )
c     &                         -( tauxo(i+1,j+1)+tauxo(i,j+1) ) )
c#  endif
c
c*       Internal points 
c        do i=2,nxto-1
c          wekto(i,j) = hxofac*(   tauyo(i+1,j+1)+tauyo(i+1,j)
c     &                         -( tauyo( i ,j+1)+tauyo( i ,j) )
c     &                          + tauxo(i+1, j )+tauxo(i, j )
c     &                         -( tauxo(i+1,j+1)+tauxo(i,j+1) ) )
c        enddo
c
c#  ifndef cyclic_ocean
c*         E boundary
c          i=nxto
c          wekto(i,j) = hxofac*(   
c     &                         -( tauyo( i ,j+1)+tauyo( i ,j) )
c     &                          + tauxo(i+1, j )+tauxo(i, j )
c     &                         -( tauxo(i+1,j+1)+tauxo(i,j+1) ) )
c#  else
c          wekto(i,j) = hxofac*(   tauyo(i+1,j+1)+tauyo(i+1,j)
c     &                         -( tauyo( i ,j+1)+tauyo( i ,j) )
c     &                          + tauxo(i+1, j )+tauxo(i, j )
c     &                         -( tauxo(i+1,j+1)+tauxo(i,j+1) ) )
c#  endif
c
c      enddo
c
c!$OMP END DO NOWAIT
c
c!$OMP DO SCHEDULE (STATIC)
c
c*       S & N boundaries
c        do i=2,nxto-1
c          j=1
c#  ifndef sb_hflux
c          wekto(i,j) = hxofac*(   tauyo(i+1,j+1)+tauyo(i+1,j)
c     &                         -( tauyo( i ,j+1)+tauyo( i ,j) )
c     &                         -( tauxo(i+1,j+1)+tauxo(i,j+1) ) )
c#  else
c          wekto(i,j) = hxofac*(   tauyo(i+1,j+1)+tauyo(i+1,j)
c     &                         -( tauyo( i ,j+1)+tauyo( i ,j) )
c     &                          + tauxo(i+1, j )+tauxo(i, j )
c     &                         -( tauxo(i+1,j+1)+tauxo(i,j+1) ) )
c#  endif
c          j=nyto
c#  ifndef nb_hflux
c          wekto(i,j) = hxofac*(   tauyo(i+1,j+1)+tauyo(i+1,j)
c     &                         -( tauyo( i ,j+1)+tauyo( i ,j) )
c     &                          + tauxo(i+1, j )+tauxo(i, j )
c     &                                                        )
c#  else
c          wekto(i,j) = hxofac*(   tauyo(i+1,j+1)+tauyo(i+1,j)
c     &                         -( tauyo( i ,j+1)+tauyo( i ,j) )
c     &                          + tauxo(i+1, j )+tauxo(i, j )
c     &                         -( tauxo(i+1,j+1)+tauxo(i,j+1) ) )
c#  endif
c        enddo
c
c!$OMP END DO NOWAIT
c!$OMP END PARALLEL
c
c*       corner points, in memory order:
c*
c       i=1
c       j=1
c
c          wekto(i,j) = hxofac*(   tauyo(i+1,j+1)+tauyo(i+1,j)
c     &                         -( tauxo(i+1,j+1)+tauxo(i,j+1) ) )
c#  ifdef cyclic_ocean
c         wekto(i,j) = wekto(i,j)+hxofac*(
c     &                           -( tauyo( i ,j+1)+tauyo( i ,j) ))
c#  endif
c#  ifdef sb_hflux
c         wekto(i,j) = wekto(i,j)+hxofac*(
c     &                          tauxo(i+1, j )+tauxo(i, j ))
c#  endif
c
c       i=nxto
c       j=1
c          wekto(i,j) = hxofac*(  
c     &                         -( tauyo( i ,j+1)+tauyo( i ,j) )
c     &                         -( tauxo(i+1,j+1)+tauxo(i,j+1) ) )
c#  ifdef cyclic_ocean
c         wekto(i,j) = wekto(i,j)+hxofac*(
c     &                           +( tauyo(i+1,j+1)+tauyo(i+1,j) ))
c#  endif
c#  ifdef sb_hflux
c         wekto(i,j) = wekto(i,j)+hxofac*(
c     &                          tauxo(i+1, j )+tauxo(i, j ))
c#  endif
c
c       i=1
c       j=nyto
c          wekto(i,j) = hxofac*(   tauyo(i+1,j+1)+tauyo(i+1,j)
c     &                          + tauxo(i+1, j )+tauxo(i, j )
c     &                                                         )
c#  ifdef cyclic_ocean
c         wekto(i,j) = wekto(i,j)+hxofac*(
c     &                           -( tauyo( i ,j+1)+tauyo( i ,j) ))
c#  endif
c#  ifdef nb_hflux
c         wekto(i,j) = wekto(i,j)+hxofac*(
c     &                          -( tauxo(i+1,j+1)+tauxo(i,j+1) ))
c#  endif
c
c       i=nxto
c       j=nyto
c          wekto(i,j) = hxofac*(   
c     &                         -( tauyo( i ,j+1)+tauyo( i ,j) )
c     &                          + tauxo(i+1, j )+tauxo(i, j )
c     &                                                        )
c#  ifdef cyclic_ocean
c         wekto(i,j) = wekto(i,j)+hxofac*(
c     &                           +( tauyo(i+1,j+1)+tauyo(i+1,j) ))
c#  endif
c#  ifdef nb_hflux
c         wekto(i,j) = wekto(i,j)+hxofac*(
c     &                          -( tauxo(i+1,j+1)+tauxo(i,j+1) ))
c#  endif
c
c
c
c END OF ADDED CODE
c
c  (Sergey Kravtsov, 07/27/2020)
c
c-----------------------------------------------------------------
c
c
c!$OMP PARALLEL DEFAULT (NONE)
c!$OMP&         SHARED  (wekto,wekpo)

*     Compute oceanic Ekman velocity at p points
*     ------------------------------------------
*     Get wekpo internally by averaging wekto for consistency
!$OMP DO SCHEDULE (STATIC)
!$OMP&   PRIVATE  (io,jo)
      do jo=2,nypo-1
*       Western boundary (io = 1)
#  ifdef cyclic_ocean
*       Zonally cyclic ocean
        wekpo( 1 ,jo) = 0.25d0*(  wekto(nxto,jo-1) + wekto(nxto,jo)
     &                          + wekto(  1 ,jo-1) + wekto(  1 ,jo) )
#  else
*       Finite box ocean
        wekpo( 1 ,jo) = 0.5d0*( wekto(  1 ,jo-1) + wekto(  1 ,jo) )
#  endif
        do io=2,nxpo-1
          wekpo(io,jo) = 0.25d0*(  wekto(io-1,jo-1) + wekto(io-1,jo)
     &                           + wekto( io ,jo-1) + wekto( io ,jo) )
        enddo
*       Eastern boundary (io = nxpo)
#  ifdef cyclic_ocean
*       Zonally cyclic ocean
        wekpo(nxpo,jo) = wekpo( 1 ,jo)
#  else
*       Finite box ocean
        wekpo(nxpo,jo) = 0.5d0*( wekto(nxto,jo-1) + wekto(nxto,jo) )
#  endif
      enddo
!$OMP END DO NOWAIT

!$OMP END PARALLEL

*     Do S & N boundary values to conserve area integral
*     Western boundary (io = 1)
#  ifdef cyclic_ocean
*     Zonally cyclic ocean
      wekpo( 1 ,  1 ) = 0.5d0*( wekto(nxto,  1 ) + wekto( 1 ,  1 ) )
      wekpo( 1 ,nypo) = 0.5d0*( wekto(nxto,nyto) + wekto( 1 ,nyto) )
#  else
*     Finite box ocean
      wekpo( 1 ,  1 ) = wekto(  1 ,  1 )
      wekpo( 1 ,nypo) = wekto(  1 ,nyto)
#  endif
*     Inner points (standard case)
      do io=2,nxpo-1
        wekpo(io,  1 ) = 0.5d0*( wekto(io-1,  1 ) + wekto(io,  1 ) )
        wekpo(io,nypo) = 0.5d0*( wekto(io-1,nyto) + wekto(io,nyto) )
      enddo
*     Eastern boundary (io = nxpo)
#  ifdef cyclic_ocean
*     Zonally cyclic ocean
      wekpo(nxpo,  1 ) = wekpo( 1 ,  1 )
      wekpo(nxpo,nypo) = wekpo( 1 ,nypo)
#  else
*     Finite box ocean
      wekpo(nxpo,  1 ) = wekto(nxto,  1 )
      wekpo(nxpo,nypo) = wekto(nxto,nyto)
#  endif

#  ifdef cyclic_ocean
*     Check cyclicity of wekpo
*     ------------------------
**    print *,' xforc: wekpo cyclicity'
**    do jo=nypo,1,-1
**      if ( wekpo(1,jo).ne.wekpo(nxpo,jo) ) then
**        write (*,'(i6,1p,2d21.12)') jo,wekpo(1,jo),wekpo(nxpo,jo)
**      endif
**    enddo
#  endif

*     Check conservation of ocean Ekman velocity
c      call xintt (wekint, wekto, nxto, nyto)
c      write (*,299) '  xforc: Integ ( wekto ) doA  = ',dxo*dyo*wekint
c      call xintp (wekint, wekpo, nxpo, nypo)
c      write (*,299) '  xforc: Integ ( wekpo ) doA  = ',dxo*dyo*wekint

#  ifdef cyclic_ocean
*     Calculate contributions of integrated Ekman velocity
*     at internal p points to the momentum constraints.
*     ----------------------------------------------------
*     Compute by applying Stokes' theorem 1/2 (atmos. res)
*     gridlength in from the Southern & Northern boundaries,
*     i.e. surrounding those p-cells which are timestepped.
      txsums = 0.5d0*( tauxo(1,   1  ) + tauxo(1,  2 ) )
      txsumn = 0.5d0*( tauxo(1,nypo-1) + tauxo(1,nypo) )
      do io=2,nxpo-1
        txsums = txsums + ( tauxo(io,   1  ) + tauxo(io,  2 ) )
        txsumn = txsumn + ( tauxo(io,nypo-1) + tauxo(io,nypo) )
      enddo
      txsums = txsums + 0.5d0*( tauxo(nxpo,   1  ) + tauxo(nxpo,  2 ) )
      txsumn = txsumn + 0.5d0*( tauxo(nxpo,nypo-1) + tauxo(nxpo,nypo) )
      txisoc = 0.5d0*dxo*txsums
      txinoc = 0.5d0*dxo*txsumn

*     Compute Integ (wekpo) doAI for checking
*     ---------------------------------------
!!      wekint = 0.0d0
!!!$OMP PARALLEL DEFAULT (NONE)
!!!$OMP&         PRIVATE (io,jo,wpsumi)
!!!$OMP&         SHARED  (wekpo,wekint)
!!!$OMP DO SCHEDULE (STATIC)
!!!$OMP&   REDUCTION(+:wekint)
!!*     Inner points + 0.5d0*( W & E boundaries)
!!      do jo=2,nypo-1
!!        wpsumi = 0.5d0*wekpo(1,jo)
!!        do io=2,nxpo-1
!!          wpsumi = wpsumi + wekpo(io,jo)
!!        enddo
!!        wpsumi = wpsumi + 0.5d0*wekpo(nxpo,jo)
!!        wekint = wekint + wpsumi
!!      enddo
!!!$OMP END DO NOWAIT
!!!$OMP END PARALLEL
!!      write (*,299) '  xforc: Stokes theorem wekpo = ',
!!     &              (txisoc-txinoc)/fnot
!!      write (*,299) '  xforc: Integ ( wekpo ) doAI = ',dxo*dyo*wekint
#  endif

#endif

#if defined tau_udiff || !defined ocean_only

*     Compute atmospheric forcing
*     ===========================

*     Interpolate lagged atmospheric temperature onto oceanic grid
*     ------------------------------------------------------------
      call bilint (xta, yta, nxta, nyta, astm,
     &             xto, yto, nxto, nyto, asto, 1.0d0)

#ifdef moist

      call bilint (xta, yta, nxta, nyta, hmixam,
     &             xto, yto, nxto, nyto, hmixao, 1.0d0)

#endif


*     Specify atmospheric forcing everywhere
*     --------------------------------------
*     Land case of (7.9) + last term of (7.8)
*     These values will only be retained over land
*     MONITORING - extra section for arlaav
      arlasm = 0.0d0
!$OMP PARALLEL DEFAULT (NONE)
!$OMP&         PRIVATE (ia,ja,fsp)
!$OMP&         SHARED  (ytarel,fnetat,Dmup,Dmdown,zm,arlasm,astm)

!$OMP DO SCHEDULE (STATIC)
!$OMP&   REDUCTION(+:arlasm)
      do ja=1,nyta
*       Compute radiative forcing perturbation
        fsp = fsprim( ytarel(ja) )
        do ia=1,nxta
*         Atmospheric radiation based on astm
          fnetat(ia,ja) = -zm*fsp - (Dmup-(1.d0-zm)*Dmdown)*astm(ia,ja)
*         MONITORING - extra section for arlaav
          arlasm = arlasm + astm(ia,ja)
        enddo
      enddo
!$OMP END DO NOWAIT

!$OMP END PARALLEL

*     Reset atmospheric forcing to zero over ocean
*     Do the same with ocean temperature estimated
*     at atmo grid points
*     --------------------------------------------
      natocn = 0
      do ja=ny1,ny1+nyaooc-1
        do ia=nx1,nx1+nxaooc-1
          fnetat(ia,ja) = 0.0d0
          ssta(ia-nx1+1,ja-ny1+1)=0.d0
*     Do the same, in the moist model,
*     with the evaporation at atmo grid points

#  ifdef moist
          evapa(ia,ja)=0.d0
#  endif

          
*         MONITORING - extra section for arlaav
*         Count no. of atmos. cells over ocean
          arlasm = arlasm - astm(ia,ja)
          natocn = natocn + 1
        enddo
      enddo

*     MONITORING - extra section for arlaav
      natlan = nxta*nyta - natocn
      if (natlan.eq.0) then
        arlaav = 0.0d0
      else
        arlaav = Dmup*arlasm/dble(natlan)
      endif

*     Specify net ocean forcing (-ve of equation (7.10)),
*     and modify the net atmospheric forcing over ocean
*     points, equation (7.9) + last term of equation (7.8)
*     ----------------------------------------------------
*     fnetoc = - oF0
*     fnetat = - aFm + aF0  (not including eta terms yet)
      ocfrac = dxo*dyo/(dxa*dya)
      fmafac = Adown(1,1)*0.25d0/gpat(1)
      fmafac1 = Adown(1,2)*0.25d0/gpat(2)
      fmatop = 0.25d0*( Cmup + C1down )
      hmafac = -hmadmp - Bmup - B1down
*     MONITORING - extra section for slhfav, oradav, arocav
      slhfsm = 0.0d0
      oradsm = 0.0d0
      arocsm = 0.0d0
# ifdef moist
      evapm = 0.0d0
# endif

!$OMP PARALLEL DEFAULT (NONE)
!$OMP&       PRIVATE (i,j)
!$OMP&       SHARED  (ytorel,D0up,sstm,xlamda,asto,ssta,Dmdown,fnetoc,
!$OMP&              arocsm,Dmup,fnetat,ocfrac,slhfsm,oradsm,zm,radeta,
!$OMP&              wspdot)
!$OMP&   SHARED (fmafac,fmafac1,pam,fmatop,dtopat,hmafac,hmixam,hmat)
#  ifdef moist
!$OMP&   SHARED (hmixao,evap,evapa,rhom,Lv,cpat,Ch,Ce,
!$OMP&           tmbaro,tmbara,gamma,evapm,Flambar)
#  endif

!$OMP DO SCHEDULE (STATIC)
!$OMP&   PRIVATE (io,jo,fsp,ia,ja,ocnrad,slhf,atmrad)
#  ifdef moist
!$OMP&   PRIVATE (sstr,hsat,sstr1,hsat1,hsat2,fact)
!$OMP&   REDUCTION(+:evapm)
#  endif
!$OMP&   REDUCTION(+:arocsm,slhfsm,oradsm)
      do jo=1,nyto
        ja = ny1 + (jo-1)/ndxr
*       Compute radiative forcing perturbation
        fsp = fsprim( ytorel(jo) )
        do io=1,nxto
          ia = nx1 + (io-1)/ndxr

*         Ocean infrared radiation
          ocnrad = D0up*sstm(io,jo)

*         Sensible and latent flux
c-----------------------------------------------------------------------
c         Added wind-speed dependence
c         (Sergey Kravtsov, 09/04/2020)
c
c
#  ifndef moist
          slhf = (xlamda+4.d0*wspdot(io,jo))*(sstm(io,jo)-asto(io,jo))
#  else
          sstr=sstm(io,jo)+tmbaro         
          hsat=eeps*esat(sstr)/1000.d0          
          sstr1=asto(io,jo)+tmbara-0.5d0*(gamma-Gc)*hmat
          sstr=asto(io,jo)+tmbara-0.5d0*gamma*hmat
          fact=(sstr1/sstr)**Gtilde !pressure factor
          hsat1=eeps*esat(sstr1)/1000.d0
          hsat2=eeps*fact*esat(sstr)/1000.d0

          slhf = (xlamda+rhom*cpat*Ch*wspdot(io,jo))
     &                       *(sstm(io,jo)-asto(io,jo))
        evap(io,jo)= max(Ce*wspdot(io,jo)
     & *(hsat-hmixao(io,jo)*hsat1/hsat2),
     & 0.d0)
#  endif
c----------------------------------------------------------------------

#  ifndef atmos_only
*         Atmospheric mixed layer radiation - into ocean
          atmrad = Dmdown*asto(io,jo)

*         Oceanic mixed layer diabatic forcing
          fnetoc(io,jo) = -fsp - atmrad - ocnrad - slhf
# ifdef moist 
          fnetoc(io,jo) = fnetoc(io,jo)-evap(io,jo)*rhom*Lv
c     &                  +Flambar !cancel the mean-state heat loss
# endif
          arocsm = arocsm + atmrad
#  endif

*         Atmospheric mixed layer radiation - in atmosphere
          atmrad = (Dmdown - Dmup)*asto(io,jo)
          fnetat(ia,ja) = fnetat(ia,ja)
     &                 + ocfrac*( zm*ocnrad + atmrad + slhf )

*     Compute SST averages within the boundaries of atmospheric cells.
*     These are used in amlsubs to compute the SST contribution to 
*     atmospheric entrainments (needed when zm<1)
*
*     S. Kravtsov, 08/15/2020
*----------------------------------------------------------------
          ssta(ia-nx1+1,ja-ny1+1) = ssta(ia-nx1+1,ja-ny1+1)
     &                   + ocfrac*sstm(io,jo)
*
*
*    Also compute evaporation on the atmo. grid to step moisture
*    equation in amlsubs
*
*    S. Kravtsov 10/08/2020
*
*----------------------------------------------------------------
#  ifdef moist
          evapa(ia,ja) = evapa(ia,ja)
     &                   + ocfrac*evap(io,jo)
          evapm=evapm+evap(io,jo)
#  endif

*         MONITORING - extra section for slhfav, oradav
          slhfsm = slhfsm + slhf
          oradsm = oradsm + ocnrad
        enddo
      enddo
!$OMP END DO

*     Complete atmospheric diabatic forcing at T points by
*     adding terms involving local thickness perturbations
*     and atmospheric topography. The eta1 and topography
*     terms are averaged from the p grid to the T grid.
*     -----------------------------------------------------
*     These are the 1st, 2nd & 3rd terms of equation (7.8)
!$OMP DO SCHEDULE (STATIC)
      do j=1,nyta
        do i=1,nxta
          radeta(i,j) = - zm*fmafac*(  pam( i , j ,1) - pam( i , j ,2)
     &                           + pam(i+1, j ,1) - pam(i+1, j ,2)
     &                           + pam( i ,j+1,1) - pam( i ,j+1,2)
     &                           + pam(i+1,j+1,1) - pam(i+1,j+1,2) )
     &                 - zm*fmafac1*(  pam( i , j ,2) - pam( i , j ,3)
     &                           + pam(i+1, j ,2) - pam(i+1, j ,3)
     &                           + pam( i ,j+1,2) - pam( i ,j+1,3)
     &                           + pam(i+1,j+1,2) - pam(i+1,j+1,3) )
          fnetat(i,j) =  fnetat(i,j) +  (2.d0-zm)*radeta(i,j)
     &                 - fmatop*(  dtopat( i , j ) + dtopat(i+1, j )
     &                           + dtopat( i ,j+1) + dtopat(i+1,j+1) )
c     &                 + hmafac*( hmixam( i ,j) - hmat )
c
c The last line is commented out, as we either have a dry model with constant
c hmix, or a  moist model, in which hmixa is in fact a specific humidity in the 
c atmospheric mixed layer
c
c  Sergey Kravtsov (10.08.2020)
c
        enddo
      enddo
!$OMP END DO NOWAIT
*----------------------------------------------------------------
*    Modify atmospheric forcing over ocean for the case zm<1:
*
*    (Sergey Kravtsov, 08/15/2020)
*
!$OMP DO SCHEDULE (STATIC)
      do j=ny1,ny1+nyaooc-1
        do i=nx1,nx1+nxaooc-1
          fnetat(i,j) =  fnetat(i,j) + (zm-1.d0)*radeta(i,j)
        enddo
      enddo
!$OMP END DO NOWAIT

!$OMP END PARALLEL

# ifdef moist
      evapm=evapm*ocnorm
# endif
        

*----------------------------------------------------------------
*
* Now also modify the o.m.l. forcing for the case zm<1:
*
*    (Sergey Kravtsov, 08/15/2020)
*
*     Interpolate the appropriate part of the atmo. forcing 
*     onto oceanic grid:
*     ------------------------------------------------------------
      call bilint (xta, yta, nxta, nyta, radeta,
     &             xto, yto, nxto, nyto, radeto, 1.0d0)

!$OMP PARALLEL DEFAULT (NONE)
!$OMP&         PRIVATE (io,jo)
!$OMP&         SHARED  (fnetoc,radeto,zm)
# ifdef moist
!$OMP&         SHARED  (evapm,Lv,rhom)
# endif

!$OMP DO SCHEDULE (STATIC)
      do jo=1,nyto
*       Compute radiative forcing of o.m.l.
*       from atmospheric interior
        do io=1,nxto
         fnetoc(io,jo) = fnetoc(io,jo)
     &                  +(1.d0-zm)*radeto(io,jo)/zm
*
* for the moist model, amend the latent heat exchange to be
* an anomaly wrt mean
*
# ifdef moist
          fnetoc(io,jo) = fnetoc(io,jo)
     &                  +Lv*rhom*evapm
# endif
        enddo
      enddo
!$OMP END DO NOWAIT

!$OMP END PARALLEL


*     MONITORING - extra section for slhfav, oradav, arocav
      slhfav = slhfsm*ocnorm
      oradav = oradsm*ocnorm
      arocav = arocsm*ocnorm
#endif

*     Format for debugging printouts
c  299 format(a,1p,2d22.13)

      END SUBROUTINE xforc
c
c***********************************************************************
c
      DOUBLE PRECISION FUNCTION fsprim (yrel)
*
*     Computes perturbative radiation forcing of system.
*     fspco is peak-to-trough amplitude of perturbative forcing.
*     fspco is > 0 for Northern hemisphere, < 0 for Southern.
*     yrel is y relative to central latitude
*     Function should be chosen to have zero integral
*     over y = [0, yla], i.e. over yrel = [-yla/2, yla/2].
*     Fs is +ve upwards, so function should be -ve near
*     the equator, and become positive poleward.

*     Modules
      USE atconst, ONLY : yla
      USE radiate, ONLY : fspco

      IMPLICIT NONE

*     Subroutine arguments
      double precision yrel
*
      DOUBLE PRECISION PI
      PARAMETER ( PI=3.14159265358979324D0 )

      fsprim = fspco*0.5d0*sin( PI*yrel/yla )

      END FUNCTION fsprim
c
c***********************************************************************
c
      DOUBLE PRECISION FUNCTION esat (tabs)
*
*     Computes saturation water-vapor pressure using Bolton (1980)
*     formula. Only used, obviously, in the moist version of Q-GCM
*
*     tabs is the temperature in Kelvin, trel is the temperature in
*     degrees C. The output is in mb.
*
*     (added by Sergey Kravtsov, 10/08/2020)
*
      IMPLICIT NONE

*     Subroutine arguments
      double precision tabs,trel
*
      trel=tabs-273.15d0

      esat = 6.112d0*exp(17.67d0*trel/(trel+243.5d0))

      END FUNCTION esat
c
c***********************************************************************
c
      SUBROUTINE bilint (xa, ya, nxat, nyat, atmos,
     &                   xo, yo, nxoc, nyoc, ocean, fmult)
*
*     Performs bilinear interpolation of atmos(nxat,nyat), which
*     is tabulated at coordinates xa(nxat), ya(nyat), to fill
*     the array ocean(nxoc,nyoc), tabulated at xo(nxoc), yo(nyoc)
*     Also multiplies the interpolant by the factor fmult
*     Used to transfer data from atmospheric to oceanic grids

*     Modules
      USE parameters
      USE atconst

      IMPLICIT NONE

*     Subroutine arguments
      integer, INTENT(IN) :: nxat,nyat,nxoc,nyoc
      double precision, INTENT(IN) :: xa(nxat),ya(nyat),fmult,
     &                                xo(nxoc),yo(nyoc),atmos(nxat,nyat)
      double precision, INTENT(OUT) :: ocean(nxoc,nyoc)
*
*     Local variables
      integer io,jo,iam(nxpo),iap(nxpo),jam,jap
      double precision dxainv,dyainv,xam,wpx(nxpo),wmx(nxpo),wpy,wmy

      dxainv = 1.0d0/dxa
      dyainv = 1.0d0/dya

*     Get i-subscripts of ocean points in atmos array
*     -----------------------------------------------
*     Compute subscripts once only; store in vector.
*     Assumes fixed grid interval of dxa between xa values.
      do io=1,nxoc
        iam(io) = int( 1.0d0 + dxainv*( xo(io) - xa(1) ) )
        iap(io) = iam(io) + 1
        if ( iam(io).ge.1 ) then
          xam = xa(iam(io))
         else
          xam = xa(1) - dxa
        endif
*       Compute x-direction weights (assumes regular grid)
        wpx(io) = dxainv*( xo(io) - xam )
        wmx(io) = 1.0d0 - wpx(io)
*       Mend both pointers to give correct cyclic results for T points.
*       Results will be inaccurate for p points, but this won't matter
*       because the weight of the inaccurate value will be zero
*       (p points never involve extrapolation, but T points can).
        iam(io) = 1 + mod( iam(io)+nxat-1, nxat )
        iap(io) = 1 + mod( iap(io)+nxat-1, nxat )
**      if ( iam(io).le.0 .or. iam(io).ge.nxat ) then
**        print *,' Error in bilint; atmos. data ',
**   &            'does not cover ocean domain'
**        print *,' Atmos longitude range: ',xa(1),xa(nxat)
**        print *,' Ocean longitude range: ',xo(1),xo(nxoc)
**        print *,' Ocean gridpoint, iam: ',io,xo(io),iam(io)
**        print *,' Program terminates in bilint'
**        stop
**      endif
**      write (*,'(a,3i6,2f14.6)')
**   &         '  bilint: io, iam, iap, wmx, wpx = ',
**   &         io,iam(io),iap(io),wmx(io),wpx(io)
      enddo

*     Compute y-direction weights and perform interpolation
*     -----------------------------------------------------
*     Assumes fixed grid intervals.
!$OMP PARALLEL DO DEFAULT (NONE)
!$OMP&         PRIVATE (jo,jam,jap,wpy,wmy,io)
!$OMP&         SHARED  (nyoc,nyat,nxoc,dyainv,yo,ya,
!$OMP&                  iam,iap,wmx,wpx,atmos,ocean,fmult)
!$OMP&         SCHEDULE (STATIC)
      do jo=1,nyoc
        jam = int( 1.0d0 + dyainv*( yo(jo) - ya(1) ) )
        jap = jam + 1
*       Fix values for extrapolation.
*       Boundary condition is no normal derivative.
        jam = max(jam,  1 )
        jap = min(jap,nyat)
**      if ( jam.le.0 .or. jam.gt.nyat ) then
**        print *,' Error in bilint; atmos. data ',
**   &            'does not cover ocean domain'
**        print *,' Atmos latitude range: ',ya(1),ya(nyat)
**        print *,' Ocean latitude range: ',yo(1),yo(nyoc)
**        print *,' Ocean gridpoint, jam: ',jo,yo(jo),jam
**        print *,' Program terminates in bilint'
**        stop
**      endif
*       Compute y-direction weights (assumes regular grid)
        wpy = dyainv*( yo(jo) - ya(jam) )
        wmy = 1.0d0 - wpy
**      write (*,'(a,3i6,2f14.6)')
**   &         '  bilint: jo, jam, jap, wmy, wpy = ',
**   &         jo,jam,jap,wmy,wpy
        do io=1,nxoc
          ocean(io,jo) = fmult*(  wmx(io)*wmy*atmos(iam(io),jam)
     &                          + wpx(io)*wmy*atmos(iap(io),jam)
     &                          + wmx(io)*wpy*atmos(iam(io),jap)
     &                          + wpx(io)*wpy*atmos(iap(io),jap) )
        enddo
      enddo
!$OMP END PARALLEL DO

      END SUBROUTINE bilint
c
c***********************************************************************
c
      SUBROUTINE auvbcu (u1astd, u1afin, v1astd, v1afin)
*
*     Performs bicubic interpolation of standard resolution
*     atmospheric QG velocity component arrays u1astd, v1astd
*     (tabulated at the standard atmospheric p-points) onto
*     the ocean-resolution atmospheric p-grid required for
*     the velocity difference formulation for the windstress.
*     Results are returned in the arrays u1afin, v1afin.
*     Bicubic interpolation is required to give continuity
*     of derivatives across coarse resolution cell boundaries,
*     since the interpolated u and v fields are (in effect)
*     differentiated when windstress curl is computed.

*     Modules
      USE parameters, ONLY : nxta, nyta, nxpa, nypa, nxpaor, nypaor

      IMPLICIT NONE

*     Subroutine arguments
      double precision, INTENT(IN) :: u1astd(nxpa,nypa),
     &                                v1astd(nxpa,nypa)
      double precision, INTENT(OUT) :: u1afin(nxpaor,nypaor),
     &                                 v1afin(nxpaor,nypaor)
*
*     Local variables
      logical first
      integer ic,jc,icm1,icp2,ifoff,jfoff,jd,k,ii,jj
      double precision udat(16),vdat(16),usum,vsum
*
      SAVE first
*
      DATA first /.true./
*
*     On first pass, initialise matrices &
*     coefficients in the module's storage
*     ------------------------------------
      if ( first ) then
        call bcuini
        first = .false.
      endif
*
*     Interpolate in cells adjacent to southern boundary
*     --------------------------------------------------
      jc = 1
      jfoff = 1 + (jc-1)*ndxr
      do ic=1,nxta
        ifoff = 1 + (ic-1)*ndxr
*       Ensure cyclicity of subscripts
        icm1 = 1 + mod( ic-2+nxta, nxta )
        icp2 = 1 + mod( ic+1, nxta )
**      write (*,'(a,5i6)') ' ic, subscripts: ',ic,icm1,ic,ic+1,icp2
*       Construct u and v data vectors
*       Pad missing u values for jd = -1 with zeros for safety
        do k=1,4
          udat(k) = 0.0d0
        enddo
*       Pad missing v values for jd = -1 with u values on bdy
        vdat(1) = u1astd(icm1,1)
        vdat(2) = u1astd( ic ,1)
        vdat(3) = u1astd(ic+1,1)
        vdat(4) = u1astd(icp2,1)
        k = 4
        do jd=0,2
          udat(k+1) = u1astd(icm1,jc+jd)
          udat(k+2) = u1astd( ic ,jc+jd)
          udat(k+3) = u1astd(ic+1,jc+jd)
          udat(k+4) = u1astd(icp2,jc+jd)
          vdat(k+1) = v1astd(icm1,jc+jd)
          vdat(k+2) = v1astd( ic ,jc+jd)
          vdat(k+3) = v1astd(ic+1,jc+jd)
          vdat(k+4) = v1astd(icp2,jc+jd)
          k = k + 4
        enddo
*       Evaluate the bicubic interpolant on the finer resolution grid
        do jj=0,ndxr-1
          do ii=0,ndxr-1
            usum = 0.0d0
            vsum = 0.0d0
            do k=1,16
              usum = usum + udat(k)*stbus(k,ii,jj)
              vsum = vsum + vdat(k)*stbvs(k,ii,jj)
            enddo
            u1afin(ifoff+ii,jfoff+jj) = usum
            v1afin(ifoff+ii,jfoff+jj) = vsum
          enddo
        enddo
      enddo
*
*     Interpolate in cells within the channel
*     ---------------------------------------
*     For each coarse cell within which we interpolate, the
*     input data points are numbered within the vector thus:
*
*                              jd
*
*       13   14   15   16       2
*
*        9   10---11   12       1
*             |    |
*             |    |
*        5    6----7    8       0
*
*        1    2    3    4      -1
*
*
*     ic, jc are the T-grid cell subscripts on the coarse grid
*     (equivalently the subscripts of the p-gridpoint at the
*      bottom LH corner of the cell within which we interpolate)
*     ifoff, jfoff are the subscripts of the bottom LH corner
*     point of the coarse gridcell, in the fine resolution array
!$OMP PARALLEL DO DEFAULT (NONE)
!$OMP&         PRIVATE (ic,jc,ifoff,jfoff,icm1,icp2,k,
!$OMP&                  jd,udat,vdat,jj,ii,usum,vsum)
!$OMP&         SHARED  (u1astd,v1astd,stbbb,u1afin,v1afin)
!$OMP&         SCHEDULE (STATIC)
      do jc=2,nyta-1
        jfoff = 1 + (jc-1)*ndxr
        do ic=1,nxta
          ifoff = 1 + (ic-1)*ndxr
*         Ensure cyclicity of i-subscripts in coarse array
          icm1 = 1 + mod( ic-2+nxta, nxta )
          icp2 = 1 + mod( ic+1, nxta )
**        write (*,'(a,5i6)') ' ic, subscripts: ',ic,icm1,ic,ic+1,icp2
*         Construct u and v data vectors from the 16 data points
*         surrounding the coarse cell in which we are interpolating
          k = 0
          do jd=-1,2
            udat(k+1) = u1astd(icm1,jc+jd)
            udat(k+2) = u1astd( ic ,jc+jd)
            udat(k+3) = u1astd(ic+1,jc+jd)
            udat(k+4) = u1astd(icp2,jc+jd)
            vdat(k+1) = v1astd(icm1,jc+jd)
            vdat(k+2) = v1astd( ic ,jc+jd)
            vdat(k+3) = v1astd(ic+1,jc+jd)
            vdat(k+4) = v1astd(icp2,jc+jd)
            k = k + 4
          enddo
*         Evaluate the bicubic interpolant on the finer resolution grid
*         ii, jj are the subscripts (relative to ifoff, jfoff) of
*         the fine resolution gridpoints within, or on the boundary
*         of, the coarse gridcell currently being interpolated
*         The array of weights stbbb (a vector of 16 coefficients at
*         each of the fine resolution points) is precomputed in bcuini
          do jj=0,ndxr-1
            do ii=0,ndxr-1
              usum = 0.0d0
              vsum = 0.0d0
              do k=1,16
                usum = usum + udat(k)*stbbb(k,ii,jj)
                vsum = vsum + vdat(k)*stbbb(k,ii,jj)
              enddo
              u1afin(ifoff+ii,jfoff+jj) = usum
              v1afin(ifoff+ii,jfoff+jj) = vsum
            enddo
          enddo
        enddo
      enddo
!$OMP END PARALLEL DO
*
*     Interpolate in cells adjacent to northern boundary
*     --------------------------------------------------
      jc = nyta
      jfoff = 1 + (jc-1)*ndxr
      do ic=1,nxta
        ifoff = 1 + (ic-1)*ndxr
*       Ensure cyclicity of subscripts
        icm1 = 1 + mod( ic-2+nxta, nxta )
        icp2 = 1 + mod( ic+1, nxta )
**      write (*,'(a,5i6)') ' ic, subscripts: ',ic,icm1,ic,ic+1,icp2
*       Construct u and v data vectors
        k = 0
        do jd=-1,1
          udat(k+1) = u1astd(icm1,jc+jd)
          udat(k+2) = u1astd( ic ,jc+jd)
          udat(k+3) = u1astd(ic+1,jc+jd)
          udat(k+4) = u1astd(icp2,jc+jd)
          vdat(k+1) = v1astd(icm1,jc+jd)
          vdat(k+2) = v1astd( ic ,jc+jd)
          vdat(k+3) = v1astd(ic+1,jc+jd)
          vdat(k+4) = v1astd(icp2,jc+jd)
          k = k + 4
        enddo
*       Pad missing u values for jd = 2 with zeros for safety
        do k=13,16
          udat(k) = 0.0d0
        enddo
*       Pad missing v values for jd = 2 with u values on bdy
        vdat(13) = u1astd(icm1,nypa)
        vdat(14) = u1astd( ic ,nypa)
        vdat(15) = u1astd(ic+1,nypa)
        vdat(16) = u1astd(icp2,nypa)
*       Evaluate the bicubic interpolant on the finer resolution grid
        do jj=0,ndxr
          do ii=0,ndxr-1
            usum = 0.0d0
            vsum = 0.0d0
            do k=1,16
              usum = usum + udat(k)*stbun(k,ii,jj)
              vsum = vsum + vdat(k)*stbvn(k,ii,jj)
            enddo
            u1afin(ifoff+ii,jfoff+jj) = usum
            v1afin(ifoff+ii,jfoff+jj) = vsum
          enddo
        enddo
      enddo
*
*     Add cyclicity to fine res. arrays
*     ---------------------------------
      do jj=1,nypaor
        u1afin(nxpaor,jj) = u1afin(1,jj)
        v1afin(nxpaor,jj) = v1afin(1,jj)
      enddo
*
*     Optional code for checking correctness of interpolation
*     -------------------------------------------------------
**    print *,' '
**    print *,' v1astd:'
**    do jfr=nypaor,nypaor-ndxr/2,-1
**      write(*,'(i6,1p,11d14.6)') jfr,(v1astd(ifr,jfr),ifr=1,ndxr/2+1)
**    enddo
**    do jfr=ndxr/2+1,1,-1
**      write(*,'(i6,1p,11d14.6)') jfr,(v1astd(ifr,jfr),ifr=1,ndxr/2+1)
**    enddo
**    print *,' '
**    print *,' v1afin:'
**    do jfr=3*ndxr+1,1,-1
**      write(*,'(i6,1p,11d14.6)') jfr,(v1afin(ifr,jfr),ifr=1,ndxr+1)
**    enddo
**    do jfr=nypaor,nypaor-3*ndxr,-1
**      write(*,'(i6,1p,11d14.6)') jfr,(v1afin(ifr,jfr),ifr=1,ndxr+1)
**    enddo
**    print *,' '
**    print *,' check v1afin cyclicity:'
**    do jfr=nypaor,1,-1
**      write(*,'(i6,1p,2d23.13)') jfr,v1afin(1,jfr),v1afin(nxpaor,jfr)
**    enddo
*
      END SUBROUTINE auvbcu
c
c***********************************************************************
c
      SUBROUTINE bcuini
*
*     Computes transformation matrices between data values and
*     bicubic fit coefficients, and the coefficient multipliers
*     at the finer resolution gridpoints, and makes them available
*     within the module (only) via the module's private storage
*     This version also multiplies the coefficient multipliers
*     stfn and the transform matrices b** for various cases,
*     once and for all, to avoid repeated computations in auvbcu

*     Modules

      USE atconst, ONLY : dya, bccoat

      IMPLICIT NONE

*     Subroutine arguments
*
*     Local parameters
      logical stfprt,bbprnt
      parameter ( stfprt = .false., bbprnt = .false. )
*
*     Local variables
      double precision stfn(16,0:ndxr,0:ndxr),
     &                 wfcn(-1:2,-1:2,0:1,0:1),wfnx(-1:2,-1:2,0:1,0:1),
     &                 wfny(-1:2,-1:2,0:1,0:1),wfxy(-1:2,-1:2,0:1,0:1),
     &                 bbb(16,16),bus(16,16),bun(16,16),
     &                 bvs(16,16),bvn(16,16),ss,tt,stbsum
      integer id,jd,ip,jp,i,j,k,m,ii,jj
*
*     Setup the coefficient multipliers on the refined grid
*     =====================================================
*     ss and tt are normalised coordinates (range = [0,1])
*     within/on the coarse gridcell being interpolated
*     Loop through the fine mesh gridpoints (ii, jj)
      do jj=0,ndxr
        tt = dble(jj)/dble(ndxr)
        do ii=0,ndxr
          ss = dble(ii)/dble(ndxr)
*         Loop through the bicubic terms, saving the results
*         as a 16-element vector at each fine res. gridpoint
          m = 0
          do j=0,3
            do i=0,3
              m = m + 1
              stfn(m,ii,jj) = ss**(i)*tt**(j)
            enddo
          enddo
        enddo
      enddo
      if ( stfprt ) then
*       Check stfn values
        print *,' '
        print *,' bcuini: stfn vector(s) for jj=1:'
        do m=1,16
          write (*,'(i5,9f10.6)') m,(stfn(m,ii,1),ii=0,ndxr)
        enddo
      endif
*
*     (ip, jp) is the vertex of the coarse resolution gridbox
*     at which the function value or derivative is being computed
*     (id, jd) is the location of the data value being accessed
*     to compute the function value or derivative required
*     wfcn contains the weights of the function values
*     wfnx contains the finite-difference contributions to dx*fx
*     wfny contains the finite-difference contributions to dy*fy
*     wfxy contains the finite-difference contributions to dx*dy*fxy
*
*     General case (u & v within channel)
*     ===================================
*
*     Specify the finite-difference operators for derivatives
*     -------------------------------------------------------
      do jp=0,1
        do ip=0,1
*         Initialise weight arrays to zero
          do jd=-1,2
            do id=-1,2
              wfcn(id,jd,ip,jp) = 0.0d0
              wfnx(id,jd,ip,jp) = 0.0d0
              wfny(id,jd,ip,jp) = 0.0d0
              wfxy(id,jd,ip,jp) = 0.0d0
            enddo
          enddo
*         Overwrite with nonzero values where appropriate
          wfcn( ip , jp ,ip,jp) =  1.0d0
          wfnx(ip+1, jp ,ip,jp) =  0.5d0
          wfnx(ip-1, jp ,ip,jp) = -0.5d0
          wfny( ip ,jp+1,ip,jp) =  0.5d0
          wfny( ip ,jp-1,ip,jp) = -0.5d0
          wfxy(ip+1,jp+1,ip,jp) =  0.25d0
          wfxy(ip-1,jp+1,ip,jp) = -0.25d0
          wfxy(ip+1,jp-1,ip,jp) = -0.25d0
          wfxy(ip-1,jp-1,ip,jp) =  0.25d0
        enddo
      enddo
*
*     Derive the B matrix from the weights wf**
*     -----------------------------------------
      call wts2bb (wfcn, wfnx, wfny, wfxy, bbb)
      if ( bbprnt ) then
*       Check B values
        print *,' '
        print *,' bcuini: B matrix (general case):'
        do i=1,16
          write (*,'(i5,16f6.2)') i,(bbb(i,j),j=1,16)
        enddo
      endif
*
*     Multiply the B matrix by the weights stfn
*     -----------------------------------------
      do jj=0,ndxr-1
        do ii=0,ndxr-1
          do k=1,16
            stbsum = 0.0d0
            do m=1,16
              stbsum = stbsum + bbb(m,k)*stfn(m,ii,jj)
            enddo
            stbbb(k,ii,jj) = stbsum
          enddo
        enddo
      enddo
*
*     u adjacent to southern boundary
*     ===============================
*
*     Specify the finite-difference operators for derivatives
*     -------------------------------------------------------
      do jp=0,1
        do ip=0,1
*         Initialise weight arrays to zero
          do jd=-1,2
            do id=-1,2
              wfcn(id,jd,ip,jp) = 0.0d0
              wfnx(id,jd,ip,jp) = 0.0d0
              wfny(id,jd,ip,jp) = 0.0d0
              wfxy(id,jd,ip,jp) = 0.0d0
            enddo
          enddo
*         Overwrite with nonzero values where appropriate
          wfcn( ip , jp ,ip,jp) =  1.0d0
          wfnx(ip+1, jp ,ip,jp) =  0.5d0
          wfnx(ip-1, jp ,ip,jp) = -0.5d0
          if ( jp.eq.0 ) then
*           Use mixed pressure BC
            wfny( ip , jp ,ip,jp) = (bccoat/dya)*wfcn( ip , jp ,ip,jp)
            wfxy(ip+1, jp ,ip,jp) = (bccoat/dya)*wfnx(ip+1, jp ,ip,jp)
            wfxy(ip-1, jp ,ip,jp) = (bccoat/dya)*wfnx(ip-1, jp ,ip,jp)
           else if ( jp.eq.1 ) then
*           Regular finite-difference
            wfny( ip ,jp+1,ip,jp) =  0.5d0
            wfny( ip ,jp-1,ip,jp) = -0.5d0
            wfxy(ip+1,jp+1,ip,jp) =  0.25d0
            wfxy(ip-1,jp+1,ip,jp) = -0.25d0
            wfxy(ip+1,jp-1,ip,jp) = -0.25d0
            wfxy(ip-1,jp-1,ip,jp) =  0.25d0
          endif
        enddo
      enddo
*
*     Derive the B matrix from the weights wf**
*     -----------------------------------------
      call wts2bb (wfcn, wfnx, wfny, wfxy, bus)
      if ( bbprnt ) then
*       Check B values
        print *,' '
        print *,' bcuini: B matrix (u near southern bdy):'
        do i=1,16
          write (*,'(i5,16f6.2)') i,(bus(i,j),j=1,16)
        enddo
      endif
*
*     Multiply the B matrix by the weights stfn
*     -----------------------------------------
      do jj=0,ndxr-1  !include 0 on the left to open the boundary
        do ii=0,ndxr-1
          do k=1,16
            stbsum = 0.0d0
            do m=1,16
              stbsum = stbsum + bus(m,k)*stfn(m,ii,jj)
            enddo
            stbus(k,ii,jj) = stbsum
          enddo
        enddo
      enddo
*
*     v adjacent to southern boundary
*     ===============================
*
*     Specify the finite-difference operators for derivatives
*     -------------------------------------------------------
      do jp=0,1
        do ip=0,1
*         Initialise weight arrays to zero
          do jd=-1,2
            do id=-1,2
              wfcn(id,jd,ip,jp) = 0.0d0
              wfnx(id,jd,ip,jp) = 0.0d0
              wfny(id,jd,ip,jp) = 0.0d0
              wfxy(id,jd,ip,jp) = 0.0d0
            enddo
          enddo
*         Overwrite with nonzero values where appropriate
          wfcn( ip , jp ,ip,jp) =  1.0d0
          wfnx(ip+1, jp ,ip,jp) =  0.5d0
          wfnx(ip-1, jp ,ip,jp) = -0.5d0
          if ( jp.eq.0 ) then
*           Use the u field on the boundary, since vy = -ux, vxy = -uxx
*           The extra u values are located at jd = -1
*           (which would otherwise be empty)
            wfny(ip+1,jp-1,ip,jp) = -wfnx(ip+1, jp ,ip,jp)
            wfny(ip-1,jp-1,ip,jp) = -wfnx(ip-1, jp ,ip,jp)
            wfxy(ip+1,jp-1,ip,jp) = -1.0d0
            wfxy( ip ,jp-1,ip,jp) =  2.0d0
            wfxy(ip-1,jp-1,ip,jp) = -1.0d0
           else if ( jp.eq.1 ) then
*           Regular finite-difference
            wfny( ip ,jp+1,ip,jp) =  0.5d0
            wfny( ip ,jp-1,ip,jp) = -0.5d0
            wfxy(ip+1,jp+1,ip,jp) =  0.25d0
            wfxy(ip-1,jp+1,ip,jp) = -0.25d0
            wfxy(ip+1,jp-1,ip,jp) = -0.25d0
            wfxy(ip-1,jp-1,ip,jp) =  0.25d0
          endif
        enddo
      enddo
*
*     Derive the B matrix from the weights wf**
*     -----------------------------------------
      call wts2bb (wfcn, wfnx, wfny, wfxy, bvs)
      if ( bbprnt ) then
*       Check B values
        print *,' '
        print *,' bcuini: B matrix (v near southern bdy):'
        do i=1,16
          write (*,'(i5,16f6.2)') i,(bvs(i,j),j=1,16)
        enddo
      endif
*
*     Multiply the B matrix by the weights stfn
*     -----------------------------------------
      do jj=0,ndxr-1
        do ii=0,ndxr-1
          do k=1,16
            stbsum = 0.0d0
            do m=1,16
              stbsum = stbsum + bvs(m,k)*stfn(m,ii,jj)
            enddo
            stbvs(k,ii,jj) = stbsum
          enddo
        enddo
      enddo
*
*     u adjacent to northern boundary
*     ===============================
*
*     Specify the finite-difference operators for derivatives
*     -------------------------------------------------------
      do jp=0,1
        do ip=0,1
*         Initialise weight arrays to zero
          do jd=-1,2
            do id=-1,2
              wfcn(id,jd,ip,jp) = 0.0d0
              wfnx(id,jd,ip,jp) = 0.0d0
              wfny(id,jd,ip,jp) = 0.0d0
              wfxy(id,jd,ip,jp) = 0.0d0
            enddo
          enddo
*         Overwrite with nonzero values where appropriate
          wfcn( ip , jp ,ip,jp) =  1.0d0
          wfnx(ip+1, jp ,ip,jp) =  0.5d0
          wfnx(ip-1, jp ,ip,jp) = -0.5d0
          if ( jp.eq.0 ) then
*           Regular finite-difference
            wfny( ip ,jp+1,ip,jp) =  0.5d0
            wfny( ip ,jp-1,ip,jp) = -0.5d0
            wfxy(ip+1,jp+1,ip,jp) =  0.25d0
            wfxy(ip-1,jp+1,ip,jp) = -0.25d0
            wfxy(ip+1,jp-1,ip,jp) = -0.25d0
            wfxy(ip-1,jp-1,ip,jp) =  0.25d0
           else if ( jp.eq.1 ) then
*           Use mixed pressure BC
            wfny( ip , jp ,ip,jp) = -(bccoat/dya)*wfcn( ip , jp ,ip,jp)
            wfxy(ip+1, jp ,ip,jp) = -(bccoat/dya)*wfnx(ip+1, jp ,ip,jp)
            wfxy(ip-1, jp ,ip,jp) = -(bccoat/dya)*wfnx(ip-1, jp ,ip,jp)
          endif
        enddo
      enddo
*
*     Derive the B matrix from the weights wf**
*     -----------------------------------------
      call wts2bb (wfcn, wfnx, wfny, wfxy, bun)
      if ( bbprnt ) then
*       Check B values
        print *,' '
        print *,' bcuini: B matrix (u near northern bdy):'
        do i=1,16
          write (*,'(i5,16f6.2)') i,(bun(i,j),j=1,16)
        enddo
      endif
*
*     Multiply the B matrix by the weights stfn
*     -----------------------------------------
      do jj=0,ndxr   !comment out -1 to open the boundary
        do ii=0,ndxr-1
          do k=1,16
            stbsum = 0.0d0
            do m=1,16
              stbsum = stbsum + bun(m,k)*stfn(m,ii,jj)
            enddo
            stbun(k,ii,jj) = stbsum
          enddo
        enddo
      enddo
*
*     v adjacent to northern boundary
*     ===============================
*
*     Specify the finite-difference operators for derivatives
*     -------------------------------------------------------
      do jp=0,1
        do ip=0,1
*         Initialise weight arrays to zero
          do jd=-1,2
            do id=-1,2
              wfcn(id,jd,ip,jp) = 0.0d0
              wfnx(id,jd,ip,jp) = 0.0d0
              wfny(id,jd,ip,jp) = 0.0d0
              wfxy(id,jd,ip,jp) = 0.0d0
            enddo
          enddo
*         Overwrite with nonzero values where appropriate
          wfcn( ip , jp ,ip,jp) =  1.0d0
          wfnx(ip+1, jp ,ip,jp) =  0.5d0
          wfnx(ip-1, jp ,ip,jp) = -0.5d0
          if ( jp.eq.0 ) then
*           Regular finite-difference
            wfny( ip ,jp+1,ip,jp) =  0.5d0
            wfny( ip ,jp-1,ip,jp) = -0.5d0
            wfxy(ip+1,jp+1,ip,jp) =  0.25d0
            wfxy(ip-1,jp+1,ip,jp) = -0.25d0
            wfxy(ip+1,jp-1,ip,jp) = -0.25d0
            wfxy(ip-1,jp-1,ip,jp) =  0.25d0
           else if ( jp.eq.1 ) then
*           Use the u field on the boundary, since vy = -ux, vxy = -uxx
*           The extra u values are located at jd = 2
*           (which would otherwise be empty)
            wfny(ip+1,jp+1,ip,jp) = -wfnx(ip+1, jp ,ip,jp)
            wfny(ip-1,jp+1,ip,jp) = -wfnx(ip-1, jp ,ip,jp)
            wfxy(ip+1,jp+1,ip,jp) = -1.0d0
            wfxy( ip ,jp+1,ip,jp) =  2.0d0
            wfxy(ip-1,jp+1,ip,jp) = -1.0d0
          endif
        enddo
      enddo
*
*     Derive the B matrix from the weights wf**
*     -----------------------------------------
      call wts2bb (wfcn, wfnx, wfny, wfxy, bvn)
      if ( bbprnt ) then
*       Check B values
        print *,' '
        print *,' bcuini: B matrix (v near northern bdy):'
        do i=1,16
          write (*,'(i5,16f6.2)') i,(bvn(i,j),j=1,16)
        enddo
      endif
*
*     Multiply the B matrix by the weights stfn
*     -----------------------------------------
      do jj=0,ndxr-1
        do ii=0,ndxr-1
          do k=1,16
            stbsum = 0.0d0
            do m=1,16
              stbsum = stbsum + bvn(m,k)*stfn(m,ii,jj)
            enddo
            stbvn(k,ii,jj) = stbsum
          enddo
        enddo
      enddo
*
      END SUBROUTINE bcuini
c
c***********************************************************************
c
      SUBROUTINE wts2bb (wfcn, wfnx, wfny, wfxy, bbb)
*
*     Given sets of weights for the function values (wfcn), and the
*     finite-difference approximations to the x, y and mixed derivative
*     terms (wfnx, wfny and wfxy), this routine derives the transfor-
*     mation matrix bbb between data values and bicubic fit coefficients
*     The weights are computed in the calling routine bcuini
*
      IMPLICIT NONE
*
*     Subroutine arguments
      double precision, INTENT(IN) ::
     &                 wfcn(-1:2,-1:2,0:1,0:1),wfnx(-1:2,-1:2,0:1,0:1),
     &                 wfny(-1:2,-1:2,0:1,0:1),wfxy(-1:2,-1:2,0:1,0:1)
      double precision, INTENT(OUT) :: bbb(16,16)
*
*     Local parameters
      logical u2fprt
      parameter ( u2fprt = .false. )
*
*     Local variables
      integer ip,jp,kp,id,jd,kd,i,j
      double precision u2f(16,16),stinv(16,16),wfsum
*
      SAVE stinv
*
      data stinv/  1,  0, -3,  2, 4*0, -3,  0,  9, -6,  2,  0, -6,  4,
     &             0,  0,  3, -2, 6*0, -9,  6,  0,  0,  6, -4,
     &           8*0,  3,  0, -9,   6, -2,  0,  6, -4,
     &          10*0,  9, -6,  0,   0, -6,  4,
     &             0,  1, -2,  1, 5*0, -3,  6, -3,  0,  2, -4,  2,
     &             0,  0, -1,  1, 6*0,  3, -3,  0,  0, -2,  2,
     &           9*0,  3, -6,  3,   0, -2,  4, -2,
     &          10*0, -3,  3,  0,   0,  2, -2,
     &           4*0,  1,  0, -3,   2, -2,  0,  6, -4,  1,  0, -3,  2,
     &           6*0,  3, -2,  0,   0, -6,  4,  0,  0,  3, -2,
     &           8*0, -1,  0,  3,  -2,  1,  0, -3,  2,
     &          10*0, -3,  2,  0,   0,  3, -2,
     &           5*0,  1, -2,  1,   0, -2,  4, -2,  0,  1, -2,  1,
     &           6*0, -1,  1,  0,   0,  2, -2,  0,  0, -1,  1,
     &           9*0, -1,  2, -1,   0,  1, -2,  1,
     &          10*0,  1, -1,  0,   0, -1,  1 /
*
*     Derive the u2f matrix from the supplied weights
*     -----------------------------------------------
*     u2f transforms a vector of 16 u or v data values into the
*     16-element vector containing {f, dx*fx, dy*fy, dx*dy*fxy},
*     where f is a vector of function values and fx, fy and fxy
*     are finite-difference approximations to its derivatives.
*     Within each vector/subvector we use standard Fortran ordering
*     kd, kp are subscripts within the (16 element) data
*     and successive (4 element) vertex lists respectively
*     i.e. function-vector(kp) = Sigma [ u2f(kp,kd)*data(kd) ]
*     These loops span all the elements of u2f, and the various wf**
*     arrays are fully zero-padded, so no need to zero u2f initially
      kp = 0
      do jp=0,1
        do ip=0,1
          kp = kp + 1
          kd = 0
          do jd=-1,2
            do id=-1,2
              kd = kd + 1
*             Function values
              u2f( kp ,kd) = wfcn(id,jd,ip,jp)
*             x-derivatives
              u2f(kp+4,kd) = wfnx(id,jd,ip,jp)
*             y-derivatives
              u2f(kp+8,kd) = wfny(id,jd,ip,jp)
*             mixed derivatives
              u2f(kp+12,kd) = wfxy(id,jd,ip,jp)
            enddo
          enddo
        enddo
      enddo
      if ( u2fprt ) then
*       Check u2f values
        print *,' '
        print *,' wts2bb: u2f matrix:'
        do i=1,16
          write (*,'(i5,16f6.2)') i,(u2f(i,j),j=1,16)
        enddo
      endif
*
*     Derive the B matrix = stinv*u2f
*     -------------------------------
*     B transforms a vector of 16 u or v values into
*     the 16-element vector containing the coeffts
*     Cij strung out in standard Fortran order
*     The second subscript of B operates on the data vector,
*     the first refers to the term in the bicibic series
*     i.e. coefficient-vector(ij) = Sigma [ bbb(ij,kd)*data(kd) ]
      do kd=1,16
        do i=1,16
          wfsum = 0.0d0
          do j=1,16
            wfsum = wfsum + stinv(i,j)*u2f(j,kd)
          enddo
          bbb(i,kd) = wfsum
        enddo
      enddo
*
      END SUBROUTINE wts2bb
c
c***********************************************************************
c
      END MODULE xfosubs
c
c***********************************************************************
