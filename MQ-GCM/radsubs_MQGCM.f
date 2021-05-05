c***********************************************************************
c     Q-GCM Version 1.5.0 : last modified 26/08/2013; 07/11-18/2020 (SK)
c                                                     08/10/2020    (SK)
c                                                     08/15/2020    (SK)
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
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c
c     A new version, including mean-state determination from radiative-
c     convective balance. From input parameters, zopt (optical depths) are
c     now replaced by layer emissivities. OML is assumed to be
c     a black-body radiators, the atmospheric layers (including a.m.l.) can
c     have emissivities <1. Some of the input parameters, fixed in the 
c     previous version of the Q-GCM code, are re-computed and overwritten
c     in this version; notably, this includes mean layer temperatures
c     and,optionally, reduced gravities. Constant AML thickness is assumed.
c     Topography contributions to radiation are neglected. 
c
c     Sergey Kravtsov 07/11/2020 - 07/18/2020
c
c     This version also includes a vertical diffusion parameterization
c     and non-zero entrainments at both of the atmospheric interfaces.
c     
c     Sergey Kravtsov 08/11/2020
c
c     The version of 08/15/2020 relaxes the assumption of a.m.l. emiissivity
c     being zm=1. Now zm can be < 1 (SK).
c
c
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------

c***********************************************************************

      MODULE radsubs

*     Contains subroutines radiat and trapin. radiat computes
*     and prints the values of various radiation coefficients,
*     derived as explained in Appendix A of the "Formulation
*     and users' guide for Q-GCM". trapin is used to compute
*     the integrals required for some of the coefficients.

c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c
c     The layer radiation is determined here based on the estimates 
c     using mean layer temperatures, so the computation of integrals
c     is no longer required. The trapin subroutine is, accordingly, 
c     removed.
c
c     Sergey Kravtsov 07/11/2020
c
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c


*     Modules

      IMPLICIT NONE

      PRIVATE

      PUBLIC  :: radiat

      CONTAINS

c***********************************************************************
c
      SUBROUTINE radiat
*
*     Computes the mean state radiative balance, including the
*     atmosphere and ocean mixed layer temperatures that ensure
*     equilibrium. Also computes the perturbation (linearised)
*     radiation coefficients A, B, C and D, the entrainment
*     coefficients derived from them, and the radiative balance
*     initialisation coefficients.

*     Modules
      USE parameters, ONLY : nxta,nyta, nla, nyto, nlo, fnot
      USE atconst, ONLY : yta,ytarel,hat,hta,
     &                tabsat,tmbara,tat,gpat,
     &                rhm,rh,rhmc,rhc,qatm,qat,dqat,Lv,
     &                rhom,rho
      USE occonst, ONLY : yto, ytorel, tabsoc, tmbaro, toc
      USE intrfac, ONLY : sstbar, astbar, hmat, hmoc, tnbdy, tsbdy, evapa
      USE xfosubs, ONLY : fsprim,esat
      USE radiate

      IMPLICIT NONE

*     Subroutine arguments
*
*     Local parameters
      logical prtval
      double precision stefan
*     stefan is the Stefan-Boltzmann constant (usually denoted by sigma)
      parameter ( stefan = 5.67040D-8,prtval = .false. )
*------------------------------------------------------------------------
*     added by SK (11/07/2020)
*
      integer Nsteps
      double precision tstep,Gc,xConv
*
*     Nsteps is the number of time steps to integrate layers temps. to
*     radiative-convective equilibrium, tstep is the time step (s), Gc is the
*     critical lapse rate for convection to occur (K/m), xConv (W/m^2K)
*     is the convective heat exchange coefficient 
*
      parameter ( Nsteps=12*24*365*100, tstep=300.d0,
     &                       Gc=6.5d-3, xConv=200.d0)
*
*     Local variables
      integer i,j,k,l,iter
      double precision F0upbar,Fmupbar,
     &                 Fupbar(nla),Fmdnbar,Fdnbar(nla),rrcpdt,
     &                 Tobar,Tmbar,Tbar(nla),Thmean,
     &        Ao,Am,A1,A2,A3,dTodt,dTmdt,dTdt(nla),
     &        Bm,BB(nla),c11,c12,c21,c22,ck,fmm
*----------------------------------------------------------------------
*----------------------------------------------------------------------
*  last three lines: Mean entrainment (convective) fluxes and layer actual 
*                (not potential) temp. mean in equilibrium, Thmean is the
*                mean potential temperature of the atmosphere, variables
*                A* and dT*dt are auxiliary, used in time integration,
*                Bm and BB are auxiliary workspace for coefficients in 
*                linearized radiation formulas, c11,c12,c21,c22,ck,fmm
*                are auxiliary coefficients for computing the radiative balance
*                initial conditions
*
*     added 07/11/2020-08/15/2020 (SK)
*----------------------------------------------------------------------
*----------------------------------------------------------------------
*
*------------------------------------------------------------------------
*     added by SK (11/08/2020)
*
      double precision xmu1_star,xmu2_star
      parameter(xmu1_star=1.0d-4,xmu2_star=1.0d-4)
*
*
*     xmu1_star,xmu2_star (m/s) - vertical diffusion entrainment rate
*     at interfaces 1/2 and 2/3, respectively
*
*------------------------------------------------------------------------

*
*     Variables for calculating radiation balance coefficients
      double precision rbtmat,rbtmoc,rbtmat1,xlamda1
      parameter(xlamda1=35.d0)
*
*     This still turns out to be too high. tsbdy is, instead, specified manually
*     to match observed SST distribution
*


*
*    xlamda1 is used here instead of (a smaller) xlambda to define
*    a more reasonable southern boundary radiative-convective equilibrium
*    temperature
*

*
*     Effective temperature differences between layers in the moist model
*     (for updated computation of interior entrainments there)
*     (added by S. Kravtsov, 10/21/2020)
      double precision dt1,dt2

      double precision g,RR,RV,eeps,Gtilde
      parameter(g=9.82d0,RR=287.0d0,RV=461.d0,eeps=RR/RV,Gtilde=g/(RR*Gc))

      double precision zz,TT,TTs,fact,es

*
*     the last three lines are used in the moist model to 
*     compute the climatological humidity distribution
*


*----------------------------------------------------------------------
*     The code starts here:
*----------------------------------------------------------------------


*
*
*
*    specify fixed evaporation (to be recomputed over ocean at each time step)

      do j=1,nyta
        do i=1,nxta
         evapa(i,j)=1.d0*(1.d0-0.9d0*dble(j-1)/(nyta-1))
     &            *1000.d0/(24.d0*365.d0*3600.d0) !m/s
        enddo
      enddo

*     Mean state radiative-convective balance
*
*     (assumes an aquaplanet, which makes no difference for zm=1,
*      but is a relevant statement for zm<1)
*
*     ============================

c
c     Initialize all equilibrium temps close to the "observed" values
c     (all in Kelvins)
c

      Tobar=288.d0
      Tmbar=Tobar-Gc*hmat/2
      Tbar(1)=Tobar-Gc*(hmat+hat(1)/2)
      Tbar(2)=Tobar-Gc*(hmat+hat(1)+hat(2)/2)
      Tbar(3)=Tobar-Gc*(hmat+hat(1)+hat(2)+hat(3)/2)
c
      do iter=1,Nsteps
c
c       radiative fluxes via temps and emissivities:
c
        Ao=stefan*Tobar**4
        Am=stefan*Tmbar**4
        A1=stefan*Tbar(1)**4
        A2=stefan*Tbar(2)**4
        A3=stefan*Tbar(3)**4
     
        F0upbar=Ao
        Fmupbar=F0upbar*(1-zm)+zm*Am
        Fupbar(1)=Fmupbar*(1-zopt(1))+zopt(1)*A1
        Fupbar(2)=Fupbar(1)*(1-zopt(2))+zopt(2)*A2
        Fupbar(3)=Fupbar(2)*(1-zopt(3))+zopt(3)*A3
        Fdnbar(3)=-zopt(3)*A3
        Fdnbar(2)=Fdnbar(3)*(1-zopt(2))-zopt(2)*A2
        Fdnbar(1)=Fdnbar(2)*(1-zopt(1))-zopt(1)*A1
        Fmdnbar=Fdnbar(1)*(1-zm)-zm*Am
c
c       convective fluxes:
c
c
        if ((Tobar-Tmbar)/(hmat/2).gt.Gc) then
          Flambar=xConv*(Tobar-Tmbar-Gc*hmat/2)
        else
          Flambar=0.d0
        endif
c
        if ((Tmbar-Tbar(1))/((hmat+hat(1))/2).gt.Gc) then
          Fmebar=xConv*(Tmbar-Tbar(1)-Gc*(hmat+hat(1))/2)
        else
          Fmebar=0.d0
        endif
c
        if ((Tbar(1)-Tbar(2))/((hat(1)+hat(2))/2).gt.Gc) then
          Febar(1)=xConv*(Tbar(1)-Tbar(2)-Gc*(hat(1)+hat(2))/2)
        else
          Febar(1)=0.d0
        endif
c
        if ((Tbar(2)-Tbar(3))/((hat(2)+hat(3))/2).gt.Gc) then
          Febar(2)=xConv*(Tbar(2)-Tbar(3)-Gc*(hat(2)+hat(3))/2)
        else
          Febar(2)=0.d0
        endif
c
c       temperature tendencies:
c
c
        dTodt=(-fsbar-F0upbar-Fmdnbar-Flambar)*rrcpoc/hmoc
        dTmdt=(F0upbar+Fmdnbar-Fmupbar-Fdnbar(1)
     &         +Flambar-Fmebar)*rrcpat/hmat
        dTdt(1)=(Fdnbar(1)-Fupbar(1)+Fmupbar-Fdnbar(2)
     &         +Fmebar-Febar(1))*rrcpat/hat(1)
        dTdt(2)=(Fdnbar(2)-Fupbar(2)+Fupbar(1)-Fdnbar(3)
     &         +Febar(1)-Febar(2))*rrcpat/hat(2)
        dTdt(3)=
     &     (Fdnbar(3)-Fupbar(3)+Fupbar(2)+Febar(2))*rrcpat/hat(3)
c
c       time step:
c
       Tobar=Tobar+dTodt*tstep 
       Tmbar=Tmbar+dTmdt*tstep
       Tbar(1)=Tbar(1)+dTdt(1)*tstep
       Tbar(2)=Tbar(2)+dTdt(2)*tstep
       Tbar(3)=Tbar(3)+dTdt(3)*tstep
c
      enddo

c

*
*     Recompute all layer (potential!) temps accordingly
*     ============================

      tmbaro=Tobar
      tmbara=Tmbar+gamma*hmat/2
      tabsat(1)=Tbar(1)+gamma*(hmat+hat(1)/2)
      tabsat(2)=Tbar(2)+gamma*(hmat+hat(1)+hat(2)/2)
      tabsat(3)=Tbar(3)+gamma*(hmat+hat(1)+hat(2)+hat(3)/2)
      Thmean=(tmbara*hmat+tabsat(1)*hat(1)
     &                   +tabsat(2)*hat(2)
     &                   +tabsat(3)*hat(3))/(hmat+hta)

*     Recompute reduced gravities
*     ============================

      gpat(1)=9.82d0*(tabsat(2)-tabsat(1))/Thmean
      gpat(2)=9.82d0*(tabsat(3)-tabsat(2))/Thmean

c--------------------------------------------------------------
c     Ratio of radiatively driven entrainment rates of the
c     1/2 vs. 2/3 interfaces. This is chosen to ensure comparable
c     vertical velocity shears in these layers (per thermal-wind
c     balance).
c 
      f2=gpat(1)/gpat(2) 
c
c     added by SK 11/08/2020
c--------------------------------------------------------------

c--------------------------------------------------------------
c
c     Oceanic temperatures
c
c-------------------------------------------------------------

      tabsoc(1)=Tobar-2.d0
      tabsoc(2)=Tobar-10.d0
      tabsoc(3)=Tobar-14.d0


*     Derive temperature anomalies
      do k=1,nlo
        toc(k) = tabsoc(k) - tmbaro
      enddo
      do k=1,nla
        tat(k) = tabsat(k) - tmbara
      enddo
c
c-------------------------------------------------------------
c
c  climatological relative humidity distribution and
c  specific humidity differences between 1/2 and 2/3
c  atmopsheric layers
c  (only used in the moist version of the model)
c
c    Sergey Kravtsov (10/27/2020)
c
c  Start with the relative humidity linear z-profile
c  mimicking the observed distribution
c  rh=0.8-0.35*z/H, with z in km and H=10km
c
c
      rhm=0.8d0-0.35d0*(hmat/2.d0)/(hmat+hta)
      rh(1)=0.8d0-0.35d0*(hmat+hat(1)/2.d0)/(hmat+hta)
      rh(2)=0.8d0-0.35d0*(hmat+hat(1)+hat(2)/2.d0)/(hmat+hta)
      rh(3)=0.8d0-0.35d0*(hmat+hat(1)+hat(2)+hat(3)/2.d0)/(hmat+hta)
c
c Now use the atmo. temperature profile computed  above
c  (with constant lapse rate ~ Gc) and the relative humidities above
c  to compute the climatological moisture distribution
c

      TTs=tmbara-0.5d0*(gamma-Gc)*hmat
      zz = hmat/2.d0
      TT=TTs-Gc*zz
      fact=(TTs/TT)**Gtilde
      es=esat(Tmbar)
      qatm=0.001d0*fact*eeps*rhm*es
      rhom=(1.d5/(fact*RR*TT))

      zz = hmat+hat(1)/2.d0
      TT=TTs-Gc*zz
      fact=(TTs/TT)**Gtilde
      es=esat(Tbar(1))
      qat(1)=0.001d0*fact*eeps*rh(1)*es
      rho(1)=(1.d5/(fact*RR*TT))


      zz = hmat+hat(1)+hat(2)/2.d0
      TT=TTs-Gc*zz
      fact=(TTs/TT)**Gtilde
      es=esat(Tbar(2))
      qat(2)=0.001d0*fact*eeps*rh(2)*es
      rho(2)=(1.d5/(fact*RR*TT))

      zz = hmat+hat(1)+hat(2)+hat(3)/2.d0
      TT=TTs-Gc*zz
      fact=(TTs/TT)**Gtilde
      es=esat(Tbar(3))
      qat(3)=0.001d0*fact*eeps*rh(3)*es
      rho(3)=(1.d5/(fact*RR*TT))

      dqat(1)=qat(1)-qat(2)  !2.5d-3
      dqat(2)=qat(2)-qat(3)  !1.25d-3

      rhmc=rhm+0.2d0 !1.d0
      rhc(1)=rh(1)+0.2d0 !1.d0
      rhc(2)=rh(2)+0.2d0 !1.d0
      rhc(3)=rh(3)+0.2d0 !1.d0

c------------------------------------------------------------

      print *,' '
      print *,' Mean state radiation parameters:'
      print *,' --------------------------------'
      write(*,205) '  Mean forcing fsbar (W m^-2) = ',fsbar
      write(*,205) '  Pert. ampl. fspamp (W m^-2) = ',fspamp
*     Derive signed perturbation coefficient
      fspco = sign( fspamp, fnot)
      if ( fspamp.lt.0.0d0 ) then
        print *,' Sign of fspamp appears incorrect'
        print *,' Program terminates in radiat'
        stop
      endif
      write(*,205) '  Pert. coefft fspco (W m^-2) = ',fspco
      write(*,204) '  Emissivity in aml. zm (dmls) = ',zm
      do k=1,nla
        write(*,224)'  Emissivity zopt in layer',k,' = ',zopt(k)
      enddo
      write(*,214) '  Lapse rate gamma   (K m^-1) = ',gamma
      print *,' '
      write(*,251) '  Atmos. mixed layer T (K, C) = ',
     &             tmbara,tmbara-273.15d0
      write(*,251) '  Ocean  mixed layer T (K, C) = ',
     &             tmbaro,tmbaro-273.15d0
      write(*,207) '  Rel. ocean  layer temp. toc = ',(toc(k),k=1,nlo)
      write(*,207) '  Rel. atmos. layer temp. tat = ',(tat(k),k=1,nla)
      write(*,207) '  Reduced gravities  (m s^-2) = ',
     &             (gpat(k),k=1,nla-1)
      print *,' '
      write(*,251) '  Ocean m.l. radiat.  F0upbar = ',F0upbar
      write(*,251) '  Mixed layer  Fmbar up, down = ',Fmupbar,Fmdnbar
      do k=1,nla
        write(*,271) '  Layer', k ,'       Fbar up, down = ',
     &               Fupbar(k),Fdnbar(k)
      enddo
      write(*,251) '  Convective flux o-a  Flambar = ',Flambar
      write(*,251) '  Convective flux a-1  Fmebar = ',Fmebar
      write(*,251) '  Convective flux 1-2  Febar(1) = ',Febar(1)
      write(*,251) '  Convective flux 2-3  Febar(2) = ',Febar(2)
      write(*,214) '  Fractional error in OLR     = ',
     &             abs( Fupbar(nla) + fsbar )/abs(fsbar)


      print *,' '
      print *,' Mean state humidity parameters and thresholds:'
      print *,' --------------------------------'
      
      write(*,204) '  AML  density (kg m^-3) = ',rhom
      do k=1,nla
        write(*,224)'  Air density in layer',k,' = ',rho(k)
      enddo

      write(*,204) '  AML  relative humidity = ',rhm
      do k=1,nla
        write(*,224)'  Relative humidity in layer',k,' = ',rh(k)
      enddo
      write(*,204) '  AML relative humidity threshold = ',rhmc
      do k=1,nla
        write(*,224)'  RH threshold in layer',k,' = ',rhc(k)
      enddo
      print *,' '

      write(*,204) '  AML specific humidity (kg kg^-1) = ',qatm
      do k=1,nla
        write(*,224)'  Specific humidity in layer',k,' = ',qat(k)
      enddo

      write(*,207) '  Specific hum. diff.  (kg kg^-1) = ',
     &             (dqat(k),k=1,nla-1)
      print *,' '


*     Perturbation (linearised) radiation parameters
*     ==============================================
*     Initialise Aup and Adown arrays
      do l=1,nla-1
        do k=1,nla
          Aup(k,l) = 0.0d0
          Adown(k,l) = 0.0d0
        enddo
      enddo
c
c--------------------------------------------------------------
c
c     Auxiliary coefficients for linearized radiation formulas
c
      Bm=4.0d0*stefan*Tmbar**3
      do l=1,nla
         BB(l)=4.0d0*stefan*Tbar(l)**3
      enddo
c
*     Upgoing
*     -------
*     Ocean mixed layer
c
      D0up = 4.0d0*stefan*tmbaro**3

*     Atmospheric mixed layer

      Bmup = 0.d0
      Cmup = Bmup
      Dmup = zm*Bm
      Emup = (1.d0-zm)*D0up

*     Layer 1

      Aup(1,1) = -zopt(1)*BB(1)*(tat(2)-tat(1))/hta 
      Aup(1,2) = -zopt(1)*BB(1)*(tat(3)-tat(2))/hta
      Bup(1) =  Bmup 
      Cup(1) =  Cmup 
      Dup(1) = Dmup*(1.d0-zopt(1))
      Eup(1) = Emup*(1.d0-zopt(1))

*     Upper layers

      Aup(2,1)=-((1-zopt(2))*zopt(1)*BB(1)+zopt(2)*BB(2))
     &         *(tat(2)-tat(1))/hta
      Aup(2,2)=-((1-zopt(2))*zopt(1)*BB(1)+zopt(2)*BB(2))
     &         *(tat(3)-tat(2))/hta
      Aup(3,1)=-((1-zopt(3))*((1-zopt(2))*zopt(1)*BB(1)
     &         +zopt(2)*BB(2))+zopt(3)*BB(3))
     &         *(tat(2)-tat(1))/hta
      Aup(3,2)=-((1-zopt(3))*((1-zopt(2))*zopt(1)*BB(1)
     &         +zopt(2)*BB(2))+zopt(3)*BB(3))
     &         *(tat(3)-tat(2))/hta


      do k=2,nla
        Bup(k) = Bup(k-1)
        Cup(k) = Cup(k-1)
        Dup(k) = Dup(k-1)*(1.d0-zopt(k))
        Eup(k) = Eup(k-1)*(1.d0-zopt(k))
      enddo


*     Downgoing
*     ---------
*     Top layer
      
      Adown(nla,1) = zopt(nla)*BB(nla)*(tat(2)-tat(1))/hta
      Adown(nla,2) = zopt(nla)*BB(nla)*(tat(3)-tat(2))/hta

*     Intermediate layers

      Adown(2,1)=((1-zopt(2))*zopt(3)*BB(3)+zopt(2)*BB(2))
     &         *(tat(2)-tat(1))/hta
      Adown(2,2)=((1-zopt(2))*zopt(3)*BB(3)+zopt(2)*BB(2))
     &         *(tat(3)-tat(2))/hta
      Adown(1,1)=((1-zopt(1))*((1-zopt(2))*zopt(3)*BB(3)
     &         +zopt(2)*BB(2))+zopt(1)*BB(1))
     &         *(tat(2)-tat(1))/hta
      Adown(1,2)=((1-zopt(1))*((1-zopt(2))*zopt(3)*BB(3)
     &         +zopt(2)*BB(2))+zopt(1)*BB(1))
     &         *(tat(3)-tat(2))/hta


      B1down = 0.d0
      C1down = B1down

*     Atmospheric mixed layer

      Dmdown = - Dmup

      print *,' '
      print *,' Linearised radiation coeffts:'
      print *,' -----------------------------'
      print *,' QG internal interface eta coeffts:'
      do k=1,nla
        write(*,235) '  Layer', k ,' coeffts    Aup(k,l) = ',
     &               (Aup(k,l),l=1,nla-1)
      enddo
      do k=1,nla
        write(*,235) '  Layer', k ,' coeffts  Adown(k,l) = ',
     &               (Adown(k,l),l=1,nla-1)
      enddo
      print *,' Mixed layer eta coeffts:'
      write(*,253) '  Mixed layer coefft     Bmup = ',Bmup
      write(*,253) '  Layer 1 coefft       B1down = ',B1down
      do k=1,nla
        write(*,273) '  Layer', k ,' coefft          Bup = ',Bup(k)
      enddo
      print *,' Atmos. topography coeffts:'
      write(*,253) '  Mixed layer coefft     Cmup = ',Cmup
      write(*,253) '  Layer 1 coefft       C1down = ',C1down
      do k=1,nla
        write(*,273) '  Layer', k ,' coefft          Cup = ',Cup(k)
      enddo
      print *,' Temperature perturbation coeffts:'
      write(*,253) '  Oceanic coefft         D0up = ',D0up
      write(*,253) '  Mixed layer coefft     Dmup = ',Dmup
      write(*,253) '  Mixed layer coefft   Dmdown = ',Dmdown
      do k=1,nla
        write(*,273) '  Layer', k ,' coefft          Dup = ',Dup(k)
      enddo
      write(*,253) '  Mixed layer coefft     Emup = ',Emup
      do k=1,nla
        write(*,273) '  Layer', k ,' coefft          Eup = ',Eup(k)
      enddo

*     Entrainment factors: aface(l), bface, cface and dface such that
*     e(1) = Sum(l)[ aface(l)*eta(l) ] + bface*etam + cface*aD + dface*aTm'
*     We are assuming all entrainment is across interface 1.

***************************************************************************
*
*    The above arrangement is changed as follows:
*
*    e(1) = aface(1)*eta(1) + aface(2)*eta(2) 
*         + bface*etam + cface*aD + dface*aTm'+eface*oTm'
*    e(2) = aface1(1)*eta(1) + aface1(2)*eta(2) + dface1*aTm'
*                                               + eface1*oTm'
*
*    Note that bface=cface=0 in this version of the code, in which we assume
*    constant hmixa and neglect contribution to radiative forcing due to presence
*    of topography. The coefficients aface and aface1 include diffusive contributions.
*
*    Sergey Kravtsov, 07/18/2020 - 08/11/2020 - 08/15/2020
***************************************************************************

      rrcpdt = rrcpat/( tat(2)-tat(1) )
      bface = rrcpdt*( B1down + Bmup - Bup(nla) )
      cface = rrcpdt*( C1down + Cmup - Cup(nla) )

        dt1=tat(2)-tat(1)
        dt2=tat(3)-tat(2)

      aface(1)=((Adown(1,1)-Aup(3,1))*rrcpat
     & +xmu2_star*dt2/hat(2)
     & -xmu1_star*dt1*(1.d0/hat(1)+1.d0/hat(2)))
     &           /(dt1+f2*dt2)
     & +xmu1_star*(1.d0/hat(1)+1.d0/hat(2))

      aface(2)=((Adown(1,2)-Aup(3,2))*rrcpat
     & +xmu1_star*dt1/hat(2)
     & -xmu2_star*dt2*(1.d0/hat(2)+1.d0/hat(3)))
     &           /(dt1+f2*dt2)
     & -xmu1_star/hat(2)

      aface1(1)=f2*((Adown(1,1)-Aup(3,1))*rrcpat
     & +xmu2_star*dt2/hat(2)
     & -xmu1_star*dt1*(1.d0/hat(1)+1.d0/hat(2)))
     &           /(dt1+f2*dt2)
     & -xmu2_star/hat(2)

      aface1(2)=f2*((Adown(1,2)-Aup(3,2))*rrcpat
     & +xmu1_star*dt1/hat(2)
     & -xmu2_star*dt2*(1.d0/hat(2)+1.d0/hat(3)))
     &           /(dt1+f2*dt2)
     & +xmu2_star*(1.d0/hat(2)+1.d0/hat(3))

      dface=(Dmup-Dup(nla))*rrcpat
     &     /(dt1+f2*dt2)
      dface1=f2*dface

      eface=(Emup-Eup(nla))*rrcpat
     &     /(dt1+f2*dt2)
      eface1=f2*eface

*
*     over land (08/16/2020):
*
*     these corrections are added *on top of* the previous
*     formulas with aface and dface
*
      epsfac=1.d0-(1.d0-zopt(1))*(1.d0-zopt(2))*(1.d0-zopt(3))

      afacel(1)=-epsfac*(1.d0-zm)**2*Adown(1,1)*rrcpat
     &           /(dt1+f2*dt2)

      afacel(2)=-epsfac*(1.d0-zm)**2*Adown(1,2)*rrcpat
     &           /(dt1+f2*dt2)

      aface1l(1)=-epsfac*f2*(1.d0-zm)**2*Adown(1,1)*rrcpat
     &           /(dt1+f2*dt2)

      aface1l(2)=-epsfac*f2*(1.d0-zm)**2*Adown(1,2)*rrcpat
     &           /(dt1+f2*dt2)

      dfacel=-epsfac*(1.d0-zm)*Dmdown*rrcpat
     &           /(dt1+f2*dt2)
      dface1l=-epsfac*f2*(1.d0-zm)*Dmdown*rrcpat
     &           /(dt1+f2*dt2)


*     Radiation balance initialisation coefficients
*     =============================================
*     We assume the mixed layer interface displacement etam
*     and atmospheric topography are both zero, and
*     that the entrainments between QG layers are all zero
*
*     Note that this balance also assumes an aquaplanet:
*
c

      c11=-(dface*aface1(2)-dface1*aface(2))
     &  /(aface(1)*aface1(2)-aface(2)*aface1(1))
      c21=-(dface1*aface(1)-dface*aface1(1))
     &  /(aface(1)*aface1(2)-aface(2)*aface1(1))
      c12=-(eface*aface1(2)-eface1*aface(2))
     &  /(aface(1)*aface1(2)-aface(2)*aface1(1))
      c22=-(eface1*aface(1)-eface*aface1(1))
     &  /(aface(1)*aface1(2)-aface(2)*aface1(1))

      fmm=-zm/(Emup*zm+xlamda1) !note xlamda1 here and elsewhere

      ck=-(zm*Dmup+(Dmdown-Dmup)-xlamda1)/(Emup*zm+xlamda1)

      rbtmat=-((Emup+Adown(1,1)*c12+Adown(1,2)*c22)*fmm+1.d0)
     & /((Emup+Adown(1,1)*c12+Adown(1,2)*c22)*ck
     &  +Dmup+Adown(1,1)*c11+Adown(1,2)*c21)

      rbtmoc=fmm+ck*rbtmat
      
      rbetat(1)=c11*rbtmat+c12*rbtmoc
      rbetat(2)=c21*rbtmat+c22*rbtmoc

*
*
*     Now over land (08/22/2020)
*
*
      c11=-((dface+dfacel)*(aface1(2)+aface1l(2))
     &    -(dface1+dface1l)*(aface(2)+afacel(2)))
     &     /((aface(1)+afacel(1))*(aface1(2)+aface1l(2))
     &    -(aface(2)+afacel(2))*(aface1(1)+aface1l(1)))

      c21=-((dface1+dface1l)*(aface(1)+afacel(1))
     &    -(dface+dfacel)*(aface1(1)+aface1l(1)))
     &     /((aface(1)+afacel(1))*(aface1(2)+aface1l(2))
     &    -(aface(2)+afacel(2))*(aface1(1)+aface1l(1)))

      c12=((1.d0-zm)*epsfac*rrcpat
     &       /((tat(2)-tat(1))+f2*(tat(3)-tat(2))))
     &    *((aface1(2)+aface1l(2))
     &    -f2*(aface(2)+afacel(2)))
     &     /((aface(1)+afacel(1))*(aface1(2)+aface1l(2))
     &    -(aface(2)+afacel(2))*(aface1(1)+aface1l(1)))

      c22=((1.d0-zm)*epsfac*rrcpat
     &     /((tat(2)-tat(1))+f2*(tat(3)-tat(2))))
     &    *(f2*(aface(1)+afacel(1))
     &    -(aface1(1)+aface1l(1)))
     &     /((aface(1)+afacel(1))*(aface1(2)+aface1l(2))
     &    -(aface(2)+afacel(2))*(aface1(1)+aface1l(1)))


      rbtmat1=-zm*(1.d0+(2.d0-zm)*(Adown(1,1)*c12+Adown(1,2)*c22))
     &      /(zm*(2.d0-zm)*(Adown(1,1)*c11+Adown(1,2)*c21)
     &      +(Dmup-(1.d0-zm)*Dmdown))


      print *,' '
      print *,' Radiative balance initialisation coeffts:'
      print *,' -----------------------------------------'
      print *,' QG internal interface eta coeffts:'
      do k=1,nla-1
        write(*,271) '  Layer', k ,' coefft       rbetat = ',rbetat(k)
      enddo
      write(*,251) '  Atmos. m.l. T coefft rbtmat = ',rbtmat
      write(*,251) '  Ocean  m.l. T coefft rbtmoc = ',rbtmoc
      write(*,251) '  Atmos. m.l. T ovrlnd rbtmat1 = ',rbtmat1

*     Perturbed state radiative balance
*     =================================
*     Initialise atmos. mixed layer rel. temp.
      do j=1,nyta
        astbar(j) = rbtmat*fsprim( ytarel(j) )
      enddo
      if ( prtval ) then
        print *,' '
        print *,' Atmos. m. l. temperature initialisation:'
        print *,'    j     yta(km)          fsprim             astbar'
        do j=nyta,1,-1
          write(*,'(i6,f13.4,2f19.12)')
     &          j,1.0d-3*yta(j),fsprim(ytarel(j)),astbar(j)
        enddo
        print *,' '
      endif
      write(*,251) '  A.m.l. astbar range  S -> N = ',
     &             astbar(1),astbar(nyta)

*     Initialise oceanic mixed layer rel. temp.
      do j=1,nyto
        sstbar(j) = rbtmoc*fsprim( ytorel(j) )
      enddo
      if ( prtval ) then
        print *,' '
        print *,' Ocean m. l. temperature initialisation:'
        print *,'    j     yto(km)          fsprim             sstbar'
        do j=nyto,1,-1
          write(*,'(i6,f13.4,2f19.12)')
     &          j,1.0d-3*yto(j),fsprim(ytorel(j)),sstbar(j)
        enddo
        print *,' '
      endif
      write(*,251) '  O.m.l. sstbar range  S -> N = ',
     &             sstbar(1),sstbar(nyto)

      print *,' '
      print *,' Radiative balance boundary temperatures:'
*     Compute Northern boundary temp. to maintain equilibrium
*     Tends to keep 1st internal T point at radiative balance
      tnbdy = sstbar(nyto)
      write(*,253) '  Rel. o.m.l. N. bndry temp.  = ',tnbdy
*     Compute Southern boundary temp. to maintain equilibrium
*     Tends to keep 1st internal T point at radiative balance
*      tsbdy = sstbar(1)
      tsbdy = 12.d0
      write(*,253) '  Rel. o.m.l. S. bndry temp.  = ',tsbdy

*     Entrainment factors:

      print *,' '
      print *,' Entrainment coefficients:'
      print *,' -------------------------'
      print *,' QG internal interface coeffts:'
      write(*,215) '  eta1  coefficients  aface(l) = ',
     &             (aface(l),l=1,nla-1)
      write(*,215) '  eta2  coefficients  aface(l) = ',
     &             (aface1(l),l=1,nla-1)
      write(*,215) '  eta1  coefficients  afacel(l) = ',
     &             (afacel(l),l=1,nla-1)
      write(*,215) '  eta2  coefficients  afacel(l) = ',
     &             (aface1l(l),l=1,nla-1)
      write(*,260) '  etam   coefficient    bface = ',bface
      write(*,260) '  topog. coefficient    cface = ',cface
      write(*,260) '  aTm1   coefficient    dface = ',dface
      write(*,260) '  aTm2   coefficient    dface1 = ',dface1
      write(*,260) '  aTm1   coefficient    dfacel = ',dfacel
      write(*,260) '  aTm2   coefficient    dface1l = ',dface1l
      write(*,260) '  oTm   coefficient1    eface = ',eface
      write(*,260) '  oTm   coefficient2    eface1 = ',eface1



  204 format(a,9f13.4)
  205 format(a,9f13.5)
* 206 format(a,9f13.6)
  207 format(a,9f13.7)
* 208 format(a,9f13.8)
* 209 format(a,9f13.9)
  214 format(a,1p,9d13.4)
  215 format(a,1p,9d13.5)
  224 format(a,i2,a,9f13.4)
* 226 format(a,i2,a,9f13.6)
* 228 format(a,i2,a,9f13.8)
* 229 format(a,i2,a,9f13.9)
  235 format(a,i2,a,1p,9d13.5)
* Revised format statements with extra precision and spacing
  251 format(a,2f18.11)
  253 format(a,2f18.13)
  260 format(a,1p,2d18.10)
  271 format(a,i2,a,2f18.11)
  272 format(a,i2,a,2f18.12)
  273 format(a,i2,a,2f18.13)

      END SUBROUTINE radiat
c
c***********************************************************************
c
      END MODULE radsubs
c
c***********************************************************************
