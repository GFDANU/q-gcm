c***********************************************************************
c     Q-GCM Version 1.5.0 : last modified 18/02/2013
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
      MODULE omlsubs

*     Contains subroutines oml and omladf which step the oceanic
*     mixed layer temperature evolution equations (7.11), and update
*     the entrainment between oceanic layers given by (7.12),
*     with possible convective correction given by (7.13).

*     Modules

      IMPLICIT NONE

      PRIVATE

#ifndef atmos_only
      PUBLIC  :: oml
      PRIVATE :: omladf

      CONTAINS
#endif

c***********************************************************************
c
#ifndef atmos_only

      SUBROUTINE oml
*
*     Timestep oceanic mixed layer temperature, eqn (7.11)
*     Also compute entrainment between oceanic layers (7.12-7.13)

*     Modules
      USE parameters
      USE occonst
      USE ochomog
      USE ocstate
      USE intrfac, ONLY : sst, sstm, fnetoc, hmoc
      USE radiate, ONLY : rrcpoc
      USE monitor, ONLY : cfraoc, centoc
      USE intsubs, ONLY : xintp

      IMPLICIT NONE

*     Subroutine arguments
*
*     Local variables
      integer i,j
      double precision hmoinv,dtoinv,entfac,rhs(nxto,nyto),
     &                 xfosum,cfrasm,centsm,xfsi,xfo(nxto,nyto),
     &                 xfoent,diabat,sstnew,coneno,dtonew
#  ifdef cyclic_ocean
      double precision ensums,ensumn
#  endif

      hmoinv = 1.0d0/hmoc
      dtoinv = 1.0d0/( toc(1)-toc(2) )
      entfac = hmoc*dtoinv/tdto

*     Initialise rhs with advective and diffusive
*     terms; first 4 terms in equation (7.11)
*     -------------------------------------------
      call omladf (rhs, po(1,1,1))

*     Add forcing term and diabatic effect - fifth and
*     sixth terms in equation (7.11), then timestep sst
*     -------------------------------------------------
      xfosum = 0.0d0
*     MONITORING - extra section for cfraoc, centoc
      cfrasm = 0.0d0
      centsm = 0.0d0

!$OMP PARALLEL DEFAULT (NONE)
!$OMP&         PRIVATE (i,j,xfsi)
!$OMP&         SHARED  (wekto,sstm,toc,tdto,rhs,hmoinv,rrcpoc,fnetoc,
!$OMP&                  dtoinv,entfac,xfo,cfrasm,centsm,sst,xfosum,
!$OMP&                  entoc)

!$OMP DO SCHEDULE (STATIC)
!$OMP&   PRIVATE (diabat,sstnew,xfoent,dtonew,coneno)
!$OMP&   REDUCTION(+:cfrasm) REDUCTION(-:centsm)
      do j=1,nyto
        do i=1,nxto
*         Predict new sst - equation (7.11)
          diabat = 0.5d0*wekto(i,j)*( sstm(i,j) + toc(1) )
          sstnew = sstm(i,j) + tdto*( rhs(i,j) + hmoinv*
     &                               ( rrcpoc*fnetoc(i,j) + diabat ) )

*         Find layer 1/2 entrainment at T points - equation (7.12)
          xfoent = -(0.5d0*dtoinv)*wekto(i,j)*( sstm(i,j) - toc(1) )

*         Check for convection & if necessary correct
*         layer 1/2 entrainment and sst - equation (7.13)
*         dtonew should be <= 0 (stable case)
*         Correction is nonzero only if dtonew > 0
          dtonew = toc(1) - sstnew
          coneno = entfac*max( 0.0d0, dtonew )
          xfo(i,j) = xfoent - coneno
          sstnew = sstnew + max( 0.0d0, dtonew )

*         MONITORING - extra section for cfraoc, centoc
          cfrasm = cfrasm + ( 0.5d0 - sign( 0.5d0, -dtonew ) )
          centsm = centsm - coneno

          sstm(i,j) = sst(i,j)
          sst(i,j) = sstnew
        enddo
      enddo
!$OMP END DO

*     To appease the Intel v8.1 compiler/Itanium 2 processors,
*     compute xfosum in a separate loop. In principle it should
*     be possible to include the code of the following loop in
*     the previous loop, but this gives the wrong answer when
*     running in parallel on some systems. The following separate
*     loop should thus be regarded purely as a temporary fix.
!$OMP DO SCHEDULE (STATIC)
!$OMP&   REDUCTION(+:xfosum)
      do j=1,nyto
        xfsi = 0.0d0
        do i=1,nxto
          xfsi = xfsi + xfo(i,j)
        enddo
        xfosum = xfosum + xfsi
      enddo
!$OMP END DO

*     Correct xfo values by subtracting mean to
*     ensure net entrainment is zero, implying
*     no net heat flux into the deep ocean.
!$OMP DO SCHEDULE (STATIC)
      do j=1,nyto
        do i=1,nxto
          xfo(i,j) = xfo(i,j) - xfosum*ocnorm
        enddo
      enddo
!$OMP END DO

*     Average entrainment (computed in xfo at T points) onto p points
*     ---------------------------------------------------------------
*     Internal points
!$OMP DO SCHEDULE (STATIC)
      do j=2,nypo-1
        do i=2,nxpo-1
          entoc(i,j) = 0.25d0*(  xfo(i-1,j-1) + xfo(i,j-1)
     &                         + xfo(i-1, j ) + xfo(i, j ) )
        enddo
      enddo
!$OMP END DO NOWAIT

*     Do edge and corner values to conserve area integral
!$OMP DO SCHEDULE (STATIC)
      do i=2,nxpo-1
        entoc(i,  1 ) = 0.5d0*( xfo(i-1,  1 ) + xfo(i,  1 ) )
        entoc(i,nypo) = 0.5d0*( xfo(i-1,nyto) + xfo(i,nyto) )
      enddo
!$OMP END DO NOWAIT
#  ifdef cyclic_ocean
*     Zonally cyclic ocean
!$OMP DO SCHEDULE (STATIC)
      do j=2,nypo-1
        entoc(  1 ,j) = 0.25d0*(  xfo(nxto,j-1) + xfo(1,j-1)
     &                          + xfo(nxto, j ) + xfo(1, j ) )
        entoc(nxpo,j) = entoc(  1 ,j)
      enddo
!$OMP END DO NOWAIT
!$OMP SINGLE
      entoc(  1 ,  1 ) = 0.5d0*( xfo(nxto,  1 ) + xfo(1,  1 ) )
      entoc(  1 ,nypo) = 0.5d0*( xfo(nxto,nyto) + xfo(1,nyto) )
      entoc(nxpo,  1 ) = entoc(  1 ,  1 )
      entoc(nxpo,nypo) = entoc(  1 ,nypo)
!$OMP END SINGLE NOWAIT
#  else
*     Finite box ocean
!$OMP DO SCHEDULE (STATIC)
      do j=2,nypo-1
        entoc(  1 ,j) = 0.5d0*( xfo(  1 ,j-1) + xfo(  1 ,j) )
        entoc(nxpo,j) = 0.5d0*( xfo(nxto,j-1) + xfo(nxto,j) )
      enddo
!$OMP END DO NOWAIT
!$OMP SINGLE
      entoc(  1 ,  1 ) = xfo(  1 ,  1 )
      entoc(nxpo,  1 ) = xfo(nxto,  1 )
      entoc(  1 ,nypo) = xfo(  1 ,nyto)
      entoc(nxpo,nypo) = xfo(nxto,nyto)
!$OMP END SINGLE NOWAIT
#  endif

!$OMP END PARALLEL

*     MONITORING - extra section for cfraoc, centoc
      cfraoc = cfrasm*ocnorm
      centoc = centsm*dxo*dyo

*     Compute area integral of entrainment between
*     layers 1 and 2 to get net diabatic effect
*     --------------------------------------------
      call xintp (xon(1), entoc, nxpo, nypo)
      xon(1) = xon(1)*dxo*dyo
**    write (*,'(a,1p,8d17.8)') '  xon = ',xon

#  ifdef cyclic_ocean
*     Compute boundary integrals of entrainment
*     -----------------------------------------
      ensums = 0.5d0*entoc(1,  1 )
      ensumn = 0.5d0*entoc(1,nypo)
      do i=2,nxpo-1
        ensums = ensums + entoc(i,  1 )
        ensumn = ensumn + entoc(i,nypo)
      enddo
      ensums = ensums + 0.5d0*entoc(nxpo,  1 )
      ensumn = ensumn + 0.5d0*entoc(nxpo,nypo)
      enisoc(1) = dxo*ensums
      eninoc(1) = dxo*ensumn
#  endif

      END SUBROUTINE oml

#endif /* not atmos_only */
c
c***********************************************************************
c
#ifndef atmos_only

      SUBROUTINE omladf (rhs, po)
*
*     Computes (in rhs) the advective and diffusive
*     contributions to the RHS of the evolution equation
*     for the oceanic mixed layer temperature.
*     Uses time-lagged sea-surface temperature.
*     This version includes the coupling coefficient ycexp.
*     po is the pressure in the top oceanic layer (layer 1).
*     Diffusion implemented using dummy points west & east in del2t.
*     Apply no-flux boundary condition to Del-4th diffusion
*     across solid boundaries, equivalent to T'_{nnn} = 0

*     Modules
      USE parameters
      USE occonst
      USE intrfac, ONLY : sst, sstm, tauxo, tauyo, hmoc,
     &                    st2d, st4d, tsbdy, tnbdy
      USE monitor, ONLY : ttmads, vfmads, ttmdfs, ttmadn, vfmadn, ttmdfn

      IMPLICIT NONE

*     Subroutine arguments
      double precision, INTENT(OUT) :: rhs(nxto,nyto)
      double precision, INTENT(IN) :: po(nxpo,nypo)
*
*     Local variables
      integer i,j
      double precision uvgfac,rhf0hm,d2tfac,d4tfac,vfsmsb,tasmsb,
     &                 tdsmsb,vfsmnb,tasmnb,tdsmnb,um,tm,up,tp,
     &                 hxadv,vm,vp,hyadv,del2t(0:nxto+1,nyto)

      uvgfac = ycexp*rdxof0
      rhf0hm = 0.5d0/(fnot*hmoc)
      d2tfac = st2d*dxom2
      d4tfac = st4d*dxom2**2

*     C-grid advection scheme, second order accurate

*     Inner points + meridional boundaries
*     ====================================

!$OMP PARALLEL DEFAULT (NONE)
!$OMP&         PRIVATE (i,j,um,tm,up,tp,hxadv,vm,vp,hyadv)
!$OMP&         SHARED  (uvgfac,po,rhf0hm,tauyo,sst,hdxom1,
!$OMP&                  tauxo,rhs,del2t,sstm,tsbdy,tnbdy)

!$OMP DO SCHEDULE (STATIC)
      do j=2,nyto-1

*       Western boundary (i=1)
*       ----------------
#  ifdef cyclic_ocean
*       Zonally cyclic ocean
        um = -uvgfac*(    po(1,j+1) -    po(1,j) )
     &       +rhf0hm*( tauyo(1,j+1) + tauyo(1,j) )
        tm = sst(1,j) + sst(nxto,j)
        del2t(1,j) =   sstm( 1 ,j-1) + sstm(nxto,j) + sstm(2,j)
     &               + sstm( 1 ,j+1) - 4.0d0*sstm(1,j)
#  else
*       Finite box ocean
*       Advection (no normal mass flux)
        um = 0.0d0
        tm = 0.0d0
*       Diffusion (no normal heat flux)
        del2t(1,j) =   sstm( 1 ,j-1) + sstm(2,j)
     &               + sstm( 1 ,j+1) - 3.0d0*sstm(1,j)
#  endif
        up = -uvgfac*(    po(2,j+1) -    po(2,j) )
     &       +rhf0hm*( tauyo(2,j+1) + tauyo(2,j) )
        tp = sst(1,j) + sst(2,j)
        hxadv = hdxom1*( up*tp - um*tm )
        vm =  uvgfac*(    po(2, j ) -    po(1, j ) )
     &       -rhf0hm*( tauxo(2, j ) + tauxo(1, j ) )
        vp =  uvgfac*(    po(2,j+1) -    po(1,j+1) )
     &       -rhf0hm*( tauxo(2,j+1) + tauxo(1,j+1) )
        hyadv = hdxom1*( vp*( sst(1,j+1)+sst(1, j ) )
     &                  -vm*( sst(1, j )+sst(1,j-1) ) )
        rhs(1,j) = -( hxadv + hyadv )

*       Inner points (standard case)
*       ------------
        do i=2,nxto-1
          um = up
          tm = tp
          up = -uvgfac*(    po(i+1,j+1) -    po(i+1,j) )
     &         +rhf0hm*( tauyo(i+1,j+1) + tauyo(i+1,j) )
          tp = sst(i,j) + sst(i+1,j)
          hxadv = hdxom1*( up*tp - um*tm )
          vm =  uvgfac*(    po(i+1, j ) -    po(i, j ) )
     &         -rhf0hm*( tauxo(i+1, j ) + tauxo(i, j ) )
          vp =  uvgfac*(    po(i+1,j+1) -    po(i,j+1) )
     &         -rhf0hm*( tauxo(i+1,j+1) + tauxo(i,j+1) )
          hyadv = hdxom1*( vp*( sst(i,j+1)+sst(i, j ) )
     &                    -vm*( sst(i, j )+sst(i,j-1) ) )
          rhs(i,j) = -( hxadv + hyadv )
          del2t(i,j) =   sstm( i ,j-1) + sstm(i-1,j) + sstm(i+1,j)
     &                 + sstm( i ,j+1) - 4.0d0*sstm(i,j)
        enddo

*       Eastern boundary (i=nxto)
*       ----------------
        um = up
        tm = tp
#  ifdef cyclic_ocean
*       Zonally cyclic ocean
        up = -uvgfac*(    po(nxto+1,j+1) -    po(nxto+1,j) )
     &       +rhf0hm*( tauyo(nxto+1,j+1) + tauyo(nxto+1,j) )
        tp = sst(1,j) + sst(nxto,j)
        del2t(nxto,j) =   sstm(nxto,j-1) + sstm(nxto-1,j) + sstm(1,j)
     &                  + sstm(nxto,j+1) - 4.0d0*sstm(nxto,j)
#  else
*       Finite box ocean
*       Advection (no normal mass flux)
        up = 0.0d0
        tp = 0.0d0
*       Diffusion (no normal heat flux)
        del2t(nxto,j) =   sstm(nxto,j-1) + sstm(nxto-1,j)
     &                  + sstm(nxto,j+1) - 3.0d0*sstm(nxto,j)
#  endif
        hxadv = hdxom1*( up*tp - um*tm )
        vm =  uvgfac*(    po(nxto+1, j ) -    po(nxto, j ) )
     &       -rhf0hm*( tauxo(nxto+1, j ) + tauxo(nxto, j ) )
        vp =  uvgfac*(    po(nxto+1,j+1) -    po(nxto,j+1) )
     &       -rhf0hm*( tauxo(nxto+1,j+1) + tauxo(nxto,j+1) )
        hyadv = hdxom1*( vp*( sst(nxto,j+1)+sst(nxto, j ) )
     &                  -vm*( sst(nxto, j )+sst(nxto,j-1) ) )
        rhs(nxto,j) = -( hxadv + hyadv )

*       Fill dummy points of del2t array
*       --------------------------------
#  ifdef cyclic_ocean
*       Zonally cyclic ocean
        del2t(   0  ,j) = del2t(nxto,j)
        del2t(nxto+1,j) = del2t(  1 ,j)
#  else
*       Finite box ocean (no diffusive flux)
        del2t(   0  ,j) = del2t(  1 ,j)
        del2t(nxto+1,j) = del2t(nxto,j)
#  endif

      enddo
!$OMP END DO NOWAIT

*     Zonal boundaries, inner points (choice of BCs)
*     ==============================================
*     (N.B. can't do in recurrence form as above because //ising on i)
!$OMP DO SCHEDULE (STATIC)
      do i=2,nxto-1

*       Southern boundary (j=1)
*       -----------------
*       Advection
        um = -uvgfac*(    po( i ,2) -    po( i ,1) )
     &       +rhf0hm*( tauyo( i ,2) + tauyo( i ,1) )
        up = -uvgfac*(    po(i+1,2) -    po(i+1,1) )
     &       +rhf0hm*( tauyo(i+1,2) + tauyo(i+1,1) )
        hxadv = hdxom1*( up*( sst(i+1,1)+sst( i ,1) )
     &                  -um*( sst( i ,1)+sst(i-1,1) ) )
        vp =  uvgfac*(    po(i+1,2) -    po( i ,2) )
     &       -rhf0hm*( tauxo(i+1,2) + tauxo( i ,2) )
        tp = sst(i,1) + sst(i,2)
#  ifdef sb_hflux
*       Advection consistent with an outflow across the
*       southern boundary equal to the Ekman transport,
*       carrying fluid of a specified temperature tsbdy.
*       p contribution to vm vanishes because p is uniform along bdy.
        vm = -rhf0hm*( tauxo(i+1,1) + tauxo(i,1) )
        tm = sst(i,1) + tsbdy
        hyadv = hdxom1*( vp*tp - vm*tm )
*       Diffusion (option of specified S Bdy temperature)
        del2t(i,1) =   sstm(i-1,1) + sstm(i+1,1) + sstm(i,2)
     &               + tsbdy - 4.0d0*sstm(i,1)
#  else
*       Advection (no normal mass flux => vm = 0)
        hyadv = hdxom1*( vp*tp )
*       Diffusion (no normal heat flux)
        del2t(i,1) =   sstm(i-1,1) + sstm(i+1,1)
     &               + sstm(i,2) - 3.0d0*sstm(i,1)
#  endif
        rhs(i,1) = -( hxadv + hyadv )

*       Northern boundary (j=nyto)
*       -----------------
*       Advection
        um = -uvgfac*(    po( i ,nyto+1) -    po( i ,nyto) )
     &       +rhf0hm*( tauyo( i ,nyto+1) + tauyo( i ,nyto) )
        up = -uvgfac*(    po(i+1,nyto+1) -    po(i+1,nyto) )
     &       +rhf0hm*( tauyo(i+1,nyto+1) + tauyo(i+1,nyto) )
        hxadv = hdxom1*( up*( sst(i+1,nyto)+sst( i ,nyto) )
     &                  -um*( sst( i ,nyto)+sst(i-1,nyto) ) )
        vm =  uvgfac*(    po(i+1,nyto) -    po( i ,nyto) )
     &       -rhf0hm*( tauxo(i+1,nyto) + tauxo( i ,nyto) )
        tm = sst(i,nyto-1) + sst(i,nyto)
#  ifdef nb_hflux
*       Advection consistent with an outflow across the
*       northern boundary equal to the Ekman transport,
*       carrying fluid of a specified temperature tnbdy.
*       p contribution to vp vanishes because p is uniform along bdy.
        vp = -rhf0hm*( tauxo(i+1,nyto+1) + tauxo(i,nyto+1) )
        tp = sst(i,nyto) + tnbdy
        hyadv = hdxom1*( vp*tp - vm*tm )
*       Diffusion (option of specified N Bdy temperature)
        del2t(i,nyto) =   sstm(i,nyto-1) + sstm(i-1,nyto) + tnbdy
     &                  + sstm(i+1,nyto) - 4.0d0*sstm(i,nyto)
#  else
*       Advection (no normal mass flux => vp = 0)
        hyadv = hdxom1*( -vm*tm )
*       Diffusion (no normal heat flux)
        del2t(i,nyto) =   sstm(i,nyto-1) + sstm(i-1,nyto)
     &                  + sstm(i+1,nyto) - 3.0d0*sstm(i,nyto)
#  endif
        rhs(i,nyto) = -( hxadv + hyadv )

      enddo
!$OMP END DO NOWAIT

!$OMP END PARALLEL

*     Corner points, in memory order
*     ==============================

*     Southwestern (i=1, j=1)
*     ------------
#  ifdef cyclic_ocean
*     Zonally cyclic ocean
*     Advection
      um = -uvgfac*(    po(1,2) -    po(1,1) )
     &     +rhf0hm*( tauyo(1,2) + tauyo(1,1) )
      tm = sst(1,1) + sst(nxto,1)
#    ifdef sb_hflux
*     Diffusion (option of specified S Bdy temperature)
      del2t(1,1) =   sstm(nxto,1) + sstm(2,1) + sstm(1,2)
     &             + tsbdy - 4.0d0*sstm(1,1)
#    else
*     Diffusion (no normal heat flux)
      del2t(1,1) =   sstm(nxto,1) + sstm(2,1)
     &             + sstm(1,2) - 3.0d0*sstm(1,1)
#    endif
      del2t(nxto+1,1) = del2t(1,1)
#  else
*     Finite box ocean
*     Advection (no normal mass flux)
      um = 0.0d0
      tm = 0.0d0
#    ifdef sb_hflux
*     Diffusion (option of specified S Bdy temperature)
*     Option of specified S Bdy temperature
      del2t(1,1) =   sstm(2,1) + sstm(1,2)
     &             + tsbdy - 3.0d0*sstm(1,1)
#    else
*     Diffusion (no normal heat flux)
      del2t(1,1) = sstm(2,1) + sstm(1,2) - 2.0d0*sstm(1,1)
#    endif
      del2t(0,1) = del2t(1,1)
#  endif
      up = -uvgfac*(    po(2,2) -    po(2,1) )
     &     +rhf0hm*( tauyo(2,2) + tauyo(2,1) )
      tp = sst(1,1) + sst(2,1)
      hxadv = hdxom1*( up*tp - um*tm )
      vp =  uvgfac*(    po(2,2) -    po(1,2) )
     &     -rhf0hm*( tauxo(2,2) + tauxo(1,2) )
      tp = sst(1,1) + sst(1,2)
#  ifdef sb_hflux
*     Advection consistent with an outflow across the
*     southern boundary equal to the Ekman transport,
*     carrying fluid of a specified temperature tsbdy.
*     p contribution to vm vanishes because p is uniform along bdy.
      vm = -rhf0hm*( tauxo(2,1) + tauxo(1,1) )
      tm = sst(1,1) + tsbdy
      hyadv = hdxom1*( vp*tp - vm*tm )
#  else
*     Advection (no normal mass flux => vm = 0)
      hyadv = hdxom1*( vp*tp )
#  endif
      rhs(1,1) = -( hxadv + hyadv )

*     Southeastern (i=nxto, j=1)
*     ------------
      um = -uvgfac*(    po(nxto,2) -    po(nxto,1) )
     &     +rhf0hm*( tauyo(nxto,2) + tauyo(nxto,1) )
      tm = sst(nxto-1,1) + sst(nxto,1)
#  ifdef cyclic_ocean
*     Zonally cyclic ocean
*     Advection
      up = -uvgfac*(    po(nxto+1,2) -    po(nxto+1,1) )
     &     +rhf0hm*( tauyo(nxto+1,2) + tauyo(nxto+1,1) )
      tp = sst(1,1) + sst(nxto,1)
#    ifdef sb_hflux
*     Diffusion (option of specified S Bdy temperature)
      del2t(nxto,1) =   sstm(nxto-1,1) + sstm(1,1) + sstm(nxto,2)
     &                + tsbdy - 4.0d0*sstm(nxto,1)
#    else
*     Diffusion (no normal heat flux)
      del2t(nxto,1) =   sstm(nxto-1,1) + sstm(1,1)
     &                + sstm(nxto,2) - 3.0d0*sstm(nxto,1)
#    endif
      del2t(  0 ,1) = del2t(nxto,1)
#  else
*     Finite box ocean
*     Advection (no normal mass flux)
      up = 0.0d0
      tp = 0.0d0
#    ifdef sb_hflux
*     Diffusion (option of specified S Bdy temperature)
      del2t(nxto,1) =   sstm(nxto-1,1) + sstm(nxto,2)
     &                + tsbdy - 3.0d0*sstm(nxto,1)
#    else
*     Diffusion (no normal heat flux)
      del2t(nxto,1) =   sstm(nxto-1,1)
     &                + sstm(nxto,2) - 2.0d0*sstm(nxto,1)
#    endif
      del2t(nxto+1,1) = del2t(nxto,1)
#  endif
      hxadv = hdxom1*( up*tp - um*tm )
      vp =  uvgfac*(    po(nxto+1,2) -    po(nxto,2) )
     &     -rhf0hm*( tauxo(nxto+1,2) + tauxo(nxto,2) )
      tp = sst(nxto,1) + sst(nxto,2)
#  ifdef sb_hflux
*     Advection consistent with an outflow across the
*     southern boundary equal to the Ekman transport,
*     carrying fluid of a specified temperature tsbdy.
*     p contribution to vm vanishes because p is uniform along bdy.
      vm = -rhf0hm*( tauxo(nxto+1,1) + tauxo(nxto,1) )
      tm = sst(nxto,1) + tsbdy
      hyadv = hdxom1*( vp*tp - vm*tm )
#  else
*     Advection (no normal mass flux => vm = 0)
      hyadv = hdxom1*( vp*tp )
#  endif
      rhs(nxto,1) = -( hxadv + hyadv )

*     Northwestern (i=1, j=nyto)
*     ------------
#  ifdef cyclic_ocean
*     Zonally cyclic ocean
*     Advection
      um = -uvgfac*(    po(1,nyto+1) -    po(1,nyto) )
     &     +rhf0hm*( tauyo(1,nyto+1) + tauyo(1,nyto) )
      tm = sst(1,nyto) + sst(nxto,nyto)
#    ifdef nb_hflux
*     Diffusion (option of specified N Bdy temperature)
      del2t(1,nyto) =   sstm(1,nyto-1) + sstm(nxto,nyto) + tnbdy
     &                + sstm(2,nyto) - 4.0d0*sstm(1,nyto)
#    else
*     Diffusion (no normal heat flux)
      del2t(1,nyto) =   sstm(1,nyto-1) + sstm(nxto,nyto)
     &                + sstm(2,nyto) - 3.0d0*sstm(1,nyto)
#    endif
      del2t(nxto+1,nyto) = del2t(1,nyto)
#  else
*     Finite box ocean
*     Advection (no normal mass flux)
      um = 0.0d0
      tm = 0.0d0
#    ifdef nb_hflux
*     Diffusion (option of specified N Bdy temperature)
      del2t(1,nyto) =   sstm(1,nyto-1) + tnbdy
     &                + sstm(2,nyto) - 3.0d0*sstm(1,nyto)
#    else
*     Diffusion (no normal heat flux)
      del2t(1,nyto) = sstm(1,nyto-1) + sstm(2,nyto) - 2.0d0*sstm(1,nyto)
#    endif
      del2t(0,nyto) = del2t(1,nyto)
#  endif
      up = -uvgfac*(    po(2,nyto+1) -    po(2,nyto) )
     &     +rhf0hm*( tauyo(2,nyto+1) + tauyo(2,nyto) )
      tp = sst(1,nyto) + sst(2,nyto)
      hxadv = hdxom1*( up*tp - um*tm )
      vm =  uvgfac*(    po(2,nyto) -    po(1,nyto) )
     &     -rhf0hm*( tauxo(2,nyto) + tauxo(1,nyto) )
      tm = sst(1,nyto-1) + sst(1,nyto)
#  ifdef nb_hflux
*     Advection consistent with an outflow across the
*     northern boundary equal to the Ekman transport,
*     carrying fluid of a specified temperature tnbdy.
*     p contribution to vp vanishes because p is uniform along bdy.
      vp = -rhf0hm*( tauxo(2,nyto+1) + tauxo(1,nyto+1) )
      tp = sst(1,nyto) + tnbdy
      hyadv = hdxom1*( vp*tp - vm*tm )
#  else
*     Advection (no normal mass flux => vp = 0)
      hyadv = hdxom1*( -vm*tm )
#  endif
      rhs(1,nyto) = -( hxadv + hyadv )

*     Northeastern (i=nxto, j=nyto)
*     ------------
      um = -uvgfac*(    po(nxto,nyto+1) -    po(nxto,nyto) )
     &     +rhf0hm*( tauyo(nxto,nyto+1) + tauyo(nxto,nyto) )
      tm = sst(nxto-1,nyto) + sst(nxto,nyto)
#  ifdef cyclic_ocean
*     Zonally cyclic ocean
*     Advection
      up = -uvgfac*(    po(nxto+1,nyto+1) -    po(nxto+1,nyto) )
     &     +rhf0hm*( tauyo(nxto+1,nyto+1) + tauyo(nxto+1,nyto) )
      tp = sst(1,nyto) + sst(nxto,nyto)
#    ifdef nb_hflux
*     Diffusion (option of specified N Bdy temperature)
      del2t(nxto,nyto) =   sstm(nxto,nyto-1) + sstm(nxto-1,nyto)
     &                   + sstm(1,nyto) - 4.0d0*sstm(nxto,nyto) + tnbdy
#    else
*     Diffusion (no normal heat flux)
      del2t(nxto,nyto) =   sstm(nxto,nyto-1) + sstm(nxto-1,nyto)
     &                   + sstm(1,nyto) - 3.0d0*sstm(nxto,nyto)
#    endif
      del2t(0,nyto) = del2t(nxto,nyto)
#  else
*     Finite box ocean
*     Advection (no normal mass flux)
      up = 0.0d0
      tp = 0.0d0
#    ifdef nb_hflux
*     Diffusion (option of specified N Bdy temperature)
      del2t(nxto,nyto) =   sstm(nxto,nyto-1) + sstm(nxto-1,nyto)
     &                   + tnbdy - 3.0d0*sstm(nxto,nyto)
#    else
*     Diffusion (no normal heat flux)
      del2t(nxto,nyto) =   sstm(nxto,nyto-1) + sstm(nxto-1,nyto)
     &                   - 2.0d0*sstm(nxto,nyto)
#    endif
      del2t(nxto+1,nyto) = del2t(nxto,nyto)
#  endif
      hxadv = hdxom1*( up*tp - um*tm )
      vm =  uvgfac*(    po(nxto+1,nyto) -    po(nxto,nyto) )
     &     -rhf0hm*( tauxo(nxto+1,nyto) + tauxo(nxto,nyto) )
      tm = sst(nxto,nyto) + sst(nxto,nyto-1)
#  ifdef nb_hflux
*     Advection consistent with an outflow across the
*     northern boundary equal to the Ekman transport,
*     carrying fluid of a specified temperature tnbdy.
*     p contribution to vp vanishes because p is uniform along bdy.
      vp = -rhf0hm*( tauxo(nxto+1,nyto+1) + tauxo(nxto,nyto+1) )
      tp = sst(nxto,nyto) + tnbdy
      hyadv = hdxom1*( vp*tp - vm*tm )
#  else
*     (no normal mass flux => vp = 0)
      hyadv = hdxom1*( -vm*tm )
#  endif
      rhs(nxto,nyto) = -( hxadv + hyadv )

*     MONITORING - fluxes from modified boundary conditions
*     =====================================================
      vfsmsb = 0.0d0
      tasmsb = 0.0d0
      tdsmsb = 0.0d0
      vfsmnb = 0.0d0
      tasmnb = 0.0d0
      tdsmnb = 0.0d0
*     Sign convention: fluxes into domain are +ve
#  ifdef sb_hflux
      do i=1,nxto
*       Advection consistent with an outflow across the
*       southern boundary equal to the Ekman transport,
*       carrying fluid of a specified temperature tsbdy.
*       p contribution to vm vanishes because p is uniform along bdy.
        vm = -rhf0hm*( tauxo(i+1,1) + tauxo(i,1) )
        tm = sst(i,1) + tsbdy
        vfsmsb = vfsmsb + vm
        tasmsb = tasmsb + vm*tm
        tdsmsb = tdsmsb - ( sstm(i,1) - tsbdy )
      enddo
#  endif
#  ifdef nb_hflux
      do i=1,nxto
*       Advection consistent with an outflow across the
*       northern boundary equal to the Ekman transport,
*       carrying fluid of a specified temperature tnbdy.
*       p contribution to vp vanishes because p is uniform along bdy.
        vp = -rhf0hm*( tauxo(i+1,nyto+1) + tauxo(i,nyto+1) )
        tp = sst(i,nyto) + tnbdy
        vfsmnb = vfsmnb - vp
        tasmnb = tasmnb - vp*tp
        tdsmnb = tdsmnb + ( tnbdy - sstm(i,nyto) )
      enddo
#  endif
*     Southern bdy
      ttmads = hdxom1*tasmsb/dble( nxto )
      vfmads = vfsmsb/dble( nxto )
      ttmdfs = d2tfac*tdsmsb/dble( nxto )
*     Northern bdy
      ttmadn = hdxom1*tasmnb/dble( nxto )
      vfmadn = vfsmnb/dble( nxto )
      ttmdfn = d2tfac*tdsmnb/dble( nxto )

*     Add Del-sqd and Del-4th terms to temperature evolution term
*     ===========================================================
!$OMP PARALLEL DEFAULT (NONE)
!$OMP&         PRIVATE (i,j)
!$OMP&         SHARED  (rhs,d2tfac,del2t,d4tfac)

*     Inner points in y; complete x range
*     -----------------------------------
!$OMP DO SCHEDULE (STATIC)
      do j=2,nyto-1
        do i=1,nxto
          rhs(i,j) = rhs(i,j) + d2tfac*del2t(i,j) - d4tfac*
     &               (  del2t( i ,j-1) + del2t(i-1,j) + del2t(i+1,j)
     &                + del2t( i ,j+1) - 4.0d0*del2t(i,j) )
        enddo
      enddo
!$OMP END DO NOWAIT

*     Zonal boundaries
*     ----------------
!$OMP DO SCHEDULE (STATIC)
      do i=1,nxto
*       Southern boundary (j=1)
        rhs(i,  1 ) = rhs(i,  1 ) + d2tfac*del2t(i,  1 ) - d4tfac*
     &                 (  del2t(i-1,1) + del2t(i+1,1)
     &                  + del2t(i,2) - 3.0d0*del2t(i,1) )
*       Northern boundary (j=nyto)
        rhs(i,nyto) = rhs(i,nyto) + d2tfac*del2t(i,nyto) - d4tfac*
     &                 (  del2t(i,nyto-1) + del2t(i-1,nyto)
     &                  + del2t(i+1,nyto) - 3.0d0*del2t(i,nyto) )
      enddo
!$OMP END DO NOWAIT

!$OMP END PARALLEL

      END SUBROUTINE omladf

#endif /* not atmos_only */
c
c***********************************************************************
c
      END MODULE omlsubs
c
c***********************************************************************
