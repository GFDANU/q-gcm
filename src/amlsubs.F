c***********************************************************************
c     Q-GCM Version 1.5.0 : last modified 07/06/2013
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
      MODULE amlsubs

*     Contains subroutines aml and amladf which step the atmospheric
*     mixed layer evolution equations (7.16) and (7.17), and update
*     the entrainment between atmospheric layers given by (7.18),
*     with possible convective correction given by (7.19).

*     Modules

      IMPLICIT NONE

      PRIVATE

#ifndef ocean_only
      PUBLIC  :: aml
      PRIVATE :: amladf

      CONTAINS
#endif

c***********************************************************************
*
#ifndef ocean_only

      SUBROUTINE aml
*
*     Timestep atmospheric mixed layer height - equation (7.16)
*     and mixed layer temperature - equation (7.17).
*     Also compute entrainment between atmospheric layers
*     - equation (7.18), plus convective correction (7.19).

*     Modules
      USE parameters, ONLY : nxta, nyta, nla, nxpa, nypa, atnorm
      USE atconst, ONLY : dxa, dya, gpat, tat,
     &                    tdta, xc1ast, xcexp, dtopat
      USE athomog, ONLY : xan, enisat, eninat
      USE atstate, ONLY : pa, pam, entat, wekta
      USE intrfac, ONLY : ast, astm, fnetat, hmixa,
     &                    hmixam, hmadmp, hmat, hmamin
      USE radiate, ONLY : rrcpat, aface, bface, cface, dface
      USE monitor, ONLY : cfraat, centat
      USE intsubs, ONLY : xintp

      IMPLICIT NONE

*     Subroutine arguments
*
*     Local variables
      integer i,j,l
      double precision hmainv,hdrcdt,diabcr,entfac,xbfac,afacdp(nla-1),
     &                 tmrhs(nxta,nyta),hmrhs(nxta,nyta),cfrasm,centsm,
     &                 dhdiab,hnew,dhfix,dtfix,trhtot,astnew,xfaent,
     &                 dtanew,conena,xfa(nxta,nyta),delpm(nla-1),
     &                 adpsum,ensums,ensumn
*
      hmainv = 1.0d0/hmat
      hdrcdt = hmadmp*rrcpat*tdta
      diabcr = tat(1) - 2.0d0*hdrcdt
      entfac = 1.0d0/( tdta*(tat(2)-tat(1)) )

*     Entrainment factors: aface(l), xbfac, cface and dface such that
*     e(1) = Sum(l)[ aface(l)*eta(l) ] + xbfac*etam + cface*aD + dface*aTm'
*          = Sum(l)[ afacdp(l)*dp(l) ] + xbfac*etam + cface*aD + dface*aTm'
*     We are assuming all entrainment is across interface 1.
      xbfac = xcexp*bface
      do l=1,nla-1
        afacdp(l) = aface(l)/gpat(l)
      enddo

*     Initialise rh sides with advective, diffusive
*     and Del-4th terms - first 3 terms in equation
*     (7.16) and first four terms in equation (7.17)
*     ----------------------------------------------
      call amladf (tmrhs, hmrhs, pa(1,1,1))

*     Add forcing term and diabatic effect - the 5th and
*     6th terms in equation (7.17), then timestep ast
*     Also do last term in equation (7.16) and step that
*     --------------------------------------------------
*     MONITORING - extra section for cfraat, centat
      cfrasm = 0.0d0
      centsm = 0.0d0

!$OMP PARALLEL DEFAULT (NONE)
!$OMP&         PRIVATE (i,j,l,delpm,adpsum)
!$OMP&         SHARED  (astm,diabcr,hdrcdt,hmixam,hmat,tat,tdta,hmrhs,
!$OMP&                  hmamin,tmrhs,rrcpat,fnetat,hmainv,wekta,xbfac,
!$OMP&                  dface,xcexp,xc1ast,entfac,hmixa,xfa,cfrasm,
!$OMP&                  centsm,ast,afacdp,cface,pam,entat,dtopat)

!$OMP DO SCHEDULE (STATIC)
!$OMP&   PRIVATE (dhdiab,hnew,dhfix,dtfix,trhtot,
!$OMP&            astnew,xfaent,dtanew,conena)
!$OMP&   REDUCTION(+:cfrasm) REDUCTION(-:centsm)
      do j=1,nyta
        do i=1,nxta
*         Predict new value of hmixa - equation (7.16)
          if ( astm(i,j).le.diabcr ) then
            dhdiab = hdrcdt*( hmixam(i,j) - hmat )
     &                     /( tat(1) - astm(i,j) )
            hnew = hmixam(i,j) + tdta*hmrhs(i,j) - dhdiab
            dhfix = max( hmamin - hnew, 0.0d0 )
**          if ( dhfix.gt.0.0d0 ) write (*,'(a,2i5,1p,4d13.4)')
**   &        '  aml fix  hnew  : ',i,j,dhfix,hmixam(i,j),
**   &                              tdta*hmrhs(i,j),dhdiab
            hnew = hnew + dhfix
            dtfix = dhfix*( tat(1) - astm(i,j) )/hmixam(i,j)
           else
            hnew = hmat
            dtfix = 0.0d0
**          write (*,'(a,2i5,1p,3d13.4)')
**   &        '  aml fix diabcr : ',i,j,astm(i,j),diabcr
          endif

*         Predict new ast - equation (7.17)
          trhtot =  tmrhs(i,j) + rrcpat*fnetat(i,j)/hmixam(i,j)
     &             - hmainv*wekta(i,j)*astm(i,j)
          astnew = astm(i,j) + tdta*trhtot + dtfix

*         Find layer 1/2 entrainment at T points - equation (7.18)
*         Just terms in eta-m and Tm' (naturally at T points) here;
*         other terms (naturally at p points) are added later
          xfaent =  xbfac*( hmixam(i,j) - hmat )
     &            + dface*( xcexp*astm(i,j) + xc1ast(i,j) )

*         Check for convection & if necessary correct layer 1/2
*         entrainment and mixed layer temperature - equation (7.19)
*         dtanew should be >= 0 (stable case)
*         Correction is nonzero only if dtanew < 0
          dtanew = tat(1) - astnew
          conena = entfac*hmixa(i,j)*min( 0.0d0, dtanew )
          xfa(i,j) = xfaent - xcexp*conena
          astnew = astnew + min( 0.0d0, dtanew )

*         MONITORING - extra section for cfraat, centat
          cfrasm = cfrasm + ( 0.5d0 - sign( 0.5d0, dtanew ) )
          centsm = centsm - conena

          astm(i,j) = ast(i,j)
          ast(i,j) = astnew
          hmixam(i,j) = hmixa(i,j)
          hmixa(i,j) = hnew
        enddo
      enddo
!$OMP END DO

*     Average entrainment (computed in xfa at T points) onto p points
*     ---------------------------------------------------------------
*     Internal points
!$OMP DO SCHEDULE (STATIC)
      do j=2,nypa-1
*       Cyclicity condition
        entat(  1 ,j) = 0.25d0*(  xfa(nxta,j-1) + xfa(1,j-1)
     &                          + xfa(nxta, j ) + xfa(1, j ) )
        do i=2,nxpa-1
          entat(i,j) = 0.25d0*(  xfa(i-1,j-1) + xfa(i,j-1)
     &                         + xfa(i-1, j ) + xfa(i, j ) )
        enddo
        entat(nxpa,j) = entat(1,j)
      enddo
!$OMP END DO NOWAIT

*     Do edge and corner values to conserve area integral
!$OMP SINGLE
      entat(  1 ,  1 ) = 0.5d0*( xfa(nxta,  1 ) + xfa(1,  1 ) )
      entat(  1 ,nypa) = 0.5d0*( xfa(nxta,nyta) + xfa(1,nyta) )
      do i=2,nxpa-1
        entat(i,  1 ) = 0.5d0*( xfa(i-1,  1 ) + xfa(i,  1 ) )
        entat(i,nypa) = 0.5d0*( xfa(i-1,nyta) + xfa(i,nyta) )
      enddo
      entat(nxpa,  1 ) = entat(  1 ,  1 )
      entat(nxpa,nypa) = entat(  1 ,nypa)
!$OMP END SINGLE

*     Add eta and topography contributions
*     which evaluate naturally at p points
*     ------------------------------------
!$OMP DO SCHEDULE (STATIC)
      do j=1,nypa
        do i=1,nxpa
          do l=1,nla-1
            delpm(l) = pam(i,j,l) - pam(i,j,l+1)
          enddo
          adpsum = 0.0d0
          do l=1,nla-1
            adpsum = adpsum + afacdp(l)*delpm(l)
          enddo
          entat(i,j) = entat(i,j) + adpsum + cface*dtopat(i,j)
        enddo
      enddo
!$OMP END DO NOWAIT

!$OMP END PARALLEL

*     MONITORING - extra section for cfraat, centat
      cfraat = cfrasm*atnorm
      centat = centsm*dxa*dya

*     Compute area integral of entrainment between layers 1 and 2
      call xintp (xan(1), entat, nxpa, nypa)
      xan(1) = xan(1)*dxa*dya

*     Compute boundary integrals of entrainment
*     -----------------------------------------
      ensums = 0.5d0*entat(1,  1 )
      ensumn = 0.5d0*entat(1,nypa)
      do i=2,nxpa-1
        ensums = ensums + entat(i,  1 )
        ensumn = ensumn + entat(i,nypa)
      enddo
      ensums = ensums + 0.5d0*entat(nxpa,  1 )
      ensumn = ensumn + 0.5d0*entat(nxpa,nypa)
      enisat(1) = dxa*ensums
      eninat(1) = dxa*ensumn

      END SUBROUTINE aml

#endif /* not ocean_only */
c
c***********************************************************************
c
#ifndef ocean_only

      SUBROUTINE amladf (tmrhs, hmrhs, pa)
*
*     Computes (in tmrhs and hmrhs) the advective and diffusive
*     contributions to the RHS of the evolution equations for the
*     atmospheric mixed layer temperature and thickness.
*     Uses time-lagged mixed layer quantities.
*     pa is the pressure in the lowest atmospheric layer (layer 1).
*     Temperature equation contains both Del-sqd and Del-4th terms.
*     Diffusion implemented using dummy points west & east in del2t.
*     Apply no-flux boundary condition to Del-4th diffusion
*     across solid boundaries, equivalent to T'_{nnn} = 0
*     Version with N.E. corner bugfix included
*     (thanks to Tim Leslie for spotting that!).

*     Modules
      USE parameters, ONLY : nxta, nyta, nxpa, nypa
      USE atconst, ONLY : hdxam1, dxam2, rdxaf0
      USE intrfac, ONLY : ast, astm, uekat, vekat, hmixa,
     &                    hmixam, hmat, at2d, at4d, ahmd

      IMPLICIT NONE

*     Subroutine arguments
      double precision, INTENT(OUT) :: tmrhs(nxta,nyta),hmrhs(nxta,nyta)
      double precision, INTENT(IN) :: pa(nxpa,nypa)
*
*     Local variables
      integer i,j
      double precision um,tm,hm,up,tp,hp,xadvt,xadvh,vm,vp,yadvt,
     &                 yadvh,d2tfac,d4tfac,hmdfac,del2t(0:nxta+1,nyta)

      d2tfac = at2d*dxam2
      d4tfac = at4d*dxam2**2
      hmdfac = ahmd*dxam2
*
*     C-grid advection scheme, second order accurate
*     The Ekman components uekat, vekat of the horizontal
*     flow have been previously computed in xforc

*     Inner points + meridional boundaries
*     ------------------------------------
!$OMP PARALLEL DEFAULT (NONE)
!$OMP&         PRIVATE (i,j,um,vm,tm,hm,up,vp,tp,hp,xadvt,xadvh,yadvt,
!$OMP&                  yadvh)
!$OMP&         SHARED  (rdxaf0,uekat,vekat,pa,hdxam1,ast,hmixa,
!$OMP&                  tmrhs,del2t,astm,hmrhs,hmdfac,hmixam,hmat)

!$OMP DO SCHEDULE (STATIC)
      do j=2,nyta-1

*       Western boundary (i=1)
*       ----------------
*       Zonally cyclic atmosphere
        um = -rdxaf0*( pa(1,j+1) - pa(1,j) ) + uekat(1,j)
        tm = ast(1,j) + ast(nxta,j)
        hm = hmixa(1,j) + hmixa(nxta,j)
        up = -rdxaf0*( pa(2,j+1) - pa(2,j) ) + uekat(2,j)
        tp = ast(1,j) + ast(2,j)
        hp = hmixa(1,j) + hmixa(2,j)
        xadvt = hdxam1*( up*tp - um*tm )
        xadvh = hdxam1*( up*hp - um*hm )
        vm =  rdxaf0*( pa(2, j ) - pa(1, j ) ) + vekat(1, j )
        vp =  rdxaf0*( pa(2,j+1) - pa(1,j+1) ) + vekat(1,j+1)
        yadvt = hdxam1*( vp*( ast(1,j+1)+ast(1, j ) )
     &                  -vm*( ast(1, j )+ast(1,j-1) ) )
        yadvh = hdxam1*( vp*( hmixa(1,j+1)+hmixa(1, j ) )
     &                  -vm*( hmixa(1, j )+hmixa(1,j-1) ) )
        tmrhs( 1 ,j) = -( xadvt + yadvt )
        del2t( 1 ,j) =   astm( 1 ,j-1) + astm(nxta,j) + astm(2,j)
     &                 + astm( 1 ,j+1) - 4.0d0*astm(1,j)
        hmrhs( 1 ,j) = -( xadvh + yadvh ) + hmdfac*
     &                (  hmixam( 1 ,j-1) + hmixam(nxta,j) + hmixam(2,j)
     &                 + hmixam( 1 ,j+1) - 4.0d0*hmixam(1,j) )

*       Inner points (standard case)
*       ----------------------------
        do i=2,nxta-1
          um = up
          tm = tp
          hm = hp
          up = -rdxaf0*( pa(i+1,j+1) - pa(i+1,j) ) + uekat(i+1,j)
          tp = ast(i,j) + ast(i+1,j)
          hp = hmixa(i,j) + hmixa(i+1,j)
          xadvt = hdxam1*( up*tp - um*tm )
          xadvh = hdxam1*( up*hp - um*hm )
          vm =  rdxaf0*( pa(i+1, j ) - pa(i, j ) ) + vekat(i, j )
          vp =  rdxaf0*( pa(i+1,j+1) - pa(i,j+1) ) + vekat(i,j+1)
          yadvt = hdxam1*( vp*(ast(i,j+1)+ast(i, j ))
     &                    -vm*(ast(i, j )+ast(i,j-1)) )
          yadvh = hdxam1*( vp*(hmixa(i,j+1)+hmixa(i, j ))
     &                    -vm*(hmixa(i, j )+hmixa(i,j-1)) )
          tmrhs(i,j) = -( xadvt + yadvt )
          del2t(i,j) =   astm( i ,j-1) + astm(i-1,j) + astm(i+1,j)
     &                 + astm( i ,j+1) - 4.0d0*astm(i,j)
          hmrhs(i,j) = -( xadvh + yadvh ) + hmdfac*
     &               (  hmixam( i ,j-1) + hmixam(i-1,j) + hmixam(i+1,j)
     &                + hmixam( i ,j+1) - 4.0d0*hmixam(i,j) )
        enddo

*       Eastern boundary (i=nxta)
*       ----------------
*       Zonally cyclic atmosphere
        um = up
        tm = tp
        hm = hp
        up = -rdxaf0*( pa(nxpa,j+1) - pa(nxpa,j) ) + uekat(nxpa,j)
        tp = ast(1,j) + ast(nxta,j)
        hp = hmixa(1,j) + hmixa(nxta,j)
        xadvt = hdxam1*( up*tp - um*tm )
        xadvh = hdxam1*( up*hp - um*hm )
        vm =  rdxaf0*( pa(nxpa, j ) - pa(nxta, j ) ) + vekat(nxta, j )
        vp =  rdxaf0*( pa(nxpa,j+1) - pa(nxta,j+1) ) + vekat(nxta,j+1)
        yadvt = hdxam1*( vp*( ast(nxta,j+1)+ast(nxta, j ) )
     &                  -vm*( ast(nxta, j )+ast(nxta,j-1) ) )
        yadvh = hdxam1*( vp*( hmixa(nxta,j+1)+hmixa(nxta, j ) )
     &                  -vm*( hmixa(nxta, j )+hmixa(nxta,j-1) ) )
        tmrhs(nxta,j) = -( xadvt + yadvt )
        del2t(nxta,j) =   astm(nxta,j-1) + astm(nxta-1,j) + astm(1,j)
     &                  + astm(nxta,j+1) - 4.0d0*astm(nxta,j)
        hmrhs(nxta,j) = -( xadvh + yadvh ) + hmdfac*
     &             (  hmixam(nxta,j-1) + hmixam(nxta-1,j) + hmixam(1,j)
     &              + hmixam(nxta,j+1) - 4.0d0*hmixam(nxta,j) )

*       Fill dummy points of del2t array
*       --------------------------------
*       Zonally cyclic atmosphere
        del2t(   0  ,j) = del2t(nxta,j)
        del2t(nxta+1,j) = del2t(  1 ,j)

      enddo
!$OMP END DO NOWAIT

*     Zonal boundaries
*     ================
*     (N.B. can't do in recurrence form as above because //ising on i)
*     No normal heat flux in temperature equation
*     Normal flux allowed in thickness equation;
*     thickness set to hmat outside boundary.
*     N.B. P constant along boundaries, so Px vanishes.
*     Inner points (standard case)
*     ----------------------------
!$OMP DO SCHEDULE (STATIC)
      do i=2,nxta-1

*       Southern boundary (j=1)
        um = -rdxaf0*( pa( i ,2) - pa( i ,1) ) + uekat( i ,1)
        up = -rdxaf0*( pa(i+1,2) - pa(i+1,1) ) + uekat(i+1,1)
        xadvt = hdxam1*( up*( ast(i+1,1)+ast( i ,1) )
     &                  -um*( ast( i ,1)+ast(i-1,1) ) )
        xadvh = hdxam1*( up*( hmixa(i+1,1)+hmixa( i ,1) )
     &                  -um*( hmixa( i ,1)+hmixa(i-1,1) ) )
        vm =  vekat(i,1)
        vp =  rdxaf0*( pa(i+1,2) - pa(i,2) ) + vekat(i,2)
        yadvt = hdxam1*vp*(ast(i,2)+ast(i, 1 ))
        yadvh = hdxam1*( vp*(hmixa(i,2)+hmixa(i,1))
     &                  -vm*(hmixa(i,1)+hmat) )
        tmrhs(i,  1  ) = -( xadvt + yadvt )
        del2t(i,  1  ) =   astm(i-1,1) + astm(i+1,1)
     &                   + astm(i,2) - 3.0d0*astm(i,1)
        hmrhs(i,  1  ) = -( xadvh + yadvh ) + hmdfac*
     &                  (  hmat + hmixam(i-1,1) + hmixam(i+1,1)
     &                   + hmixam(i,2) - 4.0d0*hmixam(i,1) )

*       Northern boundary (j=nyta)
        um = -rdxaf0*( pa( i ,nypa) - pa( i ,nyta) ) + uekat( i ,nyta)
        up = -rdxaf0*( pa(i+1,nypa) - pa(i+1,nyta) ) + uekat(i+1,nyta)
        xadvt = hdxam1*( up*( ast(i+1,nyta) + ast( i ,nyta) )
     &                  -um*( ast( i ,nyta) + ast(i-1,nyta) ) )
        xadvh = hdxam1*( up*( hmixa(i+1,nyta)+hmixa( i ,nyta) )
     &                  -um*( hmixa( i ,nyta)+hmixa(i-1,nyta) ) )
        vm =  rdxaf0*( pa(i+1,nyta) - pa(i,nyta) ) + vekat(i,nyta)
        vp =  vekat(i,nyta+1)
        yadvt = hdxam1*( -vm*( ast(i, nyta )+ast(i,nyta-1) ) )
        yadvh = hdxam1*( vp*( hmat+hmixa(i,nyta) )
     &                  -vm*( hmixa(i,nyta)+hmixa(i,nyta-1) ) )
        tmrhs(i,nyta) = -( xadvt + yadvt )
        del2t(i,nyta) =   astm(i,nyta-1) + astm(i-1,nyta)
     &                  + astm(i+1,nyta) - 3.0d0*astm(i,nyta)
        hmrhs(i,nyta) = -( xadvh + yadvh ) + hmdfac*
     &                 (  hmixam( i ,nyta-1) + hmixam(i-1,nyta)
     &                  + hmixam(i+1,nyta)
     &                  + hmat - 4.0d0*hmixam(i,nyta) )

      enddo
!$OMP END DO NOWAIT

!$OMP END PARALLEL

*     Corner points, in memory order
*     ==============================
*     No normal flux for temperature
*     Normal flux allowed for thickness;
*     thickness fixed at hmat outside boundary.
*     Periodic in other direction.

*     Southwestern (i=1, j=1)
*     ------------
      um = -rdxaf0*( pa(1,2) - pa(1,1) ) + uekat(1,1)
      up = -rdxaf0*( pa(2,2) - pa(2,1) ) + uekat(2,1)
      xadvt = hdxam1*( up*( ast(2,1)+ast(  1 ,1) )
     &                -um*( ast(1,1)+ast(nxta,1) ) )
      xadvh = hdxam1*( up*( hmixa(2,1)+hmixa(  1 ,1) )
     &                -um*( hmixa(1,1)+hmixa(nxta,1) ) )
      vm =  vekat(1,1)
      vp =  rdxaf0*( pa(2,2) - pa(1,2) ) + vekat(1,2)
      yadvt = hdxam1*vp*( ast(1,2)+ast(1,1) )
      yadvh = hdxam1*( vp*( hmixa(1,2)+hmixa(1,1) )
     &                -vm*( hmixa(1,1)+hmat ) )
      tmrhs( 1 , 1 ) = -( xadvt + yadvt )
      del2t( 1 , 1 ) =   astm(nxta,1) + astm(2,1)
     &                 + astm(1,2) - 3.0d0*astm(1,1)
      del2t(nxta+1,1) = del2t(1,1)
      hmrhs( 1 , 1 ) = -( xadvh + yadvh ) + hmdfac*
     &                 (  hmat + hmixam(nxta,1) + hmixam(2,1)
     &                  + hmixam(1,2) - 4.0d0*hmixam(1,1) )

*     Southeastern (i=nxta, j=1)
*     ------------
      um = -rdxaf0*( pa(nxta,2) - pa(nxta,1) ) + uekat(nxta,1)
      up = -rdxaf0*( pa(nxpa,2) - pa(nxpa,1) ) + uekat(nxpa,1)
      xadvt = hdxam1*( up*( ast(  1 ,1)+ast( nxta ,1) )
     &                -um*( ast(nxta,1)+ast(nxta-1,1) ) )
      xadvh = hdxam1*( up*( hmixa(  1 ,1)+hmixa( nxta ,1) )
     &                -um*( hmixa(nxta,1)+hmixa(nxta-1,1) ) )
      vm =  vekat(nxta,1)
      vp =  rdxaf0*( pa(nxpa,2) - pa(nxta,2) ) + vekat(nxta,2)
      yadvt = hdxam1*vp*( ast(nxta,2)+ast(nxta, 1 ) )
      yadvh = hdxam1*( vp*(hmixa(nxta,2)+hmixa(nxta,1))
     &                -vm*(hmixa(nxta,1)+hmat) )
      tmrhs(nxta,1) = -( xadvt + yadvt )
      del2t(nxta,1) =   astm(nxta-1,1) + astm(1,1)
     &                + astm(nxta,2) - 3.0d0*astm(nxta,1)
      del2t(0,1) = del2t(nxta,1)
      hmrhs(nxta,1) = -( xadvh + yadvh ) + hmdfac*
     &               (  hmat + hmixam(nxta-1,1) + hmixam(1,1)
     &                + hmixam(nxta,2) - 4.0d0*hmixam(nxta,1) )

*     Northwestern (i=1, j=nyta)
*     ------------
      um = -rdxaf0*( pa(1,nypa) - pa(1,nyta) ) + uekat(1,nyta)
      up = -rdxaf0*( pa(2,nypa) - pa(2,nyta) ) + uekat(2,nyta)
      xadvt = hdxam1*( up*( ast( 2 ,nyta)+ast(  1 ,nyta) )
     &                -um*( ast( 1 ,nyta)+ast(nxta,nyta) ) )
      xadvh = hdxam1*( up*( hmixa( 2 ,nyta)+hmixa(  1 ,nyta) )
     &                -um*( hmixa( 1 ,nyta)+hmixa(nxta,nyta) ) )
      vm =  rdxaf0*( pa(2,nyta) - pa(1,nyta) ) + vekat(1,nyta)
      vp =  vekat(1,nyta+1)
      yadvt = hdxam1*(-vm*( ast(1, nyta)+ast(1,nyta-1) ) )
      yadvh = hdxam1*( vp*( hmat+hmixa(1,nyta) )
     &                -vm*( hmixa(1,nyta)+hmixa(1,nyta-1) ) )
      tmrhs( 1 ,nyta) = -( xadvt + yadvt )
      del2t( 1 ,nyta) =   astm(1,nyta-1) + astm(nxta,nyta)
     &                  + astm(2,nyta) - 3.0d0*astm(1,nyta)
      del2t(nxta+1,nyta) = del2t( 1 ,nyta)
      hmrhs( 1 ,nyta) = -( xadvh + yadvh ) + hmdfac*
     &                  (  hmixam(1,nyta-1) + hmixam(nxta,nyta)
     &                   + hmixam(2,nyta)
     &                   + hmat - 4.0d0*hmixam(1,nyta) )

*     Northeastern (i=nxta, j=nyta)
*     ------------
      um = -rdxaf0*( pa(nxta,nypa) - pa(nxta,nyta) ) + uekat(nxta,nyta)
      up = -rdxaf0*( pa(nxpa,nypa) - pa(nxpa,nyta) ) + uekat(nxpa,nyta)
      xadvt = hdxam1*( up*( ast(  1 ,nyta)+ast( nxta ,nyta) )
     &                -um*( ast(nxta,nyta)+ast(nxta-1,nyta) ) )
      xadvh = hdxam1*( up*( hmixa(  1 ,nyta)+hmixa( nxta ,nyta) )
     &                -um*( hmixa(nxta,nyta)+hmixa(nxta-1,nyta) ) )
      vm =  rdxaf0*( pa(nxpa,nyta) - pa(nxta,nyta) ) + vekat(nxta,nyta)
      vp =  vekat(nxta,nyta+1 )
      yadvt = hdxam1*( -vm*( ast(nxta,nyta)+ast(nxta,nyta-1) ) )
      yadvh = hdxam1*( vp*( hmat+hmixa(nxta,nyta) )
     &                -vm*( hmixa(nxta,nyta)+hmixa(nxta,nyta-1)) )
      tmrhs(nxta,nyta) = -( xadvt + yadvt )
      del2t(nxta,nyta) =   astm(nxta,nyta-1) + astm(nxta-1,nyta)
     &                   + astm(1,nyta) - 3.0d0*astm(nxta,nyta)
      del2t(0,nyta) = del2t(nxta,nyta)
      hmrhs(nxta,nyta) = -( xadvh + yadvh ) + hmdfac*
     &                   (  hmixam(nxta,nyta-1) + hmixam(nxta-1,nyta)
     &                    + hmixam(1,nyta)
     &                    + hmat - 4.0d0*hmixam(nxta,nyta) )

*     Add Del-sqd and Del-4th terms to temperature evolution term
*     ===========================================================

!$OMP PARALLEL DEFAULT (NONE)
!$OMP&         PRIVATE (i,j)
!$OMP&         SHARED  (tmrhs,d2tfac,del2t,d4tfac)

*     Inner points in y; complete x range
*     -----------------------------------
!$OMP DO SCHEDULE (STATIC)
      do j=2,nyta-1
        do i=1,nxta
          tmrhs(i,j) = tmrhs( i ,j) + d2tfac*del2t( i ,j) - d4tfac*
     &                 (  del2t( i ,j-1) + del2t(i-1,j) + del2t(i+1,j)
     &                  + del2t( i ,j+1) - 4.0d0*del2t(i,j) )
        enddo
      enddo
!$OMP END DO NOWAIT

*     Zonal boundaries
*     ----------------
!$OMP DO SCHEDULE (STATIC)
      do i=1,nxta
*       Southern boundary (j=1)
        tmrhs(i,  1 ) = tmrhs(i,  1 ) + d2tfac*del2t(i,  1 )
     &                  - d4tfac*(  del2t(i-1,1) + del2t(i+1,1)
     &                            + del2t(i,2) - 3.0d0*del2t(i,1) )
*       Northern boundary (j=nyta)
        tmrhs(i,nyta) = tmrhs(i,nyta) + d2tfac*del2t(i,nyta)
     &                  - d4tfac*(  del2t(i,nyta-1) + del2t(i-1,nyta)
     &                     + del2t(i+1,nyta) - 3.0d0*del2t(i,nyta) )
      enddo
!$OMP END DO NOWAIT

!$OMP END PARALLEL

      END SUBROUTINE amladf

#endif /* not ocean_only */
c
c***********************************************************************
c
      END MODULE amlsubs
c
c***********************************************************************
