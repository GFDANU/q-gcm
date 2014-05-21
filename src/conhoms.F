c***********************************************************************
c     Q-GCM Version 1.5.0 : last modified 20/08/2012
c***********************************************************************
c
c     Copyright 2012 Jeff Blundell, Andy Hogg and Bill Dewar.
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
      MODULE conhoms

*     Contains subroutines constr and homsol. constr computes the
*     appropriate mass and momentum constraint variables for both
*     atmosphere and ocean. homsol computes the homogeneous solutions
*     for both atmosphere and ocean, which are added to the pressure
*     fields inverted from vorticity to ensure that the full
*     solution satisfies the mass and momentum contraints.

*     Modules

      IMPLICIT NONE

      PRIVATE

      PUBLIC :: constr, homsol

      CONTAINS

c***********************************************************************
c
      SUBROUTINE constr

*     Computes values of the mass and momentum constraint variables
*     from the pressure fields (current and lagged), for both ocean
*     and atmosphere. Computes the ocean constraint variables only
*     appropriate for either box or zonally cyclic configurations.

*     Modules
      USE parameters
#ifndef ocean_only
      USE atconst
      USE athomog
      USE atstate
#endif
#ifndef atmos_only
      USE occonst
      USE ochomog
      USE ocstate
#endif
      USE intsubs, ONLY : xintp
      USE monitor, ONLY : et2moc, et2mat

      IMPLICIT NONE

*     Subroutine arguments
*
*     Local variables
      integer i,j,k
      double precision apsp,apnp,aps,apn
#ifndef ocean_only
      double precision wka(nxpa,nypa,nla),
     &                 apins(nla),apinn(nla),apinsp(nla),apinnp(nla)
#endif
#ifndef atmos_only
      double precision wko(nxpo,nypo,nlo)
#  ifdef cyclic_ocean
      double precision opins(nlo),opinn(nlo),opinsp(nlo),opinnp(nlo)
#  endif
#endif

      print *,' '
      write (*,*) ' Constraint values:'
      write (*,*) ' ------------------'

#ifndef atmos_only
*     Ocean
*     =====
      do k=1,nlo-1
!$OMP   PARALLEL DO DEFAULT (NONE)
!$OMP&           PRIVATE (i,j)
!$OMP&           SHARED  (k,wko,pom,po)
!$OMP&           SCHEDULE (STATIC)
*       Choose sign of dpioc so that +ve dpioc -> +ve eta
        do j=1,nypo
          do i=1,nxpo
            wko(i,j,1) = pom(i,j,k+1) - pom(i,j,k)
            wko(i,j,2) = po(i,j,k+1)  - po(i,j,k)
          enddo
        enddo
!$OMP   END PARALLEL DO
        call xintp (dpiocp(k), wko(1,1,1), nxpo, nypo)
        call xintp (dpioc(k),  wko(1,1,2), nxpo, nypo)
        dpiocp(k) = dpiocp(k)*dxo*dyo
        dpioc(k)  = dpioc(k)*dxo*dyo
*       Compute eta^2 integral for energy diagnostics
!$OMP   PARALLEL DO DEFAULT (NONE)
!$OMP&           PRIVATE (i,j)
!$OMP&           SHARED  (k,wko,po)
!$OMP&           SCHEDULE (STATIC)
        do j=1,nypo
          do i=1,nxpo
            wko(i,j,2) = ( po(i,j,k+1) - po(i,j,k) )**2
          enddo
        enddo
!$OMP   END PARALLEL DO
        call xintp (et2moc(k), wko(1,1,2), nxpo, nypo)
        et2moc(k) = et2moc(k)*ocnorm/gpoc(k)**2
      enddo
      write (*,240) '  dpiocp = ',dpiocp
      write (*,240) '  dpioc  = ',dpioc

#  ifdef cyclic_ocean
*     Line integrate p and dp/dy for all
*     layers along oceanic zonal boundaries
*     -------------------------------------
      do k=1,nlo
*       Line integrals of p
        opinsp(k) = 0.5d0*pom(1,  1 ,k)
        opinnp(k) = 0.5d0*pom(1,nypo,k)
        opins(k)  = 0.5d0*po (1,  1 ,k)
        opinn(k)  = 0.5d0*po (1,nypo,k)
*       Line integrals of dp/dy
        ocncsp(k) = 0.5d0*( pom(1,  2 ,k) - pom(1,   1  ,k) )
        ocncnp(k) = 0.5d0*( pom(1,nypo,k) - pom(1,nypo-1,k) )
        ocncs(k)  = 0.5d0*( po(1,  2 ,k) - po(1,   1  ,k) )
        ocncn(k)  = 0.5d0*( po(1,nypo,k) - po(1,nypo-1,k) )
        do i=2,nxpo-1
          opinsp(k) = opinsp(k) + pom(i,  1 ,k)
          opinnp(k) = opinnp(k) + pom(i,nypo,k)
          opins(k)  = opins(k)  + po(i,  1 ,k)
          opinn(k)  = opinn(k)  + po(i,nypo,k)
          ocncsp(k) = ocncsp(k) + ( pom(i,  2 ,k) - pom(i,   1  ,k) )
          ocncnp(k) = ocncnp(k) + ( pom(i,nypo,k) - pom(i,nypo-1,k) )
          ocncs(k)  = ocncs(k)  + ( po(i,  2 ,k) - po(i,   1  ,k) )
          ocncn(k)  = ocncn(k)  + ( po(i,nypo,k) - po(i,nypo-1,k) )
        enddo
        opinsp(k) = opinsp(k) + 0.5d0*pom(nxpo,  1 ,k)
        opinnp(k) = opinnp(k) + 0.5d0*pom(nxpo,nypo,k)
        opins(k)  = opins(k)  + 0.5d0* po(nxpo,  1 ,k)
        opinn(k)  = opinn(k)  + 0.5d0* po(nxpo,nypo,k)
        ocncsp(k) = ocncsp(k) + 0.5d0*(  pom(nxpo,   2  ,k)
     &                                 - pom(nxpo,   1  ,k) )
        ocncnp(k) = ocncnp(k) + 0.5d0*(  pom(nxpo, nypo ,k)
     &                                 - pom(nxpo,nypo-1,k) )
        ocncs(k)  = ocncs(k)  + 0.5d0*(   po(nxpo,   2  ,k)
     &                                 -  po(nxpo,   1  ,k) )
        ocncn(k)  = ocncn(k)  + 0.5d0*(   po(nxpo, nypo ,k)
     &                                 -  po(nxpo,nypo-1,k) )
        ocncsp(k) = ocncsp(k)*(dxo/dyo)
        ocncnp(k) = ocncnp(k)*(dxo/dyo)
        ocncs(k)  = ocncs(k)*(dxo/dyo)
        ocncn(k)  = ocncn(k)*(dxo/dyo)
        opinsp(k) = dxo*opinsp(k)
        opinnp(k) = dxo*opinnp(k)
        opins(k)  = dxo*opins(k)
        opinn(k)  = dxo*opinn(k)
      enddo

*     Add pressure integral contributions to
*     get full momentum constraint quantities
*     ---------------------------------------
      do k=1,nlo
        apsp = 0.0d0
        apnp = 0.0d0
        aps  = 0.0d0
        apn  = 0.0d0
        do j=1,nlo
          apsp = apsp + amatoc(k,j)*opinsp(j)
          apnp = apnp + amatoc(k,j)*opinnp(j)
          aps  = aps  + amatoc(k,j)*opins(j)
          apn  = apn  + amatoc(k,j)*opinn(j)
        enddo
*       Change sign of southern derivative terms
        ocncsp(k) = -ocncsp(k) + 0.5d0*dyo*fnot*fnot*apsp
        ocncnp(k) =  ocncnp(k) + 0.5d0*dyo*fnot*fnot*apnp
        ocncs(k)  = -ocncs(k)  + 0.5d0*dyo*fnot*fnot*aps
        ocncn(k)  =  ocncn(k)  + 0.5d0*dyo*fnot*fnot*apn
      enddo
      print *,' '
      write (*,240) '  ocncsp = ',ocncsp
      write (*,240) '  ocncs  = ',ocncs
      write (*,240) '  ocncnp = ',ocncnp
      write (*,240) '  ocncn  = ',ocncn
#  endif
#endif

#ifndef ocean_only
*     Atmosphere
*     ==========
      do k=1,nla-1
!$OMP   PARALLEL DO DEFAULT (NONE)
!$OMP&           PRIVATE (i,j)
!$OMP&           SHARED  (k,wka,pam,pa)
!$OMP&           SCHEDULE (STATIC)
        do j=1,nypa
          do i=1,nxpa
            wka(i,j,1) = pam(i,j,k) - pam(i,j,k+1)
            wka(i,j,2) = pa(i,j,k)  - pa(i,j,k+1)
          enddo
        enddo
!$OMP   END PARALLEL DO
        call xintp (dpiatp(k), wka(1,1,1), nxpa, nypa)
        call xintp (dpiat(k),  wka(1,1,2), nxpa, nypa)
        dpiatp(k) = dpiatp(k)*dxa*dya
        dpiat(k)  = dpiat(k)*dxa*dya
*       Compute eta^2 integral for energy diagnostics
!$OMP   PARALLEL DO DEFAULT (NONE)
!$OMP&           PRIVATE (i,j)
!$OMP&           SHARED  (k,wka,pa)
!$OMP&           SCHEDULE (STATIC)
        do j=1,nypa
          do i=1,nxpa
            wka(i,j,2) = ( pa(i,j,k) - pa(i,j,k+1) )**2
          enddo
        enddo
!$OMP   END PARALLEL DO
        call xintp (et2mat(k), wka(1,1,2), nxpa, nypa)
        et2mat(k) = et2mat(k)*atnorm/gpat(k)**2
      enddo
      print *,' '
      write (*,240) '  dpiatp = ',dpiatp
      write (*,240) '  dpiat  = ',dpiat

*     Line integrate p and dp/dy for all
*     layers along atmospheric zonal boundaries.
*     ------------------------------------------
      do k=1,nla
*       Line integrals of p
        apinsp(k) = 0.5d0*pam(1, 1 ,k)
        apinnp(k) = 0.5d0*pam(1,nypa,k)
        apins(k)  = 0.5d0*pa (1, 1 ,k)
        apinn(k)  = 0.5d0*pa (1,nypa,k)
*       Line integrals of dp/dy
        atmcsp(k) = 0.5d0*( pam(1, 2 ,k) - pam(1,  1  ,k) )
        atmcnp(k) = 0.5d0*( pam(1,nypa,k) - pam(1,nypa-1,k) )
        atmcs(k)  = 0.5d0*( pa(1, 2 ,k) - pa(1,  1  ,k) )
        atmcn(k)  = 0.5d0*( pa(1,nypa,k) - pa(1,nypa-1,k) )
        do i=2,nxpa-1
          apinsp(k) = apinsp(k) + pam(i, 1 ,k)
          apinnp(k) = apinnp(k) + pam(i,nypa,k)
          apins(k)  = apins(k)  + pa(i, 1 ,k)
          apinn(k)  = apinn(k)  + pa(i,nypa,k)
          atmcsp(k) = atmcsp(k) + ( pam(i, 2 ,k) - pam(i,  1  ,k) )
          atmcnp(k) = atmcnp(k) + ( pam(i,nypa,k) - pam(i,nypa-1,k) )
          atmcs(k)  = atmcs(k)  + ( pa(i, 2 ,k) - pa(i,  1  ,k) )
          atmcn(k)  = atmcn(k)  + ( pa(i,nypa,k) - pa(i,nypa-1,k) )
        enddo
        apinsp(k) = apinsp(k) + 0.5d0*pam(nxpa, 1 ,k)
        apinnp(k) = apinnp(k) + 0.5d0*pam(nxpa,nypa,k)
        apins(k)  = apins(k)  + 0.5d0* pa(nxpa, 1 ,k)
        apinn(k)  = apinn(k)  + 0.5d0* pa(nxpa,nypa,k)
        atmcsp(k) = atmcsp(k) + 0.5d0*(  pam(nxpa,   2  ,k)
     &                                 - pam(nxpa,   1  ,k) )
        atmcnp(k) = atmcnp(k) + 0.5d0*(  pam(nxpa, nypa ,k)
     &                                 - pam(nxpa,nypa-1,k) )
        atmcs(k)  = atmcs(k)  + 0.5d0*(   pa(nxpa,   2  ,k)
     &                                 -  pa(nxpa,   1  ,k) )
        atmcn(k)  = atmcn(k)  + 0.5d0*(   pa(nxpa, nypa ,k)
     &                                 -  pa(nxpa,nypa-1,k) )
        atmcsp(k) = atmcsp(k)*(dxa/dya)
        atmcnp(k) = atmcnp(k)*(dxa/dya)
        atmcs(k)  = atmcs(k)*(dxa/dya)
        atmcn(k)  = atmcn(k)*(dxa/dya)
        apinsp(k) = dxa*apinsp(k)
        apinnp(k) = dxa*apinnp(k)
        apins(k)  = dxa*apins(k)
        apinn(k)  = dxa*apinn(k)
      enddo

*     Add pressure integral contributions to
*     get full momentum constraint quantities
*     ---------------------------------------
      do k=1,nla
        apsp = 0.0d0
        apnp = 0.0d0
        aps  = 0.0d0
        apn  = 0.0d0
        do j=1,nla
          apsp = apsp + amatat(k,j)*apinsp(j)
          apnp = apnp + amatat(k,j)*apinnp(j)
          aps  = aps  + amatat(k,j)*apins(j)
          apn  = apn  + amatat(k,j)*apinn(j)
        enddo
*       Change sign of southern derivative terms
        atmcsp(k) = -atmcsp(k) + 0.5d0*dya*fnot*fnot*apsp
        atmcnp(k) =  atmcnp(k) + 0.5d0*dya*fnot*fnot*apnp
        atmcs(k)  = -atmcs(k)  + 0.5d0*dya*fnot*fnot*aps
        atmcn(k)  =  atmcn(k)  + 0.5d0*dya*fnot*fnot*apn
      enddo
      print *,' '
      write (*,240) '  atmcsp = ',atmcsp
      write (*,240) '  atmcs  = ',atmcs
      write (*,240) '  atmcnp = ',atmcnp
      write (*,240) '  atmcn  = ',atmcn
#endif

  240 format(a,1p,9d21.13)

      END SUBROUTINE constr
c
c***********************************************************************
c
      SUBROUTINE homsol
*
*     Computes homogeneous baroclinic modal solutions for both
*     ocean and atmosphere. These are stored in modules ochomog
*     and athomog respectively, and used in subroutines ocinvq
*     and atinvq to satisfy mass and momentum constraints.

*     Modules
      USE parameters
#ifndef ocean_only
      USE atconst
      USE athomog
      USE atisubs
#endif
#ifndef atmos_only
      USE occonst
      USE ochomog
      USE ocisubs
#endif
      USE intsubs, ONLY : xintp

      IMPLICIT NONE

*     Subroutine arguments
*
*     Common parameters and variables
*
*     Local parameters
      logical prhsoc,pheloc,psoloc,prhsat,phelat,psolat
      parameter ( prhsoc=.false., pheloc=.false., psoloc=.false.,
     &            prhsat=.false., phelat=.false., psolat=.false. )
*
*     Local variables
      integer i,j,k,m,info
      double precision pch1yn,pch2yn,pch1ys,pch2ys,
     &                 pchdet,aipch1,aipch2
#ifndef ocean_only
      double precision bat(nxta),atwk1(nxpa,nypa),atwk2(nxpa,nypa)
#endif
#ifndef atmos_only
      double precision boc(nxto)
#  ifdef cyclic_ocean
      double precision ocwk1(nxpo,nypo),ocwk2(nxpo,nypo)
#  endif
#endif

      print *,' '
      print *,' Homogeneous solutions:'
      print *,' ----------------------'

#ifndef atmos_only
*     Compute homogeneous oceanic solutions
*     =====================================
#  ifdef cyclic_ocean
*     Zonally cyclic ocean
*     --------------------
*     Barotropic mode: homogeneous solutions would need
*     to satisfy d2p/dy2 = 0, but the coefficient of
*     integral of p vanishes, so only the homogeneous
*     component with a derivative can be inferred.
*     W.L.O.G. we choose this to be 1 on the Southern
*     boundary, and to vanish on the Northern boundary.
*     We then apply the constraint at the Southern boundary
      do j=1,nypo
        pbhoc(j) = dble(nypo-j)/dble(nypo-1)
      enddo
*     hbs = integral of y-derivative along Southern bdy
      hbsioc = ylo/xlo
      aipbho = 0.5d0*xlo*ylo
      print *,' '
      write (*,*) ' Oceanic homogeneous barotropic solution:'
      write (*,240) '  aipbho         = ',aipbho
      write (*,240) '  hbsioc         = ',hbsioc
*     Check structure of homogeneous ocean barotropic solution
      if ( psoloc ) then
        print *,' '
        write (*,*) ' pbhoc:'
        do j=min(241,nypo),1,-1
          write (*,250) j,pbhoc(j)
        enddo
      endif

*     Now compute baroclinic solutions for use in ocinvq
*     Solutions are functions of y only, but use the usual
*     2-D Helmholtz solver for simplicity and consistency

      print *,' '
      write (*,*) ' Oceanic homogeneous (baroclinic) solutions:'
      do m=1,nlo-1
        print *,' '
        write (*,'(a,i2)') '  Mode: ',m
        do i=1,nxto
          boc(i) = bd2oc(i) - rdm2oc(m+1)
        enddo

*       pch = L(y) + rdm2*sol0, where L(y) is linear in y,
*       and sol0 satisfies Del-sqd(sol0) - rdm2*sol0 = L(y)
*       with the usual solid boundary condition p = 0.
*       For pch1oc, L(y) = 1.0 on S bdy; = 0.0 on N bdy
*       For pch2oc, L(y) = 0.0 on S bdy; = 1.0 on N bdy
*       Specify RHSs
!$OMP   PARALLEL DO DEFAULT (NONE)
!$OMP&           PRIVATE (i,j)
!$OMP&           SHARED  (m,pch1oc,pch2oc,ylo,ypo,ocwk1,ocwk2)
!$OMP&           SCHEDULE (STATIC)
        do j=1,nypo
          pch1oc(j,m) = ( ypo(nypo)-ypo(j) )/ylo
          pch2oc(j,m) = ( ypo(  j )-ypo(1) )/ylo
          do i=1,nxpo
            ocwk1(i,j) = pch1oc(j,m)
            ocwk2(i,j) = pch2oc(j,m)
          enddo
        enddo
!$OMP   END PARALLEL DO
*       Check r.h.s. for homogeneous ocean baroclinic solutions
        if ( prhsoc ) then
          print *,' '
          write (*,*) ' r.h.s. for pch1oc:'
          do j=min(241,nypo),1,-1
            write (*,260) j,(ocwk1(i,j),i=1,3),
     &                    (ocwk1(i,j),i=nxpo-1,nxpo)
          enddo
          print *,' '
          write (*,*) ' r.h.s. for pch2oc:'
          do j=min(241,nypo),1,-1
            write (*,260) j,(ocwk2(i,j),i=1,3),
     &                    (ocwk2(i,j),i=nxpo-1,nxpo)
          enddo
        endif

*       Invert these RHSs for baroclinic homog. solutions (sol0 above)
        call hscyoc (ocwk1, boc)
        call hscyoc (ocwk2, boc)
*       Check Helmholtz component of homogeneous ocean baroclinic solutions
        if ( pheloc ) then
          print *,' '
          write (*,*) ' Helmholtz for pch1oc:'
          do j=min(241,nypo),1,-1
            write (*,260) j,(ocwk1(i,j),i=1,3),
     &                    (ocwk1(i,j),i=nxpo-1,nxpo)
          enddo
          print *,' '
          write (*,*) ' Helmholtz for pch2oc:'
          do j=min(241,nypo),1,-1
            write (*,260) j,(ocwk2(i,j),i=1,3),
     &                    (ocwk2(i,j),i=nxpo-1,nxpo)
          enddo
        endif

*       Add Helmholtz solution to L(y) to get full solutions
*       Solutions in ocwk1, ocwk2 are functions of y only, i.e.
*       independent of i, so just save solution for one i value
!$OMP   PARALLEL DO DEFAULT (NONE)
!$OMP&           PRIVATE (i,j)
!$OMP&           SHARED  (m,pch1oc,pch2oc,ocwk1,ocwk2,rdm2oc)
!$OMP&           SCHEDULE (STATIC)
        do j=1,nypo
          do i=1,nxpo
            ocwk1(i,j) = pch1oc(j,m) + rdm2oc(m+1)*ocwk1(i,j)
            ocwk2(i,j) = pch2oc(j,m) + rdm2oc(m+1)*ocwk2(i,j)
          enddo
          pch1oc(j,m) = ocwk1(1,j)
          pch2oc(j,m) = ocwk2(1,j)
        enddo
!$OMP   END PARALLEL DO
*       Check structure of homogeneous ocean baroclinic solutions
*       Confirms that 2-D solver is giving correct 1-D solution
        if ( psoloc ) then
          print *,' '
          write (*,*) ' pch1oc:'
          do j=min(241,nypo),1,-1
            write (*,260) j,(ocwk1(i,j),i=1,3),
     &                    (ocwk1(i,j),i=nxpo-1,nxpo)
          enddo
          print *,' '
          write (*,*) ' pch2oc:'
          do j=min(241,nypo),1,-1
            write (*,260) j,(ocwk2(i,j),i=1,3),
     &                    (ocwk2(i,j),i=nxpo-1,nxpo)
          enddo
        endif
*       Compute area integrals of pch1oc and pch2oc
        call xintp (aipch1, ocwk1, nxpo, nypo)
        call xintp (aipch2, ocwk2, nxpo, nypo)
*       Both solutions should have the same area integral
        aipcho(m) = 0.5d0*(aipch1+aipch2)*dxo*dyo

*       Compute dp/dy half a gridpoint in from the north
*       and south boundaries, and "integrate" in x
*       Since these solutions are independent of x,
*       x integration means just multiply by xlo
        pch1ys = ( pch1oc(  2 ,m) - pch1oc(   1  ,m) )/dyo
        pch2ys = ( pch2oc(  2 ,m) - pch2oc(   1  ,m) )/dyo
        pch1yn = ( pch1oc(nypo,m) - pch1oc(nypo-1,m) )/dyo
        pch2yn = ( pch2oc(nypo,m) - pch2oc(nypo-1,m) )/dyo
*       Correction for baroclinic modes
        pch1ys = -pch1ys + 0.5d0*dyo*rdm2oc(m+1)*pch1oc(  1 ,m)
        pch2ys = -pch2ys + 0.5d0*dyo*rdm2oc(m+1)*pch2oc(  1 ,m)
        pch1yn =  pch1yn + 0.5d0*dyo*rdm2oc(m+1)*pch1oc(nypo,m)
        pch2yn =  pch2yn + 0.5d0*dyo*rdm2oc(m+1)*pch2oc(nypo,m)
*       Convert to line integrals
        pch1ys = xlo*pch1ys
        pch2ys = xlo*pch2ys
        pch1yn = xlo*pch1yn
        pch2yn = xlo*pch2yn
*       The above are (for each mode m) the quantities in square
*       brackets on the RHS of (B.14) and (B.15) (ocean version)
        pchdet = pch1ys*pch2yn - pch2ys*pch1yn
        hc1soc(m) = pch1ys/pchdet
        hc2soc(m) = pch2ys/pchdet
        hc1noc(m) = pch1yn/pchdet
        hc2noc(m) = pch2yn/pchdet
        print *,' '
        write (*,240) '  aipch1, aipch2 = ',aipch1,aipch2
        write (*,240) '  aipcho         = ',aipcho(m)
        write (*,240) '  pch1ys, pch1yn = ',pch1ys,pch1yn
        write (*,240) '  pch2ys, pch2yn = ',pch2ys,pch2yn
        write (*,240) '  pchdet         = ',pchdet
        write (*,240) '  hc1soc, hc2soc = ',hc1soc(m),hc2soc(m)
        write (*,240) '  hc1noc, hc2noc = ',hc1noc(m),hc2noc(m)
      enddo
#  else
*     Finite box ocean
*     ----------------
      print *,' '
      write (*,*) ' Ocean homogeneous (baroclinic) solutions:'
      do m=1,nlo-1
        print *,' '
        write (*,'(a,i2)') '  Mode: ',m
        do i=1,nxto
          boc(i) = bd2oc(i) - rdm2oc(m+1)
        enddo

*       Compute new homogeneous solution = (1 + rdm2*sol0)
*       sol0 satisfies Del-sqd(sol0) - rdm2*sol0 = 1
*       with the usual boundary condition p = 0
*       These are baroclinic solutions for use in ocinvq
*       Setup RHS.
!$OMP   PARALLEL DO DEFAULT (NONE)
!$OMP&           PRIVATE (i,j)
!$OMP&           SHARED  (ochom,m)
!$OMP&           SCHEDULE (STATIC)
        do j=1,nypo
          do i=1,nxpo
            ochom(i,j,m) = 1.0d0
          enddo
        enddo
!$OMP   END PARALLEL DO
*       Solve for sol0 in ochom.
        call hsbxoc (ochom(1,1,m), boc)
*       Add constant offset
!$OMP   PARALLEL DO DEFAULT (NONE)
!$OMP&           PRIVATE (i,j)
!$OMP&           SHARED  (ochom,m,rdm2oc)
!$OMP&           SCHEDULE (STATIC)
        do j=1,nypo
          do i=1,nxpo
            ochom(i,j,m) = 1.0d0 + rdm2oc(m+1)*ochom(i,j,m)
          enddo
        enddo
!$OMP   END PARALLEL DO
*       Check structure of homogeneous ocean baroclinic solution
        if ( psoloc ) then
          print *,' '
          write (*,*) ' ochom:'
          do j=min(241,nypo),1,-1
            write (*,250) j,(ochom(i,j,m),i=1,6)
          enddo
        endif
*       Area integral of full homogeneous solution
        call xintp (aipohs(m), ochom(1,1,m), nxpo, nypo)
        aipohs(m) = aipohs(m)*dxo*dyo
        write (*,240) '  aipohs = ',aipohs(m)
      enddo
*     Compute the matrices used in the mass constraint equation
*     dpioc(k) = Area integral of pressure diff ( po(k+1) - po(k) )
*     Choose sign of dpioc so that +ve dpioc(k) -> +ve eta(k)
*     cdiffo(m,k) is  coefficient which multiplies the mode
*     m amplitude to give its contribution to dpioc(k)
*     cdhoc(k,m) is  coefficient which multiplies a homogeneous
*     baroclinic mode coefficient to give its contribution to dpioc(k)
      do k=1,nlo-1
        do m=1,nlo
          cdiffo(m,k) = ctm2loc(m,k+1) - ctm2loc(m,k)
        enddo
        do m=1,nlo-1
          cdhoc(k,m) = ( ctm2loc(m+1,k+1) - ctm2loc(m+1,k) )*aipohs(m)
          cdhlu(k,m) = cdhoc(k,m)
        enddo
      enddo
      print *,' '
      write (*,*) ' Mass constraint matrices:'
      print *,' '
      write (*,*) ' cdiffo:'
      do k=1,nlo-1
        write (*,'(2x,i2,1x,1p,9d17.9)') k,(cdiffo(m,k),m=1,nlo)
      enddo
      print *,' '
      write (*,*) ' cdhoc:'
      do k=1,nlo-1
        write (*,'(2x,i2,1x,1p,9d17.9)') k,(cdhoc(k,m),m=1,nlo-1)
      enddo
*     Compute the LU factorization of cdhoc
*     DGETRF = NAG routine F07ADF
      call DGETRF (nlo-1, nlo-1, cdhlu, nlo-1, ipivch, info)
      if ( info.ne.0 ) then
        print *,'  DGETRF for ocean in homsol returns info = ',info
        print *,'  program terminates in homsol'
        stop
      endif
      print *,' '
      write (*,*) ' cdhlu:'
      do k=1,nlo-1
        write (*,'(2x,i2,1x,1p,9d17.9)') k,(cdhlu(k,m),m=1,nlo-1)
      enddo
      print *,' '
      write (*,*) ' ipivch:'
      write (*,'(4x,9i4)') (ipivch(k),k=1,nlo-1)
#  endif
#endif

#ifndef ocean_only
*     Compute homogeneous atmospheric channel solutions
*     =================================================
*     Barotropic mode: homogeneous solutions would need
*     to satisfy d2p/dy2 = 0, but the coefficient of
*     integral of p vanishes, so only the homogeneous
*     component with a derivative can be inferred.
*     W.L.O.G. we choose this to be 1 on the Southern
*     boundary, and to vanish on the Northern boundary.
*     We then apply the constraint at the Southern boundary
      do j=1,nypa
        pbhat(j) = dble(nypa-j)/dble(nypa-1)
      enddo
*     hbs = integral of y-derivative along Southern bdy
      hbsiat = yla/xla
      aipbha = 0.5d0*xla*yla
      print *,' '
      write (*,*) ' Atmospheric homogeneous barotropic solution:'
      write (*,240) '  aipbha         = ',aipbha
      write (*,240) '  hbsiat         = ',hbsiat
*     Check structure of homogeneous atmos. barotropic solution
      if ( psolat ) then
        print *,' '
        write (*,*) ' pbhat:'
        do j=min(121,nypa),1,-1
          write (*,250) j,pbhat(j)
        enddo
      endif

*     Now compute baroclinic solutions for use in atinvq
*     Solutions are functions of y only, but use the usual
*     2-D Helmholtz solver for simplicity and consistency

      print *,' '
      write (*,*) ' Atmospheric homogeneous (baroclinic) solutions:'
      do m=1,nla-1
        print *,' '
        write (*,'(a,i2)') '  Mode: ',m
        do i=1,nxta
          bat(i) = bd2at(i) - rdm2at(m+1)
        enddo

*       pch = L(y) + rdm2*sol0, where L(y) is linear in y,
*       and sol0 satisfies Del-sqd(sol0) - rdm2*sol0 = L(y)
*       with the usual boundary condition p = 0.
*       For pch1at, L(y) = 1.0 on S bdy; = 0.0 on N bdy
*       For pch2at, L(y) = 0.0 on S bdy; = 1.0 on N bdy
*       Specify RHSs
!$OMP   PARALLEL DO DEFAULT (NONE)
!$OMP&           PRIVATE (i,j)
!$OMP&           SHARED  (m,pch1at,pch2at,yla,ypa,atwk1,atwk2)
!$OMP&           SCHEDULE (STATIC)
        do j=1,nypa
          pch1at(j,m) = ( ypa(nypa)-ypa(j) )/yla
          pch2at(j,m) = ( ypa(  j )-ypa(1) )/yla
          do i=1,nxpa
            atwk1(i,j) = pch1at(j,m)
            atwk2(i,j) = pch2at(j,m)
          enddo
        enddo
!$OMP   END PARALLEL DO
*       Check r.h.s. for homogeneous atmos. baroclinic solutions
        if ( prhsat ) then
          print *,' '
          write (*,*) ' r.h.s. for pch1at:'
          do j=min(121,nypa),1,-1
            write (*,260) j,(atwk1(i,j),i=1,3),
     &                    (atwk1(i,j),i=nxpa-1,nxpa)
          enddo
          print *,' '
          write (*,*) ' r.h.s. for pch2at:'
          do j=min(121,nypa),1,-1
            write (*,260) j,(atwk2(i,j),i=1,3),
     &                    (atwk2(i,j),i=nxpa-1,nxpa)
          enddo
        endif

*       Invert these RHSs for baroclinic homog. solutions (sol0 above)
        call hscyat (atwk1, bat)
        call hscyat (atwk2, bat)
*       Check Helmholtz component of homogeneous atmos. baroclinic solutions
        if ( phelat ) then
          print *,' '
          write (*,*) ' Helmholtz for pch1at:'
          do j=min(121,nypa),1,-1
            write (*,260) j,(atwk1(i,j),i=1,3),
     &                    (atwk1(i,j),i=nxpa-1,nxpa)
          enddo
          print *,' '
          write (*,*) ' Helmholtz for pch2at:'
          do j=min(121,nypa),1,-1
            write (*,260) j,(atwk2(i,j),i=1,3),
     &                    (atwk2(i,j),i=nxpa-1,nxpa)
          enddo
        endif

*       Add Helmholtz solution to L(y) to get full solutions
*       Solutions in atwk1, atwk2 are functions of y only, i.e.
*       independent of i, so just save solution for one i value
!$OMP   PARALLEL DO DEFAULT (NONE)
!$OMP&           PRIVATE (i,j)
!$OMP&           SHARED  (m,pch1at,pch2at,atwk1,atwk2,rdm2at)
!$OMP&           SCHEDULE (STATIC)
        do j=1,nypa
          do i=1,nxpa
            atwk1(i,j) = pch1at(j,m) + rdm2at(m+1)*atwk1(i,j)
            atwk2(i,j) = pch2at(j,m) + rdm2at(m+1)*atwk2(i,j)
          enddo
          pch1at(j,m) = atwk1(1,j)
          pch2at(j,m) = atwk2(1,j)
        enddo
!$OMP   END PARALLEL DO
*       Check structure of homogeneous atmos. baroclinic solutions
*       Confirms that 2-D solver is giving correct 1-D solution
        if ( psolat ) then
          print *,' '
          write (*,*) ' pch1at:'
          do j=min(121,nypa),1,-1
            write (*,260) j,(atwk1(i,j),i=1,3),
     &                    (atwk1(i,j),i=nxpa-1,nxpa)
          enddo
          print *,' '
          write (*,*) ' pch2at:'
          do j=min(121,nypa),1,-1
            write (*,260) j,(atwk2(i,j),i=1,3),
     &                    (atwk2(i,j),i=nxpa-1,nxpa)
          enddo
        endif
*       Compute area integrals of pch1at and pch2at
        call xintp (aipch1, atwk1, nxpa, nypa)
        call xintp (aipch2, atwk2, nxpa, nypa)
*       Both solutions should have the same area integral
        aipcha(m) = 0.5d0*(aipch1+aipch2)*dxa*dya

*       Compute dp/dy half a gridpoint in from the north
*       and south boundaries, and "integrate" in x
*       Since these solutions are independent of x,
*       x integration means just multiply by xla
        pch1ys = ( pch1at(  2 ,m) - pch1at(   1  ,m) )/dya
        pch2ys = ( pch2at(  2 ,m) - pch2at(   1  ,m) )/dya
        pch1yn = ( pch1at(nypa,m) - pch1at(nypa-1,m) )/dya
        pch2yn = ( pch2at(nypa,m) - pch2at(nypa-1,m) )/dya
*       Correction for baroclinic modes
        pch1ys = -pch1ys + 0.5d0*dya*rdm2at(m+1)*pch1at(  1 ,m)
        pch2ys = -pch2ys + 0.5d0*dya*rdm2at(m+1)*pch2at(  1 ,m)
        pch1yn =  pch1yn + 0.5d0*dya*rdm2at(m+1)*pch1at(nypa,m)
        pch2yn =  pch2yn + 0.5d0*dya*rdm2at(m+1)*pch2at(nypa,m)
*       Convert to line integrals
        pch1ys = xla*pch1ys
        pch2ys = xla*pch2ys
        pch1yn = xla*pch1yn
        pch2yn = xla*pch2yn
*       The above are (for each mode m) the quantities in square
*       brackets on the RHS of (B.14) and (B.15) (atmos. version)
        pchdet = pch1ys*pch2yn - pch2ys*pch1yn
        hc1sat(m) = pch1ys/pchdet
        hc2sat(m) = pch2ys/pchdet
        hc1nat(m) = pch1yn/pchdet
        hc2nat(m) = pch2yn/pchdet
        print *,' '
        write (*,240) '  aipch1, aipch2 = ',aipch1,aipch2
        write (*,240) '  aipcha         = ',aipcha(m)
        write (*,240) '  pch1ys, pch1yn = ',pch1ys,pch1yn
        write (*,240) '  pch2ys, pch2yn = ',pch2ys,pch2yn
        write (*,240) '  pchdet         = ',pchdet
        write (*,240) '  hc1sat, hc2sat = ',hc1sat(m),hc2sat(m)
        write (*,240) '  hc1nat, hc2nat = ',hc1nat(m),hc2nat(m)
      enddo
#endif

  240 format(a,1p,9d21.13)
  250 format(i6,1x,6f12.7)
  260 format(i6,1x,1p,5d14.6)

      END SUBROUTINE homsol
c
c***********************************************************************
c
      END MODULE conhoms
c
c***********************************************************************
