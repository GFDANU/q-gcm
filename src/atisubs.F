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
      MODULE atisubs

*     Contains subprograms to invert atmospheric potential vorticity
*     qa to give atmospheric pressure pa, and (as part of that process)
*     to solve the modified Helmholtz equation for a given RHS,
*     using the Fourier transform and tridiagonal solver
*     method outlined in section 8 of the UserGuide.
*     This version uses the FFTPACK Fourier transform package.
*     The choice of package determines
*     the ordering of the values in bd2at

*     Modules
      USE parameters, ONLY : nxta

      IMPLICIT NONE

      PRIVATE

#ifndef ocean_only
      PUBLIC :: atinvq, hscyat
#endif

*     Storage for atmospheric FFT coefficients + workspace
      integer, PARAMETER, PUBLIC :: lwftat = 2*nxta + 15
      double precision, PUBLIC :: aftwrk(lwftat)
*
*     Coefficients for the tridiagonal equation
      double precision, PUBLIC :: aat,bd2at(nxta)

#ifndef ocean_only
      CONTAINS
#endif

c***********************************************************************
c
#ifndef ocean_only

      SUBROUTINE atinvq
*
*     Invert atmospheric potential vorticity qa to give atmospheric
*     pressure pa. i.e. calculates LHS of equation (7.15) and then
*     solves the modified Helmholtz equation for modal pressures.
*     Suitable multiples of the homogeneous solutions are added
*     to satisfy the momentum constraints at the zonal boundaries.
*     It then transforms the pressures from the modal to the
*     layer representation.
*     hscyat does the actual work of solving the Helmholtz
*     equation in a periodic channel domain, to give pa.

*     Modules
      USE parameters, ONLY : fnot, beta
      USE atconst
      USE athomog
      USE atstate, ONLY : pa, pam, qa
      USE intsubs, ONLY : xintp
      USE monitor, ONLY : ermasa, emfrat

      IMPLICIT NONE

*     Subroutine arguments
*
*     Local parameters
      double precision ecrita
      parameter ( ecrita=1.0d-13 )

*     Local variables
      integer i,j,k,m
      double precision betay,ql(nla),qm,wrk(nxpa,nypa,nla),bat(nxta),
     &                 pm(nla),pl,ainhom(nla),ayis(nla),ayin(nla),
     &                 entfac,rhss(nla),rhsn(nla),atsnew(nla),
     &                 atnnew(nla),clhss(nla),clhsn(nla),c1(nla-1),
     &                 c2(nla-1),c3,aipmod(nla),aiplay(nla),est1,est2,
     &                 edif,esum,homcor(nla)

c     PRESSURES FROM PV (ATMOSPHERE)

*     Compute vorticity RHS for each mode - equation (8.13)
*     =====================================================
*     Values not needed on zonal boundaries
!$OMP PARALLEL DO DEFAULT (NONE)
!$OMP&         PRIVATE (i,j,k,m,betay,ql,qm)
!$OMP&         SHARED  (yparel,qa,ctl2mat,wrk,ddynat)
!$OMP&         SCHEDULE (STATIC)
      do j=2,nypa-1
        betay = beta*yparel(j)
        do i=1,nxpa
*         Get layer contributions
          do k=1,nla
            ql(k) = qa(i,j,k) - betay
          enddo
          ql(1) = ql(1) - ddynat(i,j)
*         Compute modal values
          do m=1,nla
            qm = 0.0d0
            do k=1,nla
              qm = qm + ctl2mat(k,m)*ql(k)
            enddo
            wrk(i,j,m) = fnot*qm
          enddo
        enddo
      enddo
!$OMP END PARALLEL DO

*     Solve modified Helmholtz equation to invert each mode
*     =====================================================
*     m = 1 is the barotropic mode
*     m=2,nla are baroclinic modes

      do m=1,nla
*       Compute Helmholtz operator for mode m
        do i=1,nxta
          bat(i) = bd2at(i) - rdm2at(m)
        enddo
*       Overwrite wrk with new modal pressure
        call hscyat (wrk(1,1,m), bat)
*       Compute area integral of new inhomogeneous solution
        call xintp (ainhom(m), wrk(1,1,m), nxpa, nypa)
        ainhom(m) = ainhom(m)*dxa*dya
*       Compute line integrals of p_y for new modal solutions
*       Integrate along south & north boundaries for all modes
*       2-point formulation, but values on bdy are exactly zero
        ayis(m) =  0.5d0*wrk(1,   2  ,m)
        ayin(m) = -0.5d0*wrk(1,nypa-1,m)
        do i=2,nxpa-1
          ayis(m) = ayis(m) + wrk(i,   2  ,m)
          ayin(m) = ayin(m) - wrk(i,nypa-1,m)
        enddo
        ayis(m) = ayis(m) + 0.5d0*wrk(nxpa,   2  ,m)
        ayin(m) = ayin(m) - 0.5d0*wrk(nxpa,nypa-1,m)
        ayis(m) = ayis(m)*(dxa/dya)
        ayin(m) = ayin(m)*(dxa/dya)
      enddo
*     Have solved inhomogeneous modal problem with pa = 0 zonal bc

*     Solve constraint equations and add homogeneous solutions
*     ========================================================

*     Compute homogeneous solution coefficients
*     -----------------------------------------
*     Accumulate RHSs for the constraint equations for each layer
      entfac = 0.5d0*dya*fnot*fnot
*     Bottom layer
      rhss(1) = - (entfac/hat(1))*enisat(1) - (fnot/hat(1))*txisat
     &          + ajisat(1) + ap5sat(1)
      rhsn(1) = - (entfac/hat(1))*eninat(1) + (fnot/hat(1))*txinat
     &          + ajinat(1) - ap5nat(1)
*     Intermediate layers
      do k=2,nla-1
        rhss(k) = - (entfac/hat(k))*( enisat(k) - enisat(k-1) )
     &            + ajisat(k) + ap5sat(k)
        rhsn(k) = - (entfac/hat(k))*( eninat(k) - eninat(k-1) )
     &            + ajinat(k) - ap5nat(k)
      enddo
*     Top layer
      rhss(nla) =   (entfac/hat(nla))*enisat(nla-1)
     &            + ajisat(nla) + ap5sat(nla)
      rhsn(nla) =   (entfac/hat(nla))*eninat(nla-1)
     &            + ajinat(nla) - ap5nat(nla)
**    write (*,299) '  a rhss = ',rhss
**    write (*,299) '  a rhsn = ',rhsn

*     Update the constraint vectors
      do k=1,nla
        atsnew(k) = atmcsp(k) + tdta*rhss(k)
        atnnew(k) = atmcnp(k) + tdta*rhsn(k)
        atmcsp(k) = atmcs(k)
        atmcnp(k) = atmcn(k)
        atmcs(k) = atsnew(k)
        atmcn(k) = atnnew(k)
      enddo
**    write (*,299) '  atsnew = ',atsnew
**    write (*,299) '  atnnew = ',atnnew
*
*     Compute LHSs for the c1, c2 equations
      do m=1,nla
        clhss(m) = 0.0d0
        clhsn(m) = 0.0d0
        do k=1,nla
          clhss(m) = clhss(m) + ctl2mat(k,m)*atsnew(k)
          clhsn(m) = clhsn(m) + ctl2mat(k,m)*atnnew(k)
        enddo
        clhss(m) = clhss(m) + ayis(m)
        clhsn(m) = clhsn(m) - ayin(m)
      enddo
*
*     Get coefft for barotropic mode
      c3 = clhss(1)*hbsiat
*     Derive c1, c2 for baroclinic modes
      do m=1,nla-1
        c1(m) = hc2nat(m)*clhss(m+1) - hc2sat(m)*clhsn(m+1)
        c2(m) = hc1sat(m)*clhsn(m+1) - hc1nat(m)*clhss(m+1)
      enddo
*
*     Compute area integrals of pressures
*     -----------------------------------
*     Integrals of modal pressures
      aipmod(1) = ainhom(1) + c3*aipbha
      do m=2,nla
        aipmod(m) = ainhom(m) + ( c1(m-1) + c2(m-1) )*aipcha(m-1)
      enddo
*     Integrals of layer pressures
      do k=1,nla
        pl = 0.0d0
        do m=1,nla
          pl = pl + ctm2lat(m,k)*aipmod(m)
        enddo
        aiplay(k) = pl
      enddo
*
*     Check continuity is satisfied at each interface
*     Update continuity constants
      do k=1,nla-1
*       MONITORING - extra section for ermasa, emfrat
*       Compute alternative estimates of new dpiat
        est1 = aiplay(k) - aiplay(k+1)
        est2 = dpiatp(k) - tdta*gpat(k)*xan(k)
        edif = est1 - est2
        esum = abs(est1) + abs(est2)
        ermasa(k) = edif
*       Compute fractional error if entrainment is significant;
*       fraction is meaningless if est1, est2 just noisy zeros
        if ( esum.gt.(ecrita*xla*yla*tdta*gpat(k)) ) then
          emfrat(k) = 2.0d0*edif/esum
         else
          emfrat(k) = 0.0d0
        endif
**      write (*,'(a,i4,1p,3d16.6)') '  k, dpdt atm forcing = ',
**   &                               k,tdta*gpat(k)*xan(k)
**      write (*,'(a,i4,1p,3d16.6)') '  k, est1, est2, edif = ',
**   &                               k,est1,est2,edif
**      write (*,'(a,i4,1p,3d16.6)') '  k,  ermasa,  emfrat = ',
**   &                               k,ermasa(k),emfrat(k)
*       Update continuity constants
        dpiatp(k) = dpiat(k)
        dpiat(k) = aiplay(k) - aiplay(k+1)
      enddo

*     Add suitable multiple of homogeneous solutions and unpack modal
*     pressures to layer pressures - equations (7.18) and (7.19)
*     Also copy current pa to pam before overwriting with new value
*     -------------------------------------------------------------
!$OMP PARALLEL DO DEFAULT (NONE)
!$OMP&         PRIVATE (i,j,k,m,homcor,pm,pl)
!$OMP&         SHARED  (wrk,c1,pch1at,c2,pch2at,c3,pbhat,ctm2lat,pam,pa)
!$OMP&         SCHEDULE (STATIC)
      do j=1,nypa
*       Compute homogeneous corrections (indep. of i)
*       Barotropic mode
        homcor(1) = c3*pbhat(j)
*       Baroclinic modes
        do m=2,nla
          homcor(m) = c1(m-1)*pch1at(j,m-1) + c2(m-1)*pch2at(j,m-1)
        enddo
        do i=1,nxpa
          do m=1,nla
            pm(m) = wrk(i,j,m) + homcor(m)
          enddo
*         Compute layer values
          do k=1,nla
            pam(i,j,k) = pa(i,j,k)
            pl = 0.0d0
            do m=1,nla
              pl = pl + ctm2lat(m,k)*pm(m)
            enddo
            pa(i,j,k) = pl
          enddo
        enddo
      enddo
!$OMP END PARALLEL DO

      END SUBROUTINE atinvq

#endif /* not ocean_only */
c
c***********************************************************************
c
#ifndef ocean_only

      SUBROUTINE hscyat (wrk, bat)
*
*     Solves the inhomogeneous Helmholtz equation
*     for given rhs in a zonally periodic domain.
*     On entry wrk contains the rhs for a given mode.
*     On exit wrk contains the modal pressure
*     solution, including zonal boundary values.
*     aat and bat are the coefficients of the
*     Fourier transformed form of the equation.
*     Version specifically optimised for atmosphere arrays
*
*     This version uses the FFTPACK routines DRFFTF/DRFFTB

*     Modules
      USE parameters, ONLY : nxta, nxpa, nypa
*
      IMPLICIT NONE

*     Subroutine arguments
      double precision, INTENT(INOUT) :: wrk(nxpa,nypa)
      double precision, INTENT(IN) :: bat(nxta)
*
*     Local parameters
      double precision ftnorm
      parameter ( ftnorm = 1.0d0/nxta )
*
*     Local variables
      integer i,j
      double precision afttmp(lwftat)
*     Variables for tridiagonal solver
      double precision gam(nypa-1),betinv,uvec(nypa-1)
*
      do i=1,lwftat
        afttmp(i) = aftwrk(i)
      enddo
!$OMP PARALLEL DEFAULT (NONE)
!$OMP&         PRIVATE (i,j,uvec) FIRSTPRIVATE (afttmp)
!$OMP&         SHARED  (wrk,bat,aat)

*     Compute FFT of rhs along latitude lines
*     ---------------------------------------
!$OMP DO SCHEDULE (STATIC)
      do j=2,nypa-1
        call drfftf (nxta, wrk(1,j), afttmp)
      enddo
!$OMP END DO
*
*     For each wavenumber i, solve for FFT of p
*     -----------------------------------------
!$OMP DO SCHEDULE (STATIC)
!$OMP&   PRIVATE (betinv,gam)
      do i=1,nxta
*       Compute solution in vector uvec
*       Tridiagonal solver routine now inlined.
        betinv = 1.0d0/bat(i)
        uvec(2) = wrk(i,2)*betinv
*       Decomposition and forward substitution.
        do j=3,nypa-1
          gam(j) = aat*betinv
          betinv = 1.0d0/( bat(i) - aat*gam(j) )
          uvec(j) = ( wrk(i,j) - aat*uvec(j-1) )*betinv
        enddo
*       Backsubstitution.
        do j=nypa-2,2,-1
          uvec(j) = uvec(j) - gam(j+1)*uvec(j+1)
        enddo
*       Copy back solution and rescale
        do j=2,nypa-1
          wrk(i,j) = ftnorm*uvec(j)
        enddo
      enddo
!$OMP END DO
*
*     Inverse FFT solution along latitude lines
*     -----------------------------------------
*     Non-ideal cache access following previous loop
*     is spread over extensive floating-point work
!$OMP DO SCHEDULE (STATIC)
      do j=2,nypa-1
        call drfftb (nxta, wrk(1,j), afttmp)
*       Impose cyclic condition E-W
        wrk(nxpa,j) = wrk(1,j)
      enddo
!$OMP END DO NOWAIT

!$OMP END PARALLEL
*
*     Impose N & S boundary values, which are
*     implicit in the tridiagonal formulation.
      do i=1,nxpa
        wrk(i,  1 ) = 0.0d0
        wrk(i,nypa) = 0.0d0
      enddo

      END SUBROUTINE hscyat

#endif /* not ocean_only */
c
c***********************************************************************
c
      END MODULE atisubs
c
c***********************************************************************
