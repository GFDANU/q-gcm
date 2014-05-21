c***********************************************************************
c     Q-GCM Version 1.5.0 : last modified 07/03/2013
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
      MODULE ocisubs

*     Contains subprograms to invert oceanic potential vorticity
*     qo to give oceanic pressure po, and (as part of that process)
*     to solve the modified Helmholtz equation for a given RHS,
*     using the Fourier transform and tridiagonal solver
*     method outlined in section 8 of the UserGuide.
*     This version uses the FFTPACK Fourier transform package.
*     The choice of FFT package determines
*     the ordering of the values in bd2oc

*     Modules
      USE parameters, ONLY : nxto

      IMPLICIT NONE

      PRIVATE

#ifndef atmos_only
      PUBLIC :: ocinvq
#  ifdef cyclic_ocean
      PUBLIC :: hscyoc
#  else
      PUBLIC :: hsbxoc
#  endif
#endif /* not atmos_only */

*     Storage for oceanic FFT coefficients + workspace
      integer, PARAMETER, PUBLIC :: lwftoc = 3*nxto + 15
      double precision, PUBLIC :: oftwrk(lwftoc)
*
*     Coefficients for the tridiagonal equation
      double precision, PUBLIC :: aoc,bd2oc(nxto)

#ifndef atmos_only
      CONTAINS
#endif

c***********************************************************************

#ifndef atmos_only
      SUBROUTINE ocinvq
*
*     Invert oceanic potential vorticity qo to give oceanic
*     pressure po. i.e. calculates LHS of equation (7.15) and then
*     solves the modified Helmholtz equation for modal pressures.
*     Suitable multiples of the homogeneous solutions
*     are added to satisfy appropriate constraints.
*     It then unpacks the pressures from the modal to the
*     layer representation.
*     hsbxoc does the actual work of solving the Helmholtz
*     eqn in a domain with a closed boundary, to give po
*     hscyoc solves for po in the cyclic domain case

*     Modules
      USE parameters
      USE occonst
      USE ochomog
      USE ocstate
      USE intsubs, ONLY : xintp
#  ifdef cyclic_ocean
      USE monitor, ONLY : ermaso, emfroc
#  endif

      IMPLICIT NONE

*     Subroutine arguments

*     Local parameters
      double precision ecrito
      parameter ( ecrito=1.0d-13 )

*     Local variables
      integer i,j,k,m
      double precision betay,ql(nlo),qm,wrk(nxpo,nypo,nlo),boc(nxto),
     &                 xinhom(nlo),pm(nlo),pl
#  ifdef cyclic_ocean
*     Zonally cyclic ocean
      double precision entfac,rhss(nlo),rhsn(nlo),ocsnew(nlo),
     &                 ocnnew(nlo),ayis,ayin,clhss(nlo),clhsn(nlo),
     &                 c1(nlo-1),c2(nlo-1),c3,aipmod(nlo),aiplay(nlo),
     &                 est1,est2,edif,esum,homcor(nlo)
#  else
*     Finite box ocean
      integer info,iwork(nlo)
      double precision aient(nlo-1),aitmp,rhsum,rhs(nlo-1),
     &                 hclco(nlo-1),ferr,berr,work(3*nlo)
#  endif

c     PRESSURES FROM PV (OCEAN)

*     Compute vorticity RHS for each mode - equation (8.13)
*     =====================================================
*     Values not needed on zonal boundaries
!$OMP PARALLEL DO DEFAULT (NONE)
!$OMP&         PRIVATE (i,j,k,m,betay,ql,qm)
!$OMP&         SHARED  (yporel,qo,ctl2moc,wrk,ddynoc)
!$OMP&         SCHEDULE (STATIC)
      do j=2,nypo-1
        betay = beta*yporel(j)
        do i=1,nxpo
*         Get layer contributions
          do k=1,nlo
            ql(k) = qo(i,j,k) - betay
          enddo
          ql(nlo) = ql(nlo) - ddynoc(i,j)
*         Compute modal values
          do m=1,nlo
            qm = 0.0d0
            do k=1,nlo
              qm = qm + ctl2moc(k,m)*ql(k)
            enddo
            wrk(i,j,m) = fnot*qm
          enddo
        enddo
      enddo
!$OMP END PARALLEL DO

*     Solve modified Helmholtz equation to invert each mode
*     =====================================================
*     m = 1 is the barotropic mode
*     m=2,nlo are baroclinic modes

      do m=1,nlo
*       Compute Helmholtz operator for mode m
        do i=1,nxto
          boc(i) = bd2oc(i) - rdm2oc(m)
        enddo
*       Overwrite wrk with new modal pressure
#  ifdef cyclic_ocean
*       Zonally cyclic ocean
        call hscyoc (wrk(1,1,m), boc)
#  else
*       Finite box ocean
        call hsbxoc (wrk(1,1,m), boc)
#  endif
*       Compute area integral of new inhomogeneous solution
        call xintp (xinhom(m), wrk(1,1,m), nxpo, nypo)
        xinhom(m) = xinhom(m)*dxo*dyo
      enddo
**    write (*,299) '  xinhom = ',xinhom
*     Have solved inhomogeneous modal problem
*     with po = 0 on all solid boundaries

*     Solve constraint equations and add homogeneous solutions
*     ========================================================
#  ifdef cyclic_ocean
*     Zonally cyclic ocean
*     Compute homogeneous solution coefficients
*     -----------------------------------------
*     Accumulate RHSs for the constraint equations for each layer
      entfac = 0.5d0*dyo*fnot*fnot
*     Top layer
      rhss(1) =   (entfac/hoc(1))*enisoc(1) + (fnot/hoc(1))*txisoc
     &          + ajisoc(1) - ap3soc(1) + ap5soc(1)
      rhsn(1) =   (entfac/hoc(1))*eninoc(1) - (fnot/hoc(1))*txinoc
     &          + ajinoc(1) + ap3noc(1) - ap5noc(1)
*     Intermediate layers
      do k=2,nlo-1
        rhss(k) =   (entfac/hoc(k))*( enisoc(k) - enisoc(k-1) )
     &            + ajisoc(k) - ap3soc(k) + ap5soc(k)
        rhsn(k) =   (entfac/hoc(k))*( eninoc(k) - eninoc(k-1) )
     &            + ajinoc(k) + ap3noc(k) - ap5noc(k)
      enddo
*     Bottom layer
      rhss(nlo) = - (entfac/hoc(nlo))*enisoc(nlo-1)
     &            + ajisoc(nlo) - ap3soc(nlo) + ap5soc(nlo)
     &            + (fnot/hoc(nlo))*bdrins
      rhsn(nlo) = - (entfac/hoc(nlo))*eninoc(nlo-1)
     &            + ajinoc(nlo) + ap3noc(nlo) - ap5noc(nlo)
     &            - (fnot/hoc(nlo))*bdrinn
**    write (*,299) '  o rhss = ',rhss
**    write (*,299) '  o rhsn = ',rhsn

*     Update the constraint vectors
*     -----------------------------
      do k=1,nlo
        ocsnew(k) = ocncsp(k) + tdto*rhss(k)
        ocnnew(k) = ocncnp(k) + tdto*rhsn(k)
        ocncsp(k) = ocncs(k)
        ocncnp(k) = ocncn(k)
        ocncs(k) = ocsnew(k)
        ocncn(k) = ocnnew(k)
      enddo
**    write (*,299) '  ocsnew = ',ocsnew
**    write (*,299) '  ocnnew = ',ocnnew
*
*     Compute LHSs for the c1, c2 equations
*     -------------------------------------
      do m=1,nlo
*       Compute line integrals of p_y for new modal solutions
*       Integrate along south & north boundaries for all modes
*       2-point formulation, but values on bdy are exactly zero
        ayis =  0.5d0*wrk(1,   2  ,m)
        ayin = -0.5d0*wrk(1,nypo-1,m)
        do i=2,nxpo-1
          ayis = ayis + wrk(i,   2  ,m)
          ayin = ayin - wrk(i,nypo-1,m)
        enddo
        ayis = ayis + 0.5d0*wrk(nxpo,   2  ,m)
        ayin = ayin - 0.5d0*wrk(nxpo,nypo-1,m)
        ayis = ayis*(dxo/dyo)
        ayin = ayin*(dxo/dyo)
        clhss(m) = 0.0d0
        clhsn(m) = 0.0d0
        do k=1,nlo
          clhss(m) = clhss(m) + ctl2moc(k,m)*ocsnew(k)
          clhsn(m) = clhsn(m) + ctl2moc(k,m)*ocnnew(k)
        enddo
        clhss(m) = clhss(m) + ayis
        clhsn(m) = clhsn(m) - ayin
      enddo
**    write (*,299) ' o clhss = ',clhss
**    write (*,299) ' o clhsn = ',clhsn
*
*     Get coefft for barotropic mode
      c3 = clhss(1)*hbsioc
*     Derive c1, c2 for baroclinic modes
      do m=1,nlo-1
        c1(m) = hc2noc(m)*clhss(m+1) - hc2soc(m)*clhsn(m+1)
        c2(m) = hc1soc(m)*clhsn(m+1) - hc1noc(m)*clhss(m+1)
      enddo
**    write (*,299) ' o c1(m) = ',c1
**    write (*,299) ' o c2(m) = ',c2
*
*     Compute area integrals of pressures
*     -----------------------------------
*     Integrals of modal pressures
      aipmod(1) = xinhom(1) + c3*aipbho
      do m=2,nlo
        aipmod(m) = xinhom(m) + ( c1(m-1) + c2(m-1) )*aipcho(m-1)
      enddo
**    write (*,299) '  aipmod = ',aipmod
*     Integrals of layer pressures
      do k=1,nlo
        pl = 0.0d0
        do m=1,nlo
          pl = pl + ctm2loc(m,k)*aipmod(m)
        enddo
        aiplay(k) = pl
      enddo
**    write (*,299) '  aiplay = ',aiplay
*
*     Update continuity measures at each interface
*     --------------------------------------------
      do k=1,nlo-1
*       Choose sign of dpioc so that +ve dpioc -> +ve eta
*       Check continuity is satisfied at each interface
*       MONITORING - extra section for ermaso, emfroc
*       Compute alternative estimates of new dpioc
        est1 = aiplay(k+1) - aiplay(k)
        est2 = dpiocp(k) - tdto*gpoc(k)*xon(k)
        edif = est1 - est2
        esum = abs(est1) + abs(est2)
        ermaso(k) = edif
*       Compute fractional error if entrainment is significant;
*       fraction is meaningless if est1, est2 just noisy zeros
        if ( esum.gt.(ecrito*xlo*ylo*tdto*gpoc(k)) ) then
          emfroc(k) = 2.0d0*edif/esum
         else
          emfroc(k) = 0.0d0
        endif
**      write (*,'(a,i4,1p,3d16.6)') '  k, dpdt ocn forcing = ',
**   &                               k,tdto*gpoc(k)*xon(k)
**      write (*,'(a,i4,1p,3d16.6)') '  k, est1, est2, edif = ',
**   &                               k,est1,est2,edif
**      write (*,'(a,i4,1p,3d16.6)') '  k,  ermaso,  emfroc = ',
**   &                               k,ermaso(k),emfroc(k)
*       Update continuity constants
        dpiocp(k) = dpioc(k)
        dpioc(k) = aiplay(k+1) - aiplay(k)
      enddo

*     Add suitable multiple of homogeneous solutions and unpack modal
*     pressures to layer pressures - equations (7.18) and (7.19)
*     Also copy current po to pom before overwriting with new value
*     -------------------------------------------------------------
!$OMP PARALLEL DO DEFAULT (NONE)
!$OMP&         PRIVATE (i,j,k,m,homcor,pm,pl)
!$OMP&         SHARED  (wrk,c1,pch1oc,c2,pch2oc,c3,pbhoc,ctm2loc,pom,po)
!$OMP&         SCHEDULE (STATIC)
      do j=1,nypo
*       Compute homogeneous corrections (indep. of i)
*       Barotropic mode
        homcor(1) = c3*pbhoc(j)
*       Baroclinic modes
        do m=2,nlo
          homcor(m) = c1(m-1)*pch1oc(j,m-1) + c2(m-1)*pch2oc(j,m-1)
        enddo
        do i=1,nxpo
          do m=1,nlo
            pm(m) = wrk(i,j,m) + homcor(m)
          enddo
*         Compute layer values
          do k=1,nlo
            pom(i,j,k) = po(i,j,k)
            pl = 0.0d0
            do m=1,nlo
              pl = pl + ctm2loc(m,k)*pm(m)
            enddo
            po(i,j,k) = pl
          enddo
        enddo
      enddo
!$OMP END PARALLEL DO
#  else
*     Finite box ocean
*     Get multiples of homogeneous solutions to conserve thickness
*     ------------------------------------------------------------
*     aient(k) = Area integral of oceanic entrainment e(k)
      aient(1) = xon(1)
*     N.B. xon(1) is now zero by construction of entoc in oml
*     All other entrainments assumed exactly zero in the ocean
***   aient(1) = 0.0d0
      do k=2,nlo-1
        aient(k) = 0.0d0
      enddo
**    write (*,299) '  aient  = ',aient
*     Area integral of d(eta(k))/dt = - Area integral of entrainment e(k)
      do k=1,nlo-1
        aitmp = dpioc(k)
        dpioc(k) = dpiocp(k) - tdto*gpoc(k)*aient(k)
        dpiocp(k) = aitmp
        rhsum = 0.0d0
        do m=1,nlo
          rhsum = rhsum + cdiffo(m,k)*xinhom(m)
        enddo
        rhs(k) = dpioc(k) - rhsum
        hclco(k) = rhs(k)
      enddo
**    write (*,299) '  dpioc  = ',dpioc
**    write (*,299) '  o rhs  = ',rhs
*     Matrix equation is cdhoc*hclco = rhs
*     Solve equation for homogeneous solution coeffts using LAPACK
*     Solve the linear system using the LU factorised matrix cdhlu
*     DGETRS = NAG routine F07AEF
      call DGETRS ('Norm', nlo-1, 1, cdhlu, nlo-1,
     &             ipivch, hclco, nlo-1, info)
      if ( info.ne.0 ) then
        print *,'  DGETRS in ocinvq returns info = ',info
        print *,'  program terminates in ocinvq'
        stop
      endif
*     Improve the solution by iterative refinement
*     DGERFS = NAG routine F07AHF
      call DGERFS ('Norm', nlo-1, 1, cdhoc, nlo-1, cdhlu, nlo-1,
     &             ipivch, rhs, nlo-1, hclco, nlo-1, ferr, berr,
     &             work, iwork, info)
**    write (*,299) '  hclco  = ',hclco

*     Add suitable multiple of homogeneous solution, and unpack modal
*     pressures to layer pressures - equations (7.18) and (7.19)
*     Also copy current po to pom before overwriting with new value
*     ---------------------------------------------------------------
!$OMP PARALLEL DO DEFAULT (NONE)
!$OMP&         PRIVATE (i,j,k,m,pm,pl)
!$OMP&         SHARED  (wrk,hclco,ochom,ctm2loc,pom,po)
!$OMP&         SCHEDULE (STATIC)
      do j=1,nypo
        do i=1,nxpo
*         Get modal contributions
*         Barotropic mode
          pm(1) = wrk(i,j,1)
*         Baroclinic modes including homogeneous contribution
          do m=2,nlo
            pm(m) = wrk(i,j,m) + hclco(m-1)*ochom(i,j,m-1)
          enddo
*         Compute layer values
          do k=1,nlo
            pom(i,j,k) = po(i,j,k)
            pl = 0.0d0
            do m=1,nlo
              pl = pl + ctm2loc(m,k)*pm(m)
            enddo
            po(i,j,k) = pl
          enddo
        enddo
      enddo
!$OMP END PARALLEL DO
#  endif

*     Format for debugging printouts
**299 format(a,1p,8d17.8)

      END SUBROUTINE ocinvq
#endif
c
c***********************************************************************
c
#ifndef atmos_only
#  ifndef cyclic_ocean

      SUBROUTINE hsbxoc (wrk, boc)
*
*     Solves the inhomogeneous Helmholtz equation for
*     given rhs in a domain with meridional boundaries
*     On entry wrk contains the rhs for a given mode.
*     On exit wrk contains the modal pressure
*     solution, including all boundary values.
*     boc are the coefficients of the sine
*     transformed version of the equation.
*     Only needed for box ocean case, hence the #ifndef
*
*     This version uses the FFTPACK routine DSINT

*     Modules
      USE parameters, ONLY : nxpo, nypo, nxto

      IMPLICIT NONE

*     Subroutine arguments
      double precision, INTENT(INOUT) :: wrk(nxpo,nypo)
      double precision, INTENT(IN) :: boc(nxto)
*
*     Local parameters
      character (len=*), parameter :: subnam = 'hsbxoc'
      double precision ftnorm
      parameter ( ftnorm = 0.5d0/nxto )
*
*     Local variables
      integer i,j
      double precision ofttmp(lwftoc)
*     Variables for tridiagonal solver
      double precision gam(nypo-1),betinv,uvec(nypo-1)
*
*     Create a local copy of the fft coeffts + workspace
      do i=1,lwftoc
        ofttmp(i) = oftwrk(i)
      enddo
!$OMP PARALLEL DEFAULT (NONE)
!$OMP&         PRIVATE (i,j,uvec) FIRSTPRIVATE (ofttmp)
!$OMP&         SHARED  (wrk,boc,aoc)

*     Compute sine transform of rhs along latitude lines
*     --------------------------------------------------
*     N.B. uses extra element wrk(nxpo,j) as workspace.
*     Value in this location does not affect inverse.
!$OMP DO SCHEDULE (STATIC)
      do j=2,nypo-1
        call dsint (nxto-1, wrk(2,j), ofttmp)
      enddo
!$OMP END DO
*
*     For each wavenumber i, solve for sine transform of p
*     ----------------------------------------------------
!$OMP DO SCHEDULE (STATIC)
!$OMP&   PRIVATE (betinv,gam)
      do i=2,nxpo-1
*       Compute solution in vector uvec
        betinv = 1.0d0/boc(i-1)
        uvec(2) = wrk(i,2)*betinv
*       Decomposition and forward substitution.
        do j=3,nypo-1
          gam(j) = aoc*betinv
          betinv = 1.0d0/( boc(i-1) - aoc*gam(j) )
          uvec(j) = ( wrk(i,j) - aoc*uvec(j-1) )*betinv
        enddo
*       Backsubstitution.
        do j=nypo-2,2,-1
          uvec(j) = uvec(j) - gam(j+1)*uvec(j+1)
        enddo
*       Copy back solution and rescale
        do j=2,nypo-1
          wrk(i,j) = ftnorm*uvec(j)
        enddo
      enddo
!$OMP END DO
*
*     Inverse sine transform along latitude lines
*     -------------------------------------------
!$OMP DO SCHEDULE (STATIC)
      do j=2,nypo-1
        call dsint (nxto-1, wrk(2,j), ofttmp)
*       Zero pressure on Western & Eastern boundaries
        wrk(  1 ,j) = 0.0d0
        wrk(nxpo,j) = 0.0d0
      enddo
!$OMP END DO NOWAIT

!$OMP END PARALLEL

*     Impose N & S boundary values, which are
*     implicit in the tridiagonal formulation.
      do i=1,nxpo
        wrk(i,  1 ) = 0.0d0
        wrk(i,nypo) = 0.0d0
      enddo
*

      END SUBROUTINE hsbxoc
#  endif /* not cyclic_ocean */
#endif /* not atmos_only */
c
c***********************************************************************
c
#ifndef atmos_only
#  ifdef cyclic_ocean

      SUBROUTINE hscyoc (wrk, boc)
*
*     Solves the inhomogeneous Helmholtz equation
*     for given rhs in a zonally periodic domain.
*     On entry wrk contains the rhs for a given mode.
*     On exit wrk contains the modal pressure
*     solution, including zonal boundary values.
*     boc contains the coefficients of the
*     Fourier transformed form of the equation.
*     Version specifically optimised for ocean arrays
*     Only needed for cyclic ocean case, hence the #ifdef
*
*     This version uses the FFTPACK routines DRFFTF/DRFFTB

*     Modules
      USE parameters, ONLY : nxpo, nypo, nxto

      IMPLICIT NONE

*     Subroutine arguments
      double precision, INTENT(INOUT) :: wrk(nxpo,nypo)
      double precision, INTENT(IN) :: boc(nxto)
*
*     Local parameters
      character (len=*), parameter :: subnam = 'hscyoc'
      double precision ftnorm
      parameter ( ftnorm = 1.0d0/nxto )
*
*     Local variables
      integer i,j
      double precision ofttmp(lwftoc)
*     Variables for tridiagonal solver
      double precision gam(nypo-1),betinv,uvec(nypo-1)
*
*     Create a local copy of the fft coeffts + workspace
      do i=1,lwftoc
        ofttmp(i) = oftwrk(i)
      enddo
!$OMP PARALLEL DEFAULT (NONE)
!$OMP&         PRIVATE (i,j,uvec) FIRSTPRIVATE (ofttmp)
!$OMP&         SHARED  (wrk,boc,aoc)

*     Compute FFT of rhs along latitude lines
*     ---------------------------------------
!$OMP DO SCHEDULE (STATIC)
      do j=2,nypo-1
        call drfftf (nxto, wrk(1,j), ofttmp)
      enddo
!$OMP END DO
*
*     For each wavenumber i, solve for FFT of p
*     -----------------------------------------
!$OMP DO SCHEDULE (STATIC)
!$OMP&   PRIVATE (betinv,gam)
      do i=1,nxto
*       Compute solution in vector uvec
        betinv = 1.0d0/boc(i)
        uvec(2) = wrk(i,2)*betinv
*       Decomposition and forward substitution.
        do j=3,nypo-1
          gam(j) = aoc*betinv
          betinv = 1.0d0/( boc(i) - aoc*gam(j) )
          uvec(j) = ( wrk(i,j) - aoc*uvec(j-1) )*betinv
        enddo
*       Backsubstitution.
        do j=nypo-2,2,-1
          uvec(j) = uvec(j) - gam(j+1)*uvec(j+1)
        enddo
*       Copy back solution and rescale
        do j=2,nypo-1
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
      do j=2,nypo-1
        call drfftb (nxto, wrk(1,j), ofttmp)
*       Impose cyclic condition E-W
        wrk(nxpo,j) = wrk(1,j)
      enddo
!$OMP END DO NOWAIT

!$OMP END PARALLEL
*
*     Impose N & S boundary values, which are
*     implicit in the tridiagonal formulation.
      do i=1,nxpo
        wrk(i,  1 ) = 0.0d0
        wrk(i,nypo) = 0.0d0
      enddo
*

      END SUBROUTINE hscyoc

#  endif /* cyclic_ocean */
#endif /* not atmos_only */
c
c***********************************************************************
c
      END MODULE ocisubs
c
c***********************************************************************
