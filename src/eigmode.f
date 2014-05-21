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
      MODULE eigmode

*     Contains subroutine eigmod for computing the vertical
*     eigenmodes and eigenvalues (for either the atmosphere or
*     the ocean), and the mode/layer transformation matrices

*     Modules

      IMPLICIT NONE

      PRIVATE

      PUBLIC :: eigmod

      CONTAINS

c***********************************************************************
*
      SUBROUTINE eigmod (nl, gpr, h, case, aaa, cphs,
     &                   rdef, rdm2, ctl2m, ctm2l)
*
*     Derive eigenmodes: eigenvalues and left & right
*     eigenvectors, which determine the coefficients
*     for transforming between layer and modal pressures
*     Works for both oceanic and atmospheric cases
*     This version includes the Flierl normalisation of
*     the right eigenvectors (left normalisation cancels)

*     Input arguments:
*     nl   : no. of QG layers
*     gpr  : vector of reduced gravities across QG layer interfaces (m s^-2)
*     h    : vector of unperturbed thicknesses of QG layers (m)
*     case : character string for identifying results as ocean or atmosphere
*     (all the above are unchanged on exit)

*     Output arguments:
*     aaa  : A(nl,nl) matrix linking pressures and eta
*     cph  : vector of nl modal phasespeed values (m s^-1)
*            (infinite value for barotropic mode replaced by 0.0)
*     rdef : vector of nl modal deformation radii (m)
*            (infinite value for barotropic mode replaced by 0.0)
*     rdm2 : vector of nl values of 1/(modal deformation radius)**2 (m^-2)
*     In each case 1st entry = barotropic mode; others are
*     baroclinic modes in order of decreasing modal speed
*     ctl2m : transpose of matrix for converting from layer to modal amplitudes
*     ctm2l : transpose of matrix for converting from modal to layer amplitudes

*     Modules
      USE parameters, ONLY : fnot

      IMPLICIT NONE

*     Subroutine arguments
      integer nl
      double precision gpr(nl-1),h(nl),aaa(nl,nl),cphs(nl),
     &                 rdef(nl),rdm2(nl),ctl2m(nl,nl),ctm2l(nl,nl)
      character (len=*) :: case

*     Local parameters
      integer nlmax,nblk,lwork
      parameter ( nlmax=9, nblk=64, lwork=nlmax*nblk )
      character (len=*), parameter :: baljob = 'Both'
      logical prtval,eigchk
      parameter ( prtval = .false., eigchk = .true. )
*
*     Local variables
      integer i,j,k,l,m,n,ilo,ihi,info,ncomp,mwant,mgot,index(nlmax),
     &        indtmp,inm
      double precision aba(nlmax,nlmax),scale(nlmax),
     &                 hup(nlmax,nlmax),tauvec(nlmax),work(lwork),
     &                 qqq(nlmax,nlmax),uqt(nlmax,nlmax),
     &                 zzz(nlmax,nlmax),wre(nlmax),wim(nlmax),
     &                 evecl(nlmax,nlmax),evecr(nlmax,nlmax),
     &                 eval(nlmax),c2rabs(nlmax),elder(nlmax,nlmax),
     &                 cl2m(nlmax,nlmax),cm2l(nlmax,nlmax),flfac,
     &                 ccprod(nlmax,nlmax),c2temp,dotp,htotal
      logical select(nlmax)
*     Extra variables for eigensolution verification
      double precision aevec(nlmax),eevec(nlmax)
*     Extra variables for triple product test
**    double precision ac(nlmax),ca

*     Check internal storage is sufficient
*     ------------------------------------
      if ( nl.gt.nlmax ) then
        print *,' eigmod has insufficient nlmax = ',nlmax
        print *,' called with nl = ',nl
        print *,' case = ',case
        print *,' program terminates in eigmod'
        stop
      endif

*     Setup A matrix
*     --------------
*     Zero all elements explicitly in A and associated matrices
      do j=1,nl
        do i=1,nl
          aaa(i,j) = 0.0d0
        enddo
      enddo
      do j=1,nlmax
        do i=1,nlmax
          aba(i,j) = 0.0d0
          hup(i,j) = 0.0d0
          qqq(i,j) = 0.0d0
          uqt(i,j) = 0.0d0
        enddo
      enddo
*     k = 1 layer (top(ocean) or bottom(atmos.))
      k = 1
      aaa(k,k+1) = -1.0d0/( gpr( k )*h(k) )
      aaa(k, k ) = - aaa(k,k+1)
*     Intermediate layers
      do k=2,nl-1
        aaa(k,k-1) = -1.0d0/( gpr(k-1)*h(k) )
        aaa(k,k+1) = -1.0d0/( gpr( k )*h(k) )
        aaa(k, k ) = - aaa(k,k-1) - aaa(k,k+1)
      enddo
*     k = nl layer (bottom(ocean) or top(atmos.))
      k = nl
      aaa(k,k-1) = -1.0d0/( gpr(k-1)*h(k) )
      aaa(k, k ) = - aaa(k,k-1)
*     Print A matrix as a check
      if ( prtval ) then
        write (*,'(a)') ' '
        write (*,'(2x,a,a)') 'case = ',case
        write (*,'(2x,a)') 'gprime vector:'
        write (*,'(2x,1p,9d16.7)') (gpr(k),k=1,nl-1)
        write (*,'(2x,a)') 'H vector:'
        write (*,'(2x,1p,9d16.7)') (h(k),k=1,nl)
        write (*,'(2x,a)') 'A matrix:'
        do i=1,nl
          write (*,'(2x,i3,1x,1p,9d16.7)') i,(aaa(i,j),j=1,nl)
        enddo
      endif

*     Solve for the eigenvalues and L & R eigenvectors of A
*     -----------------------------------------------------
*     (real nonsymmetric eigenproblem). Use LAPACK routines
*     DGEBAL = F08NHF
*     DGEHRD = F08NEF
*     DORGHR = F08NFF
*     DHSEQR = F08PEF
*     DTREVC = F08QKF
*     DGEBAK = F08NJF
*     Set default values of ILO, IHI.
      ilo = 1
      ihi = nl
*     ------------------------------------------------------
*     Balance a real general matrix to improve accuracy
*     Use a copy of A in case the original is required
*     BALJOB passed to DGEBAK to transform eigenvectors suitably
      do j=1,nl
        do i=1,nl
          aba(i,j) = aaa(i,j)
        enddo
      enddo
      call DGEBAL (baljob, nl, aba, nlmax, ilo, ihi, scale, info)
      if ( info.ne.0 ) then
        print *,' Problem in DGEBAL, INFO = ',info
        print *,' case = ',case
        print *,' program terminates in eigmod'
        stop
      endif
*     ABA is now overwritten with the balanced matrix
*     ------------------------------------------------------
*     Reduce a real general matrix to Hessenberg form
*     N.B. overwrites matrix with its upper Hessenberg form H
*     and details of the orthogonal transformation matrix Q
      do j=1,nl
        do i=1,nl
          hup(i,j) = aba(i,j)
        enddo
      enddo
      call DGEHRD (nl, ilo, ihi, hup, nlmax,
     &             tauvec, work, lwork, info)
      if ( info.ne.0 ) then
        print *,' Problem in DGEHRD, INFO = ',info
        print *,' case = ',case
        print *,' program terminates in eigmod'
        stop
***    else
***     print *,' For DGEHRD, min. optimal LWORK = ',work(1)
      endif
*     ------------------------------------------------------
*     Generate the real orthogonal matrix Q
*     that reduces A to upper Hessenberg form
*     Provide as input a copy of the upper Hessenberg matrix
*     previously generated, which also contains details of
*     the vectors which define the elementary reflectors
      do j=1,nl
        do i=1,nl
          qqq(i,j) = hup(i,j)
        enddo
      enddo
      call DORGHR (nl, ilo, ihi, qqq, nlmax,
     &             tauvec, work, lwork, info)
      if ( info.ne.0 ) then
        print *,' Problem in DORGHR, INFO = ',info
        print *,' case = ',case
        print *,' program terminates in eigmod'
        stop
***    else
***     print *,' For DORGHR, min. optimal LWORK = ',work(1)
      endif
*     ------------------------------------------------------
*     Compute all the eigenvalues, and optionally the Schur
*     factorization, of a real upper Hessenberg matrix
*     We require the Schur form in order to get the eigenvectors
*     This means that UQT is overwritten with the upper
*     quasi-triangular matrix T from the Schur decomposition
*     Set ZZZ to contain the matrix Q on entry.
      do j=1,nl
        do i=1,nl
          uqt(i,j) = hup(i,j)
          zzz(i,j) = qqq(i,j)
        enddo
      enddo
      call DHSEQR ('Schur', 'V', nl, ilo, ihi, uqt, nlmax,
     &             wre, wim, zzz, nlmax, work, lwork, info)
      if ( info.ne.0 ) then
        print *,' Problem in DHSEQR, INFO = ',info
        print *,' case = ',case
        print *,' program terminates in eigmod'
        stop
      endif
*     Check for reality of eigenvalues.
      ncomp = 0
      do m=1,nl
        if ( wim(m).ne.0.0d0 ) ncomp = ncomp + 1
      enddo
      if ( ncomp.ne.0 ) then
        print *,' DHSEQR returns complex eigenvalues'
        print *,' case = ',case
        print *,' program terminates in eigmod'
        print *,' Input parameters were:'
        write (*,'(2x,a)') 'gprime vector:'
        write (*,'(2x,1p,9d16.7)') (gpr(k),k=1,nl-1)
        write (*,'(2x,a)') 'H vector:'
        write (*,'(2x,1p,9d16.7)') (h(k),k=1,nl)
        stop
      endif
*     ------------------------------------------------------
*     Routine computes selected left and/or right eigen-
*     vectors of a real upper quasi-triangular matrix
*     We require ALL left AND right eigenvectors
*     EVECL and EVECR need to be set on entry to the
*     matrix of Schur vectors returned by DHSEQR,
*     which is not necessarily the same as Q
      do j=1,nl
        do i=1,nl
          evecl(i,j) = zzz(i,j)
          evecr(i,j) = zzz(i,j)
        enddo
      enddo
      mwant = nl
      call DTREVC ('Both', 'Back', select, nl, uqt, nlmax, evecl,
     &             nlmax, evecr, nlmax, mwant, mgot, work, info)
      if ( info.ne.0 ) then
        print *,' Problem in DTREVC, INFO = ',info
        print *,' case = ',case
        print *,' program terminates in eigmod'
        stop
      endif
*     Left and right eigenvectors are stored in the columns of
*     EVECL, EVECR, i.e. 1st subscript = eigenvector component,
*     2nd subscript = eigenmode number
*     ------------------------------------------------------
*     Transform the eigenvectors of a balanced matrix to
*     those of the original real nonsymmetric matrix A
      call DGEBAK (baljob, 'Left', nl, ilo, ihi,
     &             scale, mgot, evecl, nlmax, info)
      call DGEBAK (baljob, 'Righ', nl, ilo, ihi,
     &             scale, mgot, evecr, nlmax, info)
*     ------------------------------------------------------
*     At this point we should have a full set of eigen-
*     values and eigenvectors of the original matrix A
*
*     Apply the normalisation of Fleirl (1978) to the oceanic
*     right eigenvectors, and ensure they are +ve at the surface
*     ----------------------------------------------------------
*     Normalisation is equation (2.5) of "Models of Vertical Structure
*     and the Calibration of Two-layer Models" by Glenn R. Flierl,
*     Dynamics of Atmospheres and Oceans vol. 2 pp 341-381 (1978)
*     Surface sign convention follows the Rossby wave papers of
*     P.D. Killworth & J.R. Blundell in J. Physical Oceanography.
*     This is only immediately relevant for the ocean case.
      if ( case.eq.'Ocean' ) then
        htotal = 0.0d0
        do k=1,nl
          htotal = htotal + h(k)
        enddo
        do m=1,nl
*         Compute integrals
          dotp = 0.0d0
          do k=1,nl
            dotp = dotp + h(k)*evecr(k,m)*evecr(k,m)
          enddo
*         Just normalising the right eigenvectors; this should suffice.
*         Also ensure all right eigenvectors are +ve at the atmosphere/ocean
*         interface, to agree with the Killworth & Blundell sign convention.
          flfac = SIGN( sqrt(htotal/dotp), evecr(1,m) )
          do k=1,nl
            evecr(k,m) = flfac*evecr(k,m)
          enddo
        enddo
*       Check orthonormality of eigenvectors, using Flierl's convention
        if ( prtval ) then
          print *,' '
          write(*,'(a)') '  Check Flierl orthonormalisation:'
          do m=1,nl
            do n=m,nl
              dotp = 0.0d0
              do k=1,nl
                dotp = dotp + h(k)*evecr(k,m)*evecr(k,n)
              enddo
              dotp = dotp/htotal
              write(*,'(a,2i4,1p,d16.7)')
     &              '  modes, Sum( H(k)*Rm(k)*Rn(k) )/Depth = ',m,n,dotp
            enddo
          enddo
        endif
      endif
*
*     Check orthogonality of eigenvectors
*     -----------------------------------
      do n=1,nl
        do m=1,nl
          dotp = 0.0d0
          do k=1,nl
            dotp = dotp + evecl(k,m)*evecr(k,n)
          enddo
          elder(m,n) = dotp
        enddo
      enddo
*
*     Optionally print eigenvalues and vectors as a check
*     ---------------------------------------------------
      if ( prtval ) then
        write (*,'(a)') ' '
        write (*,'(2x,a,a)') 'case = ',case
        write (*,'(2x,a,i9,8i16)') ' m:  ',(m,m=1,nl)
        write (*,'(2x,a,1p,9d16.7)') ' Wre',(wre(m),m=1,nl)
        write (*,'(2x,a,1p,9d16.7)') ' Wim',(wim(m),m=1,nl)
        write (*,'(2x,a)') 'Right eigenvectors:'
        do k=1,nl
          write (*,'(2x,i3,1x,1p,9d16.7)') k,(evecr(k,m),m=1,nl)
        enddo
        write (*,'(2x,a)') 'Left  eigenvectors:'
        do k=1,nl
          write (*,'(2x,i3,1x,1p,9d16.7)') k,(evecl(k,m),m=1,nl)
        enddo
        write (*,'(2x,a)') 'Dot products evecl(m).evecr(n):'
        write (*,'(2x,a,i5,8i16)') '  m  n:  ',(n,n=1,nl)
        do m=1,nl
          write (*,'(2x,i3,1x,1p,9d16.7)') m,(elder(m,n),n=1,nl)
        enddo
      endif

*     Derive quantities required for Q-GCM
*     ------------------------------------
*     Eigenvalues are 1/c^2; we want these in increasing order
*     Barotropic mode is then first in list
      do m=1,nl
        c2rabs(m) = abs( wre(m) )
        index(m) = m
      enddo
*     Index the vector c2rabs, i.e. output the vector index(1:nl) such
*     that cr2abs( index(m) ) is in ascending order for m = 1,2,...,nl
      do m=2,nl
        indtmp = index(m)
        c2temp = c2rabs(indtmp)
*       Order m-th entry w.r.t. previous (already sorted) ones
        do i=m-1,1,-1
          if ( c2rabs(index(i)).le.c2temp ) goto 100
          index(i+1) = index(i)
        enddo
        i = 0
  100   index(i+1) = indtmp
      enddo
*     Modal eigenvalues (barotropic is ~0; replace with exact 0)
*     and wavespeeds (barotropic is infinite; replace with 0)
      eval(1) = 0.0d0
      cphs(1) = 0.0d0
      do m=2,nl
        eval(m) = c2rabs(index(m))
        cphs(m) = 1.0d0/sqrt( c2rabs(index(m)) )
      enddo
*     Deformation radii
      rdef(1) = 0.0d0
      rdm2(1) = 0.0d0
      do m=2,nl
        rdef(m) = 1.0d0/sqrt( c2rabs(index(m)) )/abs(fnot)
        rdm2(m) = fnot*fnot*c2rabs(index(m))
      enddo
*     Mode/layer conversion coefficients
*     m = mode number; k = layer number
      do m=1,nl
        inm = index(m)
        do k=1,nl
          cl2m(m,k) = evecl(k,inm)/elder(inm,inm)
          cm2l(k,m) = evecr(k,inm)
          ctl2m(k,m) = cl2m(m,k)
          ctm2l(m,k) = cm2l(k,m)
        enddo
      enddo
*     Check product of conversion matrices
      do l=1,nl
        do k=1,nl
          dotp = 0.0d0
          do m=1,nl
            dotp = dotp + cm2l(k,m)*cl2m(m,l)
          enddo
          ccprod(k,l) = dotp
        enddo
      enddo

*     Write formatted results to standard output
*     ------------------------------------------
      write (*,'(a)') ' '
      write (*,'(2x,a,a)') 'Eigenmode solver; case = ',case
      write (*,'(2x,a)') 'A matrix:'
      write (*,'(2x,a,i5,8i16)') '  i  j:  ',(j,j=1,nl)
      do i=1,nl
        write (*,'(2x,i3,1x,1p,9d16.7)') i,(aaa(i,j),j=1,nl)
      enddo
      write (*,'(a)') ' '
      write (*,'(2x,a)') 'Modes:'
      write (*,'(2x,a,i9,8i16)') ' m:  ',(m,m=1,nl)
      write (*,'(2x,a,1p,9d16.7)') 'eval',(eval(m),m=1,nl)
      write (*,'(2x,a,1p,9d16.7)')
     &         'Cphs     Infinite   ',(cphs(m),m=2,nl)
      write (*,'(2x,a,1p,9d16.7)')
     &         'Rdef     Infinite   ',(rdef(m),m=2,nl)
      write (*,'(2x,a,1p,9d16.7)') 'Rdm2',(rdm2(m),m=1,nl)
      write (*,'(2x,a)') '(eigenvalue = 1/Cphs(m)^2)'
      write (*,'(2x,a)') 'right eigenvectors Rm(k):'
      do k=1,nl
        write (*,'(2x,i3,1x,1p,9d16.7)') k,(evecr(k,index(m)),m=1,nl)
      enddo
      write (*,'(2x,a)') 'left  eigenvectors Lm(k):'
      do k=1,nl
        write (*,'(2x,i3,1x,1p,9d16.7)') k,(evecl(k,index(m)),m=1,nl)
      enddo
      write (*,'(a)') ' '
      write (*,'(2x,a)') 'layer to mode coefficients cl2m(m,k):'
      write (*,'(2x,a,i5,8i16)') '  m  k:  ',(k,k=1,nl)
      do m=1,nl
        write (*,'(2x,i3,1x,1p,9d16.7)') m,(cl2m(m,k),k=1,nl)
      enddo
      write (*,'(2x,a)') 'mode to layer coefficients cm2l(k,m):'
      write (*,'(2x,a,i5,8i16)') '  k  m:  ',(m,m=1,nl)
      do k=1,nl
        write (*,'(2x,i3,1x,1p,9d16.7)') k,(cm2l(k,m),m=1,nl)
      enddo
      write (*,'(2x,a)') 'matrix product cm2l(k,m)*cl2m(m,l):'
      write (*,'(2x,a,i5,8i16)') '  k  l:  ',(l,l=1,nl)
      do k=1,nl
        write (*,'(2x,i3,1x,1p,9d16.7)') k,(ccprod(k,l),l=1,nl)
      enddo
      write (*,'(2x,a)') '(should be the identity matrix)'
*
*     Optionally verify that the (ordered) eigenvectors are correct
*     -------------------------------------------------------------
      if ( eigchk ) then
        write (*,'(a)') ' '
        write (*,'(2x,a)') 'Verify ordered eigenmodes:'
        do m=1,nl
          inm = index(m)
          do i=1,nl
            dotp = 0.0d0
            do j=1,nl
              dotp = dotp + aaa(i,j)*evecr(j,inm)
            enddo
            aevec(i) = dotp
            eevec(i) = eval(m)*evecr(i,inm)
          enddo
          write (*,'(2x,a,i2)') 'For mode m = ',m
          write (*,'(2x,a,1p,9d17.9)') '  Rm   = ',(evecr(i,inm),i=1,nl)
          write (*,'(2x,a,1p,9d17.9)') ' A *Rm = ',(aevec(i),i=1,nl)
          write (*,'(2x,a,1p,9d17.9)') 'lam*Rm = ',(eevec(i),i=1,nl)
        enddo
      endif

*     Check matrix triple product occuring in cyclic constraint equations
*     -------------------------------------------------------------------
*     Matrix is Cl2m*A*Cm2l
**    do j=1,nl
*       Work out A*Cm2l for current j
**      do k=1,nl
**        ac(k) = 0.0d0
**        do m=1,nl
**          ac(k) = ac(k) + aaa(k,m)*cm2l(m,j)
**        enddo
**      enddo
*       Work out Cl2m*(A*Cm2l) for current i
**      do i=1,nl
**        ca = 0.0d0
**        do k=1,nl
**          ca = ca + cl2m(i,k)*ac(k)
**        enddo
**        ccprod(i,j) = ca
**      enddo
**    enddo
**    write (*,'(a)') ' '
**    write (*,'(2x,a,i5,8i16)') '  i  j:  ',(j,j=1,nl)
**    write (*,'(2x,a,1p,9d16.7)') 'eval',(eval(j),j=1,nl)
**    write (*,'(2x,a,a)') 'matrix triple product ',
**   &                     'cl2m(i,k)*A(k,m)*cm2l(m,j):'
**    do i=1,nl
**      write (*,'(2x,i3,1x,1p,9d16.7)') i,(ccprod(i,j),j=1,nl)
**    enddo
*     Confirms that the matrix triple product is a diagonal
*     matrix whose entries are the eigenvalues Lambda

      END SUBROUTINE eigmod
*
c***********************************************************************
c
      END MODULE eigmode
c
c***********************************************************************
