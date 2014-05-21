c***********************************************************************
c     Q-GCM Version 1.5.0 : last modified 26/08/2013
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
      MODULE radsubs

*     Contains subroutines radiat and trapin. radiat computes
*     and prints the values of various radiation coefficients,
*     derived as explained in Appendix A of the "Formulation
*     and users' guide for Q-GCM". trapin is used to compute
*     the integrals required for some of the coefficients.

*     Modules

      IMPLICIT NONE

      PRIVATE

      PUBLIC  :: radiat
      PRIVATE :: trapin

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
      USE parameters, ONLY : nyta, nla, nyto, nlo, fnot
      USE atconst, ONLY : yta, ytarel, hat, hta, tabsat, tmbara, tat
      USE occonst, ONLY : yto, ytorel, tabsoc, tmbaro, toc
      USE intrfac, ONLY : sstbar, astbar, hmat, tnbdy, tsbdy
      USE xfosubs, ONLY : fsprim
      USE radiate

      IMPLICIT NONE

*     Subroutine arguments
*
*     Local parameters
      double precision stefan,sigov2,tmbtol
      integer nz,nitmax
      logical prtval
*     stefan is the Stefan-Boltzmann constant (usually denoted by sigma)
      parameter ( stefan = 5.67040D-8, sigov2 = 0.5d0*stefan,
     &            nz = 10001, nitmax = 200, tmbtol = 1.0d-13,
     &            prtval = .false. )
*
*     Local variables
      integer i,j,k,l,iter
      double precision taum,tauk(nla),tupmul,hbot,htop,delz,zz,fup(nz),
     &                 fdn(nz),upint,dnint,uprad(nla),dnrad(nla),
     &                 rhstat,deltm,rhstoc,tocold,F0upbar,Fmupbar,
     &                 Fupbar(nla),Fmdnbar,Fdnbar(nla),rrcpdt
*     Variables for calculating radiation balance coefficients
      integer ipivrb(nla),info,iwork(nla)
      double precision rbalar(nla,nla),rballu(nla,nla),balrhs(nla),
     &                 rbafac(nla),ferr,berr,work(3*nla),rbtmat,
     &                 rbtmoc

*     Mean state radiation balance
*     ============================

*     Find layer transmissivities
*     ---------------------------
      taum = exp( -hmat/zm )
      tauk(1) = exp( -(hat(1)-hmat)/zopt(1) )
      tupmul = tauk(1)
      do k=2,nla
        tauk(k) = exp( -hat(k)/zopt(k) )
        tupmul = tupmul*tauk(k)
      enddo

*     Layer 1 up- and downgoing mean radiation
*     ----------------------------------------
      hbot = hmat
      htop = hat(1)
      delz = ( htop - hbot )/dble(nz-1)
!$OMP PARALLEL DO DEFAULT (NONE)
!$OMP&         PRIVATE (i,zz)
!$OMP&         SHARED  (hbot,htop,delz,fup,fdn,tabsat,gamma,zopt,zm)
!$OMP&         SCHEDULE (STATIC)
      do i=1,nz
        zz = hbot + dble(i-1)*delz
        fup(i) = ( tabsat(1)-gamma*zz )**4 * exp( -(htop - zz)/zopt(1) )
        fdn(i) = ( tabsat(1)-gamma*zz )**4 * exp( (hbot - zz)/zopt(1) )
      enddo
!$OMP END PARALLEL DO
      call trapin (fup, nz, delz, upint)
      call trapin (fdn, nz, delz, dnint)
*     upint = Iup(1) as defined in Appendix A
*     dnint = Idn(1) as defined in Appendix A
      uprad(1) = sigov2*upint/zopt(1)
      dnrad(1) = sigov2*dnint/zopt(1)
      rhstat = uprad(1)

*     Upper layers up- and downgoing mean radiation
*     ---------------------------------------------
      do k=2,nla
        hbot = htop
        htop = hbot + hat(k)
        delz = hat(k)/dble(nz-1)
!$OMP   PARALLEL DO DEFAULT (NONE)
!$OMP&           PRIVATE (i,zz)
!$OMP&           SHARED  (k,hbot,htop,delz,fup,fdn,tabsat,gamma,zopt)
!$OMP&           SCHEDULE (STATIC)
        do i=1,nz
          zz = hbot + dble(i-1)*delz
          fup(i) = ( tabsat(k)-gamma*zz )**4
     &             * exp( -(htop - zz)/zopt(k) )
          fdn(i) = ( tabsat(k)-gamma*zz )**4
     &             * exp( (hbot - zz)/zopt(k) )
        enddo
!$OMP   END PARALLEL DO
        call trapin (fup, nz, delz, upint)
        call trapin (fdn, nz, delz, dnint)
*       upint = Iup(k) as defined in Appendix A
*       dnint = Idn(k) as defined in Appendix A
        uprad(k) = sigov2*upint/zopt(k)
        dnrad(k) = sigov2*dnint/zopt(k)
        rhstat = rhstat*tauk(k) + uprad(k)
      enddo

*     Compute atmosphere m. l. mean temperature
*     -----------------------------------------
      rhstat = (-rhstat - fsbar)/tupmul
      rhstat = 2.0d0*zm*rhstat/stefan
**    print *,' rhstat = ',rhstat
*     rhstat now contains required value of
*     integral for mixed layer upgoing radiation
*     Adjust tmbara to get required value of integral
      tmbara = 300.0d0
      deltm = 0.0d0
      delz = hmat/dble(nz-1)
      iter = 0
  100 continue
!$OMP   PARALLEL DO DEFAULT (NONE)
!$OMP&           PRIVATE (i,zz)
!$OMP&           SHARED  (hmat,delz,fup,tmbara,gamma,zm)
!$OMP&           SCHEDULE (STATIC)
        do i=1,nz
          zz = dble(i-1)*delz
          fup(i) = ( tmbara-gamma*zz )**4 * exp( -(hmat - zz)/zm )
        enddo
!$OMP   END PARALLEL DO
        call trapin (fup, nz, delz, upint)
*       Convert error to (approximate) temperature change
        deltm = 0.25d0*(rhstat - upint)*tmbara/upint
        tmbara = tmbara + 0.75d0*deltm
        iter = iter + 1
**      print *,' iter, deltm, tmbara = ',iter,deltm,tmbara
        if ( iter.gt.nitmax ) then
          print *,' iteration for tmbara not converged'
          print *,' max. no. of iterations nitmax = ',nitmax
          print *,' deltm, tmbara = ',deltm,tmbara
          print *,' program terminates in radiat'
          stop
        endif
      if ( abs( deltm ).gt.tmbtol ) goto 100

*     Compute ocean m. l. mean temperature
*     ------------------------------------
      rhstoc = xlamda*tmbara + sigov2*tmbara**4 - fsbar
*     Use atmos. temp. as initial guess for ocean temp.
      tmbaro = tmbara
      iter = 0
  150 continue
        tocold = tmbaro
        tmbaro = rhstoc/( xlamda + stefan*tocold**3 )
        iter = iter + 1
**      print *,' iter, tmbaro = ',iter,tmbaro
        if ( iter.gt.nitmax ) then
          print *,' iteration for tmbaro not converged'
          print *,' max. no. of iterations nitmax = ',nitmax
          print *,' tocold, tmbaro = ',tocold,tmbaro
          print *,' program terminates in radiat'
          stop
        endif
      if ( abs( tmbaro - tocold ).gt.tmbtol ) goto 150

*     Derive temperature anomalies
      do k=1,nlo
        toc(k) = tabsoc(k) - tmbaro
      enddo
      do k=1,nla
        tat(k) = tabsat(k) - tmbara
      enddo

*     Upgoing mean state fluxes
*     -------------------------
*     Ocean
      F0upbar = stefan*tmbaro**4
*     Atmos. mixed layer
      Fmupbar = sigov2*upint/zm
*     Layer 1
      Fupbar(1) = Fmupbar*tauk(1) + uprad(1)
*     Upper layers
      do k=2,nla
        Fupbar(k) = Fupbar(k-1)*tauk(k) + uprad(k)
      enddo

*     Downgoing mean state fluxes
*     ---------------------------
*     Top layer
      Fdnbar(nla) = -dnrad(nla)
*     Other QG layers
      do k=nla-1,1,-1
        Fdnbar(k) = Fdnbar(k+1)*tauk(k) - dnrad(k)
      enddo
*     Atmos. mixed layer
      Fmdnbar = -sigov2*tmbara**4

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
      write(*,204) '  Optic. depth in aml. zm (m) = ',zm
      do k=1,nla
        write(*,224) '  Optic. depth (m) in layer',k,' = ',zopt(k)
      enddo
      write(*,214) '  Lapse rate gamma   (K m^-1) = ',gamma
      print *,' '
      write(*,253) '  A.m.l. transmissivity  taum = ',taum
      do k=1,nla
        write(*,273) '  Layer', k ,'  transmissivity tau = ',
     &               tauk(k)
      enddo
      write(*,253) '  Transmissivity prod. tupmul = ',tupmul
      do k=1,nla
        write(*,272) '  Layer', k ,' integs uprad, dnrad = ',
     &               uprad(k),dnrad(k)
      enddo
      print *,' '
      write(*,214) '  Tolerance for m.l. temp (K) = ',tmbtol
      write(*,251) '  Atmos. mixed layer T (K, C) = ',
     &             tmbara,tmbara-273.15d0
      write(*,251) '  Ocean  mixed layer T (K, C) = ',
     &             tmbaro,tmbaro-273.15d0
      write(*,207) '  Rel. ocean  layer temp. toc = ',(toc(k),k=1,nlo)
      write(*,207) '  Rel. atmos. layer temp. tat = ',(tat(k),k=1,nla)
      print *,' '
      write(*,251) '  Ocean m.l. radiat.  F0upbar = ',F0upbar
      write(*,251) '  Mixed layer  Fmbar up, down = ',Fmupbar,Fmdnbar
      do k=1,nla
        write(*,271) '  Layer', k ,'       Fbar up, down = ',
     &               Fupbar(k),Fdnbar(k)
      enddo
      write(*,214) '  Fractional error in OLR     = ',
     &             abs( Fupbar(nla) + fsbar )/abs(fsbar)

*     Perturbation (linearised) radiation parameters
*     ==============================================
*     Initialise Aup and Adown arrays
      do l=1,nla-1
        do k=1,nla
          Aup(k,l) = 0.0d0
          Adown(k,l) = 0.0d0
        enddo
      enddo

*     Upgoing
*     -------
*     Ocean mixed layer
      D0up = 4.0d0*stefan*tmbaro**3
*     Atmospheric mixed layer
      Bmup = ( sigov2*( tmbara-gamma*hmat )**4 - Fmupbar )/zm
      Cmup = Bmup
      delz = hmat/dble(nz-1)
!$OMP PARALLEL DO DEFAULT (NONE)
!$OMP&         PRIVATE (i,zz)
!$OMP&         SHARED  (hmat,delz,fup,tmbara,gamma,zm)
!$OMP&         SCHEDULE (STATIC)
      do i=1,nz
        zz = dble(i-1)*delz
        fup(i) = ( tmbara-gamma*zz )**3 * exp( -(hmat - zz)/zm )
      enddo
!$OMP END PARALLEL DO
      call trapin (fup, nz, delz, upint)
      Dmup = 2.0d0*stefan*upint/zm
*     Layer 1
      hbot = hmat
      htop = hat(1)
      Aup(1,1) = ( -tauk(1)*Fmupbar - uprad(1)
     &            + sigov2*( tabsat(1)-gamma*hat(1) )**4 )/zopt(1)
      Bup(1) = tauk(1)*( Bmup + Fmupbar/zopt(1) - sigov2*
     &                    ( tabsat(1)-gamma*hmat )**4/zopt(1) )
      Cup(1) = tauk(1)*( Cmup + Fmupbar/zopt(1) - sigov2*
     &                    ( tabsat(1)-gamma*hmat )**4/zopt(1) )
      Dup(1) = Dmup*tauk(1)
*     Upper layers
      do k=2,nla
        hbot = htop
        htop = hbot + hat(k)
        Bup(k) = Bup(k-1)*tauk(k)
        Cup(k) = Cup(k-1)*tauk(k)
        Dup(k) = Dup(k-1)*tauk(k)
        do l=1,k-2
          Aup(k,l) = Aup(k-1,l)*tauk(k)
        enddo
        Aup(k,k-1) = tauk(k)*( Aup(k-1,k-1) + Fupbar(k-1)/zopt(k)
     &                    - sigov2*(tabsat(k)-gamma*hbot)**4/zopt(k) )
*       Top layer has no eta coefft at top
        if ( k.lt.nla ) then
          Aup(k, k ) = ( -tauk(k)*Fupbar(k-1) - uprad(k)
     &                  + sigov2*(tabsat(k)-gamma*htop)**4 )/zopt(k)
        endif
      enddo

*     Downgoing
*     ---------
*     Top layer
      htop = hta
      hbot = htop - hat(nla)
      Adown(nla,nla-1) = (  sigov2*( tabsat(nla)-gamma*hbot )**4
     &                    - dnrad(nla) )/zopt(nla)
*     Intermediate layers
      do k=nla-1,2,-1
        htop = hbot
        hbot = htop - hat(k)
        do l=k+1,nla-1
          Adown(k,l) = Adown(k+1,l)*tauk(k)
        enddo
        Adown(k,k-1) = (  Fdnbar(k+1)*tauk(k) - dnrad(k)
     &                  + sigov2*(tabsat(k)-gamma*hbot)**4 )/zopt(k)
        Adown(k, k ) = tauk(k)*( Adown(k+1,k) - Fdnbar(k+1)/zopt(k)
     &                 - sigov2*( tabsat(k)-gamma*htop )**4/zopt(k) )
      enddo
*     Layer 1
      do l=2,nla-1
        Adown(1,l) = Adown(2,l)*tauk(1)
      enddo
      Adown(1,1) = tauk(1)*( Adown(2,1) - Fdnbar(2)/zopt(1) -
     &                 sigov2*( tabsat(1)-gamma*hat(1) )**4/zopt(1) )
      B1down = (  Fdnbar(2)*tauk(1) - dnrad(1)
     &          + sigov2*( tabsat(1)-gamma*hmat )**4 )/zopt(1)
      C1down = B1down
*     Atmospheric mixed layer
      Dmdown = - 2.0d0*stefan*tmbara**3

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

*     Radiation balance initialisation coefficients
*     =============================================
*     We assume the mixed layer interface displacement etam
*     and atmospheric topography are both zero, and
*     that the entrainments between QG layers are all zero
*     The first nla-1 entries in the solution vector are the
*     coefficients of the internal interface displacements;
*     the last entry is the coefficient of Tm'
*     Zero the matrices
      do i=1,nla
        do k=1,nla
          rbalar(k,i) = 0.0d0
          rballu(k,i) = 0.0d0
        enddo
      enddo
*     Compute the matrices from the radiation coefficients
*     1st equation is balance for layer 1; equation (D.4)
      do i=1,nla-1
        rbalar(1,i) = Adown(1,i)
      enddo
      rbalar(1,nla) = Dmup
*     Option 1: only e(1) is nonzero
      do k=2,nla-1
        do i=1,nla-1
          rbalar(k,i) = Adown(k+1,i) + Aup(k,i)
        enddo
        rbalar(k,nla) = Dup(k)
      enddo
*     Final equation is OLR at top of atmosphere
      do i=1,nla-1
        rbalar(nla,i) = Aup(nla,i)
      enddo
      rbalar(nla,nla) = Dup(nla)
*     Setup RHS & make copy for LU factorisation
      do k=1,nla
        balrhs(k) = -1.0d0
        rbafac(k) = balrhs(k)
        do i=1,nla
          rballu(k,i) = rbalar(k,i)
        enddo
      enddo
*     Compute the LU factorization of rbalar
*     DGETRF = NAG routine F07ADF
      call DGETRF (nla, nla, rballu, nla, ipivrb, info)
      if ( info.ne.0 ) then
        print *,'  DGETRF in radiat returns info = ',info
        print *,'  program terminates in radiat'
        stop
      endif
      print *,' '
      write (*,*) ' Radiation balance matrices:'
      print *,' '
      write (*,*) ' rbalar:'
      do k=1,nla
        write (*,'(2x,i2,1x,1p,5d17.9)') k,(rbalar(k,i),i=1,nla)
      enddo
      print *,' '
      write (*,*) ' rballu:'
      do k=1,nla
        write (*,'(2x,i2,1x,1p,5d17.9)') k,(rballu(k,i),i=1,nla)
      enddo

*     Solve matrix equation for radiative balance coeffts using LAPACK
*     Solve the linear system using the LU factorised matrix rballu
*     DGETRS = NAG routine F07AEF
      call DGETRS ('Norm', nla, 1, rballu, nla,
     &             ipivrb, rbafac, nla, info)
      if ( info.ne.0 ) then
        print *,'  DGETRS in radiat returns info = ',info
        print *,'  program terminates in radiat'
        stop
      endif
*     Improve the solution in rbafac by iterative refinement
*     DGERFS = NAG routine F07AHF
      call DGERFS ('Norm', nla, 1, rbalar, nla, rballu, nla,
     &             ipivrb, balrhs, nla, rbafac, nla, ferr, berr,
     &             work, iwork, info)
      if ( info.ne.0 ) then
        print *,'  DGERFS in radiat returns info = ',info
        print *,'  program terminates in radiat'
        stop
      endif
      do k=1,nla-1
        rbetat(k) = rbafac(k)
      enddo
      rbtmat = rbafac(nla)
      rbtmoc = ( (xlamda - Dmdown)*rbtmat - 1.0d0 )/( xlamda + D0up )

      print *,' '
      print *,' Radiative balance initialisation coeffts:'
      print *,' -----------------------------------------'
      print *,' QG internal interface eta coeffts:'
      do k=1,nla-1
        write(*,271) '  Layer', k ,' coefft       rbetat = ',rbetat(k)
      enddo
      write(*,251) '  Atmos. m.l. T coefft rbtmat = ',rbtmat
      write(*,251) '  Ocean  m.l. T coefft rbtmoc = ',rbtmoc

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
      tsbdy = sstbar(1)
      write(*,253) '  Rel. o.m.l. S. bndry temp.  = ',tsbdy

*     Entrainment factors: aface(l), bface, cface and dface such that
*     e(1) = Sum(l)[ aface(l)*eta(l) ] + bface*etam + cface*aD + dface*aTm'
*     We are assuming all entrainment is across interface 1.
      rrcpdt = rrcpat/( tat(2)-tat(1) )
      do l=1,nla-1
        aface(l) = rrcpdt*( Adown(1,l) - Aup(nla,l) )
      enddo
      bface = rrcpdt*( B1down + Bmup - Bup(nla) )
      cface = rrcpdt*( C1down + Cmup - Cup(nla) )
      dface = rrcpdt*( Dmup - Dup(nla) )
      print *,' '
      print *,' Entrainment coefficients:'
      print *,' -------------------------'
      print *,' QG internal interface 1 coeffts:'
      write(*,215) '  eta  coefficients  aface(l) = ',
     &             (aface(l),l=1,nla-1)
      write(*,260) '  etam   coefficient    bface = ',bface
      write(*,260) '  topog. coefficient    cface = ',cface
      write(*,260) '  aTm    coefficient    dface = ',dface

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
      SUBROUTINE trapin (fofz, nz, delz, vinteg)
*
*     Computes the extended trapezoidal rule approximation to
*     the integral of a tabulated function. fofz is a vector
*     of length nz containing the function tabulation (including
*     at the end points), delz is the tabulation interval,
*     and vinteg is the returned value of the integral.
*
*     This version modified to use Kahan summation (a.k.a.
*     compensated summation), in order to improve accuracy.
*     We use the variant of Kahan summation where the correction from the
*     last term is added on, in case any further benefit can be achieved.

      IMPLICIT NONE

*     Subroutine arguments
      integer nz
      double precision fofz(nz),delz,vinteg
*
*     Local variables
      integer i
      double precision sum,corr,yadd,sest
*
      sum = 0.5d0*fofz(1)
      corr = 0.0d0
      do i=2,nz-1
        yadd = fofz(i) - corr
        sest = sum + yadd
        corr = (sest - sum ) - yadd
        sum = sest
      enddo
      yadd = 0.5d0*fofz(nz) - corr
      sest = sum + yadd
*     Compute and apply 'last term' correction
      corr = ( sest - sum ) - yadd
      sum = sest - corr
      vinteg = delz*sum

      END SUBROUTINE trapin
c
c***********************************************************************
c
      END MODULE radsubs
c
c***********************************************************************
