c***********************************************************************
c     Q-GCM Version 1.5.0 : last modified 22/02/2013
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

*     Modules
      USE parameters
      USE atconst
      USE occonst
      USE atstate
      USE ocstate
      USE intrfac
      USE radiate
      USE radsubs

      IMPLICIT NONE

*     Local parameters
      double precision secday,daysyr,secsyr
      parameter ( secday=86400.0d0, daysyr=365.0d0,
     &            secsyr=secday*daysyr )
*     Define I/O unit numbers:
*     ipunit is for input.params
*     odunit is for output directory
      integer ipunit, odunit
      parameter ( ipunit=44, odunit=45 )
*
      DOUBLE PRECISION PI,TWOPI,PIBY2
      PARAMETER ( PI=3.14159265358979324D0,
     &            TWOPI=6.28318530717958648D0,
     &            PIBY2=1.57079632679489662D0 )
*
*     Local variables
      integer i,j,k,nsko,nska,nstr,lenod,lename,outfloc(7),outflat(7)
      double precision tini,trun,valday,odiday,adiday,dgnday,
     &                 resday,dtavat,dtavoc,dtcovat,dtcovoc
      character (len=72) :: name
      character (len=80) :: outdir, inpbuf
*     Extra variables for topography
      character (len=80) :: topocname,topatname

*     Extra variables for OpenMP
!$    integer nprocs,OMP_GET_NUM_PROCS,nthmax,OMP_GET_MAX_THREADS,
!$   &        numthr,OMP_GET_NUM_THREADS,thrnum,OMP_GET_THREAD_NUM
!$    logical dynadj,OMP_GET_DYNAMIC,nested,OMP_GET_NESTED

      print *,' '
      write(*,*) ' Computational parameters:'
      write(*,*) ' -------------------------'

*     Specify input parameters in included file
      INCLUDE './in_param.f'
      print *,' '
      lenod = index(outdir, '   ') - 1
      print *,' outdir = ',outdir(1:lenod)
      lename = index(name, '   ') - 1
      print *,' name = ',name(1:lename)

*     Specify configuration being used

      print *,' '
      print *,' Control parameters:'
      print *,' -------------------'

*     Derive atmos gridspace and ocean timestep
*     Derive larger from smaller to ensure integer ratio
      dxa = ndxr*dxo
      dto = nstr*dta
      write(*,201) '  Atmos/ocean grid ratio ndxr = ',ndxr
      write(*,201) '  Oc/atm. timestep ratio nstr = ',nstr
      write(*,201) '  Atmos. gridcells over ocean = ',nxaooc,nyaooc
      write(*,201) '  Ocn start indices  nx1, ny1 = ',nx1,ny1
      write(*,215) '  Coriolis param.  (rad s^-1) = ',fnot
      write(*,215) '  Beta =df/dy (rad s^-1 m^-1) = ',beta

*     Atmospheric grid
*     ----------------
      dya = dxa
      hdxam1 = 0.5d0/dxa
      dxam2 = 1.0d0/(dxa*dxa)
      xla = nxta*dxa
      yla = nyta*dya
      do i=1,nxpa
        xpa(i) = (i-1)*dxa
      enddo
      do i=1,nxta
        xta(i) = xpa(i) + 0.5d0*dxa
      enddo
      do j=1,nypa
        ypa(j) = (j-1)*dya
        yparel(j) = ypa(j) - 0.5d0*yla
      enddo
      do j=1,nyta
        yta(j) = ypa(j) + 0.5d0*dya
        ytarel(j) = yta(j) - 0.5d0*yla
      enddo

*     Oceanic grid
*     ------------
      dyo = dxo
      hdxom1 = 0.5d0/dxo
      dxom2 = 1.0d0/(dxo*dxo)
      xlo = nxto*dxo
      ylo = nyto*dyo
      do i=1,nxpo
        xpo(i) = (i-1)*dxo + (nx1-1)*dxa
      enddo
      do i=1,nxto
        xto(i) = xpo(i) + 0.5d0*dxo
      enddo
      do j=1,nypo
        ypo(j) = (ny1-1)*dya + (j-1)*dyo
        yporel(j) = ypo(j) - 0.5d0*yla
      enddo
      do j=1,nyto
        yto(j) = ypo(j) + 0.5d0*dyo
        ytorel(j) = yto(j) - 0.5d0*yla
      enddo

*     Derive useful constants
*     -----------------------
      rdxaf0 = 1.0d0/(dxa*fnot)
      rdxof0 = 1.0d0/(dxo*fnot)
      rrcpat = 1.0d0/(rhoat*cpat)
      rrcpoc = 1.0d0/(rhooc*cpoc)
      raoro = rhoat/rhooc
      tdto = 2.0d0*dto
      tdta = 2.0d0*dta

*     Derive total thicknesses
*     ------------------------
      hto = 0.0d0
      do k=1,nlo
        hto = hto + hoc(k)
      enddo
      hta = 0.0d0
      do k=1,nla
        hta = hta + hat(k)
      enddo

*     Print out a few interesting numbers
*     -----------------------------------
      print *,' '
      print *,' Oceanic parameters:'
      print *,' -------------------'
      write(*,201) '  No. of ocean QG layers  nlo = ',nlo
      write(*,201) '  No. of gridcells nxto, nyto = ',nxto,nyto
      write(*,204) '  Gridlength dxo         (km) = ',1.0d-3*dxo
      write(*,203) '  Domain sizes xlo, ylo  (km) = ',
     &             1.0d-3*xlo,1.0d-3*ylo
      write(*,205) '  Timestep dto      (minutes) = ',dto/60.0d0
      write(*,201) '  No. of timesteps per day    = ',nint(secday/dto)
      write(*,203) '  Mixed layer thickness   (m) = ',hmoc
      write(*,213) '  Sp. ht. cap. (J kg^-1 K^-1) = ',cpoc
      write(*,205) '  Abs. pot. temp. tabsoc  (K) = ',
     &             (tabsoc(k),k=1,nlo)
      write(*,203) '  Layer thicknesses hoc   (m) = ',(hoc(k),k=1,nlo)
      write(*,207) '  Reduced gravities  (m s^-2) = ',
     &             (gpoc(k),k=1,nlo-1)

      print *,' '
      print *,' Atmospheric parameters:'
      print *,' -----------------------'
      write(*,201) '  No. of atmos. QG layers nla = ',nla
      write(*,201) '  No. of gridcells nxta, nyta = ',nxta,nyta
      write(*,204) '  Gridlength dxa         (km) = ',1.0d-3*dxa
      write(*,203) '  Domain sizes xla, yla  (km) = ',
     &             1.0d-3*xla,1.0d-3*yla
      write(*,205) '  Timestep dta      (minutes) = ',dta/60.0d0
      write(*,201) '  No. of timesteps per day    = ',nint(secday/dta)
      write(*,203) '  Mixed layer thickness   (m) = ',hmat
      write(*,203) '  Min. mixed layer thick. (m) = ',hmamin
      write(*,213) '  Sp. ht. cap. (J kg^-1 K^-1) = ',cpat
      write(*,203) '  Layer thicknesses hat   (m) = ',(hat(k),k=1,nla)
      write(*,203) '  Total thickness   hta   (m) = ',hta
      write(*,205) '  Abs. pot. temp. tabsat  (K) = ',
     &             (tabsat(k),k=1,nla)
      write(*,207) '  Reduced gravities  (m s^-2) = ',
     &             (gpat(k),k=1,nla-1)

      print *,' '
      print *,' Coupling parameters:'
      print *,' --------------------'
      write(*,205) '  Coefft. Lambda   (W m^-2/K) = ',xlamda
      write(*,204) '  Ast coupling  coefft  xcexp = ',xcexp
      write(*,204) '  Sst advection coefft  ycexp = ',ycexp

  201 format(a,9i13)
  203 format(a,9f13.3)
  204 format(a,9f13.4)
  205 format(a,9f13.5)
  206 format(a,9f13.6)
  207 format(a,9f13.7)
  213 format(a,1p,9d13.3)
  214 format(a,1p,9d13.4)
  215 format(a,1p,9d13.5)
  224 format(a,i2,a,9f13.4)
  225 format(a,i2,a,9f13.5)
  226 format(a,i2,a,9f13.6)
  227 format(a,i2,a,9f13.7)

*     Radiation section
*     =================
*     Compute mean state radiative balance and perturbation
*     radiation coefficients A, B, C and D. Also compute atmosphere
*     and ocean mixed layer temperatures that ensure equilibrium
      call radiat

*     Initialise pressure and temperature fields
*     ==========================================
***   if ( name.eq.'zero' ) then
***     call zeroin
***     tini = 0.0d0
***    else if ( name.eq.'rbal' ) then
        call rbalin
        tini = 0.0d0
***    else
***     print *,' initial state readin removed'
***     print *,' program terminates'
***     stop
***   endif

      stop
      end
c
c***********************************************************************
c
      SUBROUTINE ipbget (buffer, iounit)
*
*     Reads records from unit iounit until a valid one (i.e.
*     one not marked with a "comment" character) is found,
*     then returns this valid character string for processing.
*     The comment marker in Q-GCM is deemed to be
*     an exclamation mark "!" in the first column.
*
      IMPLICIT NONE
*
*     Subroutine arguments
      character (len=80) :: buffer
      integer iounit
*
*     Local variables

  100 continue
      read (iounit, err=200, fmt='(a80)') buffer
      if ( buffer(1:1).eq.'!' ) goto 100
      return

  200 continue
      print *,' Error reading character buffer from iounit = ',iounit
      print *,' Program terminates in ipbget'
      stop

      END SUBROUTINE ipbget
c
c***********************************************************************
c
      SUBROUTINE zeroin
*
*     Set initial state to zero pressure, and radiative
*     equilibrium with unperturbed (background) forcing,
*     i.e. mixed layer temperature anomalies are all zero.

*     Modules
      USE parameters
      USE atstate
      USE ocstate
      USE intrfac, ONLY : sst, sstm, ast, astm, hmixa, hmixam, hmat

      IMPLICIT NONE

*     Subroutine arguments
*
*     Local variables
      integer i,j,k

!$OMP PARALLEL DEFAULT (NONE)
!$OMP&         PRIVATE (i,j,k)
!$OMP&         SHARED  (sst,sstm,ast,astm,hmixa,hmixam,hmat)
!$OMP&         SHARED  (po,pom)
!$OMP&         SHARED  (pa,pam)

*     Initialise ocean pressure
      do k=1,nlo
!$OMP   DO SCHEDULE (STATIC)
        do j=1,nypo
          do i=1,nxpo
            po(i,j,k) = 0.0d0
            pom(i,j,k) = 0.0d0
          enddo
        enddo
!$OMP   END DO NOWAIT
      enddo

*     Initialise atmospheric pressure
      do k=1,nla
!$OMP   DO SCHEDULE (STATIC)
        do j=1,nypa
          do i=1,nxpa
            pa(i,j,k) = 0.0d0
            pam(i,j,k) = 0.0d0
          enddo
        enddo
!$OMP   END DO NOWAIT
      enddo

*     Initialise oceanic mixed layer temp.
!$OMP DO SCHEDULE (STATIC)
      do j=1,nyto
        do i=1,nxto
          sst(i,j) = 0.0d0
          sstm(i,j) = 0.0d0
        enddo
      enddo
!$OMP END DO NOWAIT

*     Initialise atmos. mixed layer rel. temp. and thickness
!$OMP DO SCHEDULE (STATIC)
      do j=1,nyta
        do i=1,nxta
          ast(i,j) = 0.0d0
          astm(i,j) = 0.0d0
          hmixa(i,j) = hmat
          hmixam(i,j) = hmat
        enddo
      enddo
!$OMP END DO NOWAIT

!$OMP END PARALLEL
      print *,' '
      write (*,*) ' Zero initialisation complete'

      END SUBROUTINE zeroin
c
c***********************************************************************
c
      SUBROUTINE rbalin
*
*     Initialise in radiative equilibrium, with
*     radiative forcing perturbation fsprim included.

*     Modules
      USE parameters
      USE atconst
      USE occonst
      USE atstate
      USE ocstate
      USE intrfac, ONLY : sst, sstm, sstbar, ast, astm, astbar,
     &                    hmixa, hmixam, hmat
      USE radiate, ONLY : rbetat

      IMPLICIT NONE

*     Subroutine arguments
*
*     Local parameters
      logical prtval
      parameter ( prtval = .false. )
*
*     Local variables
      integer i,j,k
      double precision fsprim,plfac(nla),play
**    double precision sumc,sumcp,p1off

      print *,' '
      write (*,*) ' Radiative balance initialisation:'
      write (*,*) ' ---------------------------------'

*     Derive suitable multiplier of Fs' for each atmos. layer,
*     from the eta coefficients of Fs' derived in radiat.
*     We have nla layers but only nla-1 eta coeffts,
*     and so need an extra constraint.
      plfac(1) = 0.0d0
      do k=2,nla
        plfac(k) = plfac(k-1) - gpat(k-1)*rbetat(k-1)
      enddo

*     Option 1: leave alone. This gives p(1) = 0 everywhere;
*     no pressure gradient and thus no flow in layer 1.

*     Option 2: apply offset so that barotropic p = 0, as in zeroin.
**    sumc = 0.0d0
**    sumcp = 0.0d0
**    do k=1,nla
**      sumc = sumc + ctl2mat(k,1)
**      sumcp = sumcp + ctl2mat(k,1)*plfac(k)
**    enddo
**    p1off = -sumcp/sumc
**    sumcp = 0.0d0
**    do k=1,nla
**      plfac(k) = plfac(k) + p1off
**      sumcp = sumcp + ctl2mat(k,1)*plfac(k)
**    enddo
**    print *,' Barotropic coefft. in rbalin = ',sumcp
      write(*,206) '  Layer coeffts for pa, plfac = ',
     &             (plfac(k),k=1,nla)
  206 format(a,9f13.6)

      if ( prtval ) then
        print *,' '
        write (*,*) ' Initial relative sst:'
        do j=nyto,1,-1
          if ( sstbar(j).lt.toc(1) ) then
            write (*,'(i6,f16.8,a)') j,sstbar(j),'  convect'
           else
            write (*,'(i6,f16.8,a)') j,sstbar(j)
          endif
        enddo
      endif

      if ( prtval ) then
        print *,' '
        write (*,*) ' Initial relative ast:'
        do j=nyta,1,-1
          if ( astbar(j).gt.tat(1) ) then
            write (*,'(i6,f16.8,a)') j,astbar(j),'  convect'
           else
            write (*,'(i6,f16.8,a)') j,astbar(j)
          endif
        enddo
      endif

!$OMP PARALLEL DEFAULT (NONE)
!$OMP&         PRIVATE (i,j,k,play)
!$OMP&         SHARED  (sst,sstm,sstbar,ast,astm,
!$OMP&                  astbar,hmixa,hmixam,hmat)
!$OMP&         SHARED  (po,pom)
!$OMP&         SHARED  (pa,pam,plfac,yparel)

*     Initialise ocean pressure
      do k=1,nlo
!$OMP   DO SCHEDULE (STATIC)
        do j=1,nypo
          do i=1,nxpo
            po(i,j,k) = 0.0d0
            pom(i,j,k) = 0.0d0
          enddo
        enddo
!$OMP   END DO NOWAIT
      enddo

*     Initialise atmospheric pressure
      do k=1,nla
!$OMP   DO SCHEDULE (STATIC)
        do j=1,nypa
          play = plfac(k)*fsprim( yparel(j) )
          do i=1,nxpa
            pa(i,j,k) = play
            pam(i,j,k) = play
          enddo
        enddo
!$OMP   END DO NOWAIT
      enddo

*     Initialise oceanic mixed layer rel. temp.
!$OMP DO SCHEDULE (STATIC)
      do j=1,nyto
        do i=1,nxto
          sst(i,j) = sstbar(j)
          sstm(i,j) = sstbar(j)
        enddo
      enddo
!$OMP END DO NOWAIT

*     Initialise atmos. mixed layer rel. temp. and thickness
!$OMP DO SCHEDULE (STATIC)
      do j=1,nyta
        do i=1,nxta
          ast(i,j) = astbar(j)
          astm(i,j) = astbar(j)
          hmixa(i,j) = hmat
          hmixam(i,j) = hmat
        enddo
      enddo
!$OMP END DO NOWAIT

!$OMP END PARALLEL
      print *,' '
      write (*,*) ' Rbal initialisation complete'

      END SUBROUTINE rbalin
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
