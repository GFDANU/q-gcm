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
*
*     Simplified version of q-gcm.F initially for testing of subroutine
*     eigmod, now used to test the implications of values in input.params

*     Modules
      USE parameters
      USE atconst
      USE occonst
      USE intrfac
      USE radiate
      USE eigmode

      IMPLICIT NONE

*     Local parameters
      character (len=*), parameter :: subnam = 'eigtest_main'
      double precision secday
      parameter ( secday=86400.0d0 )
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
      double precision tspind
      integer i,j,k,m,nsko,nska,nstr,outfloc(7),outflat(7)
      double precision trun,valday,odiday,adiday,dgnday,resday,
     &                 dtavat,dtavoc,dtcovat,dtcovoc
      character (len=72) :: name
      character (len=80) :: outdir, inpbuf
*     Extra variables for topography
      character (len=80) :: topocname,topatname

*     Specify input parameters in included file
      INCLUDE './in_param.f'
      print *,' '

      print *,' '
      write(*,*) ' Computational parameters:'
      write(*,*) ' ========================='
      write(*,*) ' Model version is Q-GCM v1.5.0'
*     Check some of the grid parameters
      if ( nla.lt.2 .or. nlo.lt.2 ) then
        print *,' nla, nlo = ',nla,nlo
        print *,' Inadequate nla or nlo, needs to be at least 2'
        print *,' Program terminates'
        stop
      endif

*     Derive atmos gridspace and ocean timestep
*     -----------------------------------------
*     Derive larger from smaller to ensure integer ratio
      dxa = ndxr*dxo
      dto = nstr*dta
      write(*,201) '  Atmos/ocean grid ratio ndxr = ',ndxr
      write(*,201) '  Oc/atm. timestep ratio nstr = ',nstr
      write(*,201) '  Atmos. gridcells over ocean = ',nxaooc,nyaooc
      write(*,201) '  Ocn start indices  nx1, ny1 = ',nx1,ny1
      write(*,215) '  Coriolis par. f0 (rad s^-1) = ',fnot
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
      call eigmod (nlo, gpoc, hoc, 'Ocean', amatoc,
     &             cphsoc, rdefoc, rdm2oc, ctl2moc, ctm2loc)
      print *,' '
      print *,' Oceanic parameters:'
      print *,' -------------------'
      write(*,201) '  No. of ocean QG layers  nlo = ',nlo
      write(*,201) '  No. of gridcells nxto, nyto = ',nxto,nyto
      write(*,204) '  Gridlength dxo         (km) = ',1.0d-3*dxo
      write(*,203) '  Domain sizes xlo, ylo  (km) = ',
     &             1.0d-3*xlo,1.0d-3*ylo
      write(*,206) '  Rossby number   Beta*ylo/f0 = ',beta*ylo/abs(fnot)
      write(*,215) '  f range S -> N   (rad s^-1) = ',
     &             fnot+beta*yporel(1),fnot+beta*yporel(nypo)
      write(*,215) '  Midlatitude Coriolis param  = ',
     &             fnot+beta*0.5d0*( yporel(1) + yporel(nypo) )
      write(*,205) '  Timestep dto      (minutes) = ',dto/60.0d0
      write(*,201) '  No. of timesteps per day    = ',nint(secday/dto)
      write(*,203) '  Mixed layer thickness   (m) = ',hmoc
      write(*,213) '  SST grad-2 diff  (m^2 s^-1) = ',st2d
      call diffts (2, nlo, st2d, 1, dxo, rdefoc, dto)
      write(*,213) '  SST grad-4 diff  (m^4 s^-1) = ',st4d
      call diffts (4, nlo, st4d, 1, dxo, rdefoc, dto)
      write(*,203) '  Layer thicknesses hoc   (m) = ',(hoc(k),k=1,nlo)
      write(*,203) '  Total thickness   hto   (m) = ',hto
      write(*,207) '  Reduced gravities  (m s^-2) = ',
     &             (gpoc(k),k=1,nlo-1)
      write(*,206) '  Baroclinic wavespeeds (m/s) = ',
     &             (cphsoc(k),k=2,nlo)
      write(*,206) '  Courant number(s)           = ',
     &             ( (dto/dxo)*cphsoc(k),k=2,nlo)
      write(*,205) '  Deformation radii      (km) = ',
     &             (1.0d-3*rdefoc(k),k=2,nlo)
      write(*,206) '  Gridlengths per def. radius = ',
     &             (rdefoc(k)/dxo,k=2,nlo)
      write(*,207) '  Long Rossby wavespeed (m/s) = ',
     &             (-beta*cphsoc(k)**2/fnot**2,k=2,nlo)
      write(*,213) '  Del-sqd coeffts  (m^2 s^-1) = ',(ah2oc(k),k=1,nlo)
      call diffts (2, nlo, ah2oc, nlo, dxo, rdefoc, dto)
      write(*,213) '  Del-4th coeffts  (m^4 s^-1) = ',(ah4oc(k),k=1,nlo)
      call diffts (4, nlo, ah4oc, nlo, dxo, rdefoc, dto)
      write(*,204) '  Munk b.l. width scale  (km) = ',
     &             (1.0d-3*(ah4oc(k)/beta)**0.2d0,k=1,nlo)
      write(*,204) '  Munk b.l. width scale (pts) = ',
     &             (((ah4oc(k)/beta)**0.2d0)/dxo,k=1,nlo)
      write(*,204) '  Bottom Ekm. layer thickness = ',delek
      write(*,213) '  Bottom layer Ekman number   = ',
     &             (delek/hoc(nlo))**2
      if ( delek.lt.0.0d0 ) then
        print *,' Invalid -ve value of delek'
        print *,' Program terminates'
        stop
       else if (delek.eq.0.0d0 ) then
        tspind = 0.0d0
       else
        tspind = 2.0d0*hoc(nlo)/(abs(fnot)*delek)/secday
      endif
      write(*,204) '  Spindown timescale   (days) = ',tspind

      call eigmod (nla, gpat, hat, 'Atmosphere', amatat,
     &             cphsat, rdefat, rdm2at, ctl2mat, ctm2lat)
      print *,' '
      print *,' Atmospheric parameters:'
      print *,' -----------------------'
      write(*,201) '  No. of atmos. QG layers nla = ',nla
      write(*,201) '  No. of gridcells nxta, nyta = ',nxta,nyta
      write(*,204) '  Gridlength dxa         (km) = ',1.0d-3*dxa
      write(*,203) '  Domain sizes xla, yla  (km) = ',
     &             1.0d-3*xla,1.0d-3*yla
      write(*,206) '  Rossby number   Beta*yla/f0 = ',beta*yla/abs(fnot)
      write(*,215) '  f range S -> N   (rad s^-1) = ',
     &             fnot+beta*yparel(1),fnot+beta*yparel(nypa)
      write(*,215) '  Midlatitude Coriolis param  = ',
     &             fnot+beta*0.5d0*( yparel(1) + yparel(nypa) )
      write(*,205) '  Timestep dta      (minutes) = ',dta/60.0d0
      write(*,201) '  No. of timesteps per day    = ',nint(secday/dta)
      write(*,203) '  Mixed layer thickness   (m) = ',hmat
      write(*,213) '  AST grad-2 diff  (m^2 s^-1) = ',at2d
      call diffts (2, nla, at2d, 1, dxa, rdefat, dta)
      write(*,213) '  AST grad-4 diff  (m^4 s^-1) = ',at4d
      call diffts (4, nla, at4d, 1, dxa, rdefat, dta)
      write(*,213) '  hmix diffusivity (m^2 s^-1) = ',ahmd
      call diffts (2, nla, ahmd, 1, dxa, rdefat, dta)
      write(*,213) '  hmix damping coefft  hmadmp = ',hmadmp
      write(*,203) '  Layer thicknesses hat   (m) = ',(hat(k),k=1,nla)
      write(*,203) '  Total thickness   hta   (m) = ',hta
      write(*,205) '  Abs. pot. temp. tabsat  (K) = ',
     &             (tabsat(k),k=1,nla)
      write(*,207) '  Reduced gravities  (m s^-2) = ',
     &             (gpat(k),k=1,nla-1)
      write(*,206) '  Baroclinic wavespeeds (m/s) = ',
     &             (cphsat(k),k=2,nla)
      write(*,206) '  Courant number(s)           = ',
     &             ( (dta/dxa)*cphsat(k),k=2,nla)
      write(*,205) '  Deformation radii      (km) = ',
     &             (1.0d-3*rdefat(k),k=2,nla)
      write(*,206) '  Gridlengths per def. radius = ',
     &             (rdefat(k)/dxa,k=2,nla)
      write(*,213) '  Del-4th coeffts  (m^4 s^-1) = ',(ah4at(k),k=1,nla)
      call diffts (4, nla, ah4at, nla, dxa, rdefat, dta)

  201 format(a,9i13)
  202 format(a,9i7)
  203 format(a,9f13.3)
  204 format(a,9f13.4)
  205 format(a,9f13.5)
  206 format(a,9f13.6)
  207 format(a,9f13.7)
  213 format(a,1p,9d13.3)
  214 format(a,1p,9d13.4)
  215 format(a,1p,9d13.5)
  225 format(a,i2,a,9f13.5)

      stop
      end

c***********************************************************************
c
      SUBROUTINE ipbget (buffer, iounit)
*
*     Reads records from unit iounit until a valid one (i.e. one
*     not marked with a "comment" character) is found, then
*     returns this valid character string for processing.
*     The comment marker in Q-GCM is deemed to be
*     an exclamation mark "!" in the first column.

*     Modules
      USE parameters

      IMPLICIT NONE

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
      SUBROUTINE diffts (nord, nl, coeff, ncoef, dx, rdef, dt)

*     Computes diffusive decay timescales for circular eddies whose radii
*     are the baroclinic Rossby radii, and for two-gridpoint noise.
*     See section 8.6 of the Userguide for derivation of timescales.
*     Also computes the nondimensional timestep stability criterion.

*     Input arguments:
*     nord  : order of the diffusive term
*     nl    : no. of QG layers  (=> nl-1 baroclinic modes)
*     coeff : vector of diffusion coefficients (should be .GE. 0)
*     ncoef : length of coefficient vector
*     dx    : gridlength (m)
*     rdef  : vector of nl modal deformation radii (m)
*             (infinite value for barotropic mode replaced by 0.0)
*     dt    : timestep (sec)
*     (all the above are unchanged on exit)

*     Modules
      USE parameters

      IMPLICIT NONE

*     Subroutine arguments
      integer, INTENT(IN) :: nord,nl,ncoef
      double precision, INTENT(IN) :: coeff(ncoef),dx,rdef(nl),dt
*
*     Local parameters
      DOUBLE PRECISION PIBY2
      PARAMETER ( PIBY2=1.57079632679489662D0 )
      integer nlmax
      double precision secday
      parameter ( nlmax=9, secday=86400.0d0 )
*
*     Local variables
      integer k,m
      double precision tdamp(nlmax),sinfac,dtstab(nlmax)

*     Check internal storage is sufficient
      if ( nl.gt.nlmax ) then
        print *,' diffts has insufficient nlmax = ',nlmax
        print *,' called with nl = ',nl
        print *,' program terminates in diffts'
        stop
      endif

*     Check all diffusion coefficients are non-negative
*     (need positive coeffts for damping)
      do k=1,ncoef
        if ( coeff(k).lt.0.0d0 ) then
          print *,' diffts called with -ve diffusion coefft'
          print *,' coeff vector = ',(coeff(m),m=1,ncoef)
          print *,' program terminates in diffts'
          stop
        endif
      enddo

*     Compute decay timescale(s) for a circular eddy
*     at the deformation radius for each baroclinic mode
*     Deformation radii are usually much greater than
*     the gridlength dx, but for very coarse resolution
*     cases where rdef(m) < dx, use dx instead of
*     rdef(m) to compute the modal decay timescale
      do m=2,nl
        if ( dx.le.rdef(m) ) then
          sinfac = 2.0d0*sin( PIBY2*dx/rdef(m) )/dx
         else
          sinfac = 2.0d0/dx
          write(*,225) '  NOTE: mode',m-1,
     &                 ' defrad < dx; dx used instead'
          write(*,205) '  of defrad for computing damping timescale'
        endif
*       Avoid infinities if coefft = 0
        do k=1,ncoef
          if ( coeff(k).eq.0.0d0 ) then
            tdamp(k) = 0.0d0
           else
            tdamp(k) = 1.0d0/( sinfac**nord*coeff(k)*dble(nord)*secday )
          endif
        enddo
        write(*,225) '  Mode',m-1,' damping time  (days) = ',
     &               (tdamp(k),k=1,ncoef)
      enddo

*     Compute decay timescale for 1-D two-gridpoint noise
*     for each coefft, avoiding infinities if coefft = 0
*     The timescale for 2-D (checkerboard) noise differs only
*     by a factor of nord. This has the shortest decay timescale,
*     and thus determines the overall timestep stability limit.
*     Also compute nondimensional timestep stability factor.
*     This should be < 1 for stable leapfrog timestepping.
      do k=1,ncoef
        if ( coeff(k).eq.0.0d0 ) then
          tdamp(k) = 0.0d0
          dtstab(k) = 0.0d0
         else
          tdamp(k) = (0.5d0*dx)**nord/coeff(k)
          dtstab(k) = nord*dt/tdamp(k)
*         Rescale damping timescale to hours
          tdamp(k) = tdamp(k)/3600.0d0
        endif
      enddo
      write(*,205) '  1-D grid timescale  (hours) = ',
     &             (tdamp(k),k=1,ncoef)
      write(*,205) '  2-D grid timescale  (hours) = ',
     &             (tdamp(k)/nord,k=1,ncoef)
      write(*,205) '  Timestep stability factor   = ',
     &             (dtstab(k),k=1,ncoef)

  205 format(a,9f13.5)
  225 format(a,i2,a,9f13.5)

      END SUBROUTINE diffts
c
c***********************************************************************
