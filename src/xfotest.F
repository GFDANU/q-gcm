c***********************************************************************
c     Q-GCM Version 1.5.0 : last modified 20/09/2013
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
*     Q-GCM v1.5.0
*
*     Simplified version of q-gcm.F for use in testing xfosubs.F
*     Code following first call to xforc (almost all) removed.

*     Modules
      USE parameters
      USE atconst
      USE occonst
      USE timinfo
      USE radiate
      USE monitor
#if defined tau_udiff || !defined ocean_only
      USE athomog
      USE atstate
#endif
#ifndef atmos_only
      USE ochomog
      USE ocstate
#endif
      USE intrfac
      USE eigmode
      USE radsubs
      USE topsubs
#ifdef use_netcdf
      USE nc_subs
#endif /* use_netcdf */

      IMPLICIT NONE

#ifdef use_netcdf
      INCLUDE 'netcdf.inc'
#endif /* use_netcdf */
*
*     Local parameters
      character (len=*), parameter :: subnam = 'xfotest_main'
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
*     Arrays passed as arguments
*
*     Local variables
      integer i,j,k,iorest,nsko,nska,nt,lenod,lename,
     &        oceanonly,atmosonly,getcovar,outfloc(7),
     &        outflat(7),cyclicoc,hflxsb,hflxnb,tauudiff
      character (len=72) :: name
      character (len-80) :: outdir, inpbuf
      logical resdmp,solnok
*     Extra variables for topography
      character (len=80) :: topocname,topatname

#ifdef use_netcdf
*     Extra variables for netCDF input/output
      integer status
#else
#  ifdef atmos_only
      character (len=*), parameter :: atforav='./atforce.avg.binary'
#  endif
#  ifdef ocean_only
      character (len=*), parameter :: ocforav='./ocforce.avg.binary'
#  endif
      logical fexist
#endif

*     Extra variables for OpenMP
!$    integer nprocs,OMP_GET_NUM_PROCS,nthmax,OMP_GET_MAX_THREADS
!$    logical dynadj,OMP_GET_DYNAMIC,nested,OMP_GET_NESTED

      print *,' '
      write(*,*) ' Computational parameters:'
      write(*,*) ' ========================='
      write(*,*) ' Model version is Q-GCM v1.5.0'

*     Examine OpenMP environment
!$    nprocs = OMP_GET_NUM_PROCS()
!$    nthmax = OMP_GET_MAX_THREADS()
!$    dynadj = OMP_GET_DYNAMIC()
!$    nested = OMP_GET_NESTED()
!$    write(*,*) ' '
!$    write(*,*) ' OpenMP parallelism activated'
!$    write(*,'(a,i5)') '  No. of processors available = ',nprocs
!$    write(*,'(a,i3)') '  Max. no. of threads available = ',nthmax
!$    write(*,*) ' Dynamic adjustment = ',dynadj
!$    write(*,*) ' Nested parallelism = ',nested

*     Read input parameters via included file
      INCLUDE './in_param.f'

*     Specify configuration being used
*     --------------------------------
*     Specify as a set of integers so that the config-
*     uration information can be conveniently written
*     to a Matlab-compatible dump via out_param.f
#ifdef ocean_only
      oceanonly = 1
#else
      oceanonly = 0
#endif
#ifdef atmos_only
      atmosonly = 1
#else
      atmosonly = 0
#endif
      getcovar = 0
#ifdef cyclic_ocean
      cyclicoc = 1
#else
      cyclicoc = 0
#endif
#ifdef sb_hflux
      hflxsb = 1
#else
      hflxsb = 0
#endif
#ifdef nb_hflux
      hflxnb = 1
#else
      hflxnb = 0
#endif
#ifdef tau_udiff
      tauudiff = 1
#else
      tauudiff = 0
#endif

      print *,' '
      print *,' Control parameters:'
      print *,' ==================='
      lenod = index(outdir, '   ') - 1
      lename = index(name, '   ') - 1
      print *,' outdir = ',outdir(1:lenod)
      print *,' name = ',name(1:lename)
      if ( oceanonly.eq.1 .and. atmosonly.eq.1 ) then
        print *,' '
        print *,' Invalid model configuration: ocean_only and'
        print *,' atmos_only options cannot both be selected'
        print *,' Program terminates'
        stop
       else if ( oceanonly.eq.1 ) then
        print *,' Model running in ocean_only configuration'
       else if ( atmosonly.eq.1 ) then
        print *,' Model running in atmos_only configuration'
       else
        print *,' Model running in full coupled configuration'
      endif
      if ( atmosonly.eq.0 ) then
#ifdef cyclic_ocean
        print *,' Model running in cyclic ocean configuration'
#else
        print *,' Model running in finite box ocean configuration'
#endif
        if ( oceanonly.eq.0 ) then
*         tau-udiff option only applies in coupled mode
*         i.e. neither atmos-only nor ocean_only are set
#ifdef tau_udiff
          print *,' Windstress depends on atmos.-ocean vel. diff.'
#else
          print *,' Windstress depends on atmos. velocity only'
#endif
        endif
      endif
#ifdef get_areav
      print *,' Model computes area average statistics'
#else
      print *,' Model does not compute area average stats'
#endif
      print *,' Model does not accumulate covariance stats'
      if ( hflxsb.eq.1 .and. hflxnb.eq.1 ) then
        print *,' '
        print *,' Invalid model configuration: sb_hflux and'
        print *,' nb_hflux options cannot both be selected'
        print *,' Program terminates'
        stop
       else if ( hflxsb.eq.1 ) then
        print *,' Model running with modified o.m.l. southern b.c.'
       else if ( hflxnb.eq.1 ) then
        print *,' Model running with modified o.m.l. northern b.c.'
       else
        print *,' Model running with no heat flux N & S'
      endif
*     Check some of the grid parameters
      if ( nla.lt.2 .or. nlo.lt.2 ) then
        print *,' '
        print *,' nla, nlo = ',nla,nlo
        print *,' Inadequate nla or nlo, needs to be at least 2'
        print *,' Program terminates'
        stop
      endif
#ifdef cyclic_ocean
*     Check cyclic ocean case is properly set up
      if ( nxta.ne.nxaooc ) then
        print *,' '
        print *,' nxta, nxaooc = ',nxta,nxaooc
        print *,' For cyclic ocean nxta should equal nxaooc'
        print *,' Program terminates'
        stop
      endif
#else
*     Check finite box ocean case is properly set up
      if ( nxta.lt.nxaooc ) then
        print *,' '
        print *,' nxta, nxaooc = ',nxta,nxaooc
        print *,' Inadequate nxta needs to be at least nxaooc'
        print *,' Program terminates'
        stop
      endif
#endif
      if ( nyta.lt.nyaooc ) then
        print *,' '
        print *,' nyta, nyaooc = ',nyta,nyaooc
        print *,' Inadequate nyta needs to be at least nyaooc'
        print *,' Program terminates'
        stop
      endif
*     Check heat-flux boundary conditions
#ifdef sb_hflux
      if ( fnot.lt.0.0d0 ) then
        print *,' '
        print *,' Southern boundary o.m.l. heat-flux activated'
        print *,' Sign of fnot -> running in southern hemisphere'
        print *,' These are inconsistent choices'
        print *,' Program terminates'
        stop
      endif
#endif
#ifdef nb_hflux
      if ( fnot.gt.0.0d0 ) then
        print *,' '
        print *,' Northern boundary o.m.l. heat-flux activated'
        print *,' Sign of fnot -> running in northern hemisphere'
        print *,' These are inconsistent choices'
        print *,' Program terminates'
        stop
      endif
#endif

*     Derive atmos gridspace and ocean timestep
*     -----------------------------------------
*     Derive larger from smaller to ensure integer ratio
      dxa = ndxr*dxo
      dto = nstr*dta
      write(*,201) '  Atmos/ocean grid ratio ndxr = ',ndxr
      write(*,201) '  Oc/atm. timestep ratio nstr = ',nstr
      write(*,201) '  Atmos. gridcells over ocean = ',nxaooc,nyaooc
      write(*,201) '  Ocn start indices  nx1, ny1 = ',nx1,ny1
      write(*,214) '  Coriolis par. f0 (rad s^-1) = ',fnot
      write(*,214) '  Beta =df/dy (rad s^-1 m^-1) = ',beta

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
      write(*,205) '  Rossby number   Beta*ylo/f0 = ',beta*ylo/abs(fnot)
      write(*,214) '  f range S -> N   (rad s^-1) = ',
     &             fnot+beta*yporel(1),fnot+beta*yporel(nypo)
      write(*,214) '  Midlatitude Coriolis param  = ',
     &             fnot+beta*0.5d0*( yporel(1) + yporel(nypo) )
      write(*,205) '  Timestep dto      (minutes) = ',dto/60.0d0
      write(*,201) '  No. of timesteps per day    = ',nint(secday/dto)
      write(*,203) '  Mixed layer thickness   (m) = ',hmoc
      write(*,213) '  Sp. ht. cap. (J kg^-1 K^-1) = ',cpoc
      write(*,213) '  SST grad-2 diff  (m^2 s^-1) = ',st2d
      call diffts (2, nlo, st2d, 1, dxo, rdefoc)
      write(*,213) '  SST grad-4 diff  (m^4 s^-1) = ',st4d
      call diffts (4, nlo, st4d, 1, dxo, rdefoc)
      write(*,205) '  Abs. pot. temp. tabsoc  (K) = ',
     &             (tabsoc(k),k=1,nlo)
      write(*,203) '  Layer thicknesses hoc   (m) = ',(hoc(k),k=1,nlo)
      write(*,203) '  Total thickness   hto   (m) = ',hto
      write(*,207) '  Reduced gravities  (m s^-2) = ',
     &             (gpoc(k),k=1,nlo-1)
      write(*,206) '  Baroclinic wavespeeds (m/s) = ',
     &             (cphsoc(k),k=2,nlo)
      write(*,206) '  Courant number(s)           = ',
     &             ( (dto/dxo)*cphsoc(k),k=2,nlo)
      write(*,204) '  Deformation radii      (km) = ',
     &             (1.0d-3*rdefoc(k),k=2,nlo)
      write(*,205) '  Gridlengths per radius      = ',
     &             (rdefoc(k)/dxo,k=2,nlo)
      write(*,213) '  Del-sqd coeffts  (m^2 s^-1) = ',(ah2oc(k),k=1,nlo)
      call diffts (2, nlo, ah2oc, nlo, dxo, rdefoc)
      write(*,213) '  Del-4th coeffts  (m^4 s^-1) = ',(ah4oc(k),k=1,nlo)
      call diffts (4, nlo, ah4oc, nlo, dxo, rdefoc)
      write(*,204) '  Munk b.l. width scale  (km) = ',
     &             (1.0d-3*(ah4oc(k)/beta)**0.2d0,k=1,nlo)
      write(*,203) '  Bottom Ekm. layer thickness = ',delek
      write(*,213) '  Bottom layer Ekman number   = ',
     &             (delek/hoc(nlo))**2
      write(*,203) '  Spindown timescale   (days) = ',
     &             2.0d0*hoc(nlo)/(abs(fnot)*delek)/secday
      write(*,213) '  Mixed BC coeff. bccooc (nd) = ',bccooc

      call eigmod (nla, gpat, hat, 'Atmosphere', amatat,
     &             cphsat, rdefat, rdm2at, ctl2mat, ctm2lat)
      print *,' '
      print *,' Atmospheric parameters:'
      print *,' -----------------------'
      write(*,201) '  No. of atmos. QG layers nla = ',nla
      write(*,201) '  No. of gridcells nxta, nyta = ',nxta,nyta
      write(*,201) '  At ocean res., no. of cells = ',nxtaor,nytaor
      write(*,204) '  Gridlength dxa         (km) = ',1.0d-3*dxa
      write(*,203) '  Domain sizes xla, yla  (km) = ',
     &             1.0d-3*xla,1.0d-3*yla
      write(*,205) '  Rossby number   Beta*yla/f0 = ',beta*yla/abs(fnot)
      write(*,214) '  f range S -> N   (rad s^-1) = ',
     &             fnot+beta*yparel(1),fnot+beta*yparel(nypa)
      write(*,214) '  Midlatitude Coriolis param  = ',
     &             fnot+beta*0.5d0*( yparel(1) + yparel(nypa) )
      write(*,205) '  Timestep dta      (minutes) = ',dta/60.0d0
      write(*,201) '  No. of timesteps per day    = ',nint(secday/dta)
      write(*,203) '  Mixed layer thickness   (m) = ',hmat
      write(*,203) '  Min. mixed layer thick. (m) = ',hmamin
      write(*,213) '  Quad. drag coefft Cd (nond) = ',cdat
      write(*,213) '  Sp. ht. cap. (J kg^-1 K^-1) = ',cpat
      write(*,213) '  AST grad-2 diff  (m^2 s^-1) = ',at2d
      call diffts (2, nla, at2d, 1, dxa, rdefat)
      write(*,213) '  AST grad-4 diff  (m^4 s^-1) = ',at4d
      call diffts (4, nla, at4d, 1, dxa, rdefat)
      write(*,213) '  hmix diffusivity (m^2 s^-1) = ',ahmd
      call diffts (2, nla, ahmd, 1, dxa, rdefat)
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
      write(*,204) '  Deformation radii      (km) = ',
     &             (1.0d-3*rdefat(k),k=2,nla)
      write(*,205) '  Gridlengths per radius      = ',
     &             (rdefat(k)/dxa,k=2,nla)
      write(*,213) '  Del-4th coeffts  (m^4 s^-1) = ',(ah4at(k),k=1,nla)
      call diffts (4, nla, ah4at, nla, dxa, rdefat)
      write(*,213) '  Mixed BC coeff. bccoat (nd) = ',bccoat

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
  225 format(a,i2,a,9f13.5)
  226 format(a,i2,a,9f13.6)

*     Topography section
*     ==================
      call topset (topocname, topatname, outdir, lenod, .true.)

*     Radiation section
*     =================
*     Compute mean state radiative balance and perturbation
*     radiation coefficients A, B, C and D. Also compute atmosphere
*     and ocean mixed layer temperatures that ensure equilibrium
      call radiat

*     Initialise pressure and temperature fields
*     ==========================================
      if ( name.eq.'zero' ) then
        call zeroin
        tini = 0.0d0
       else if ( name.eq.'rbal' ) then
        call rbalin
        tini = 0.0d0
       else
#ifdef use_netcdf
!! Read netCDF restart dump
        status = nf_open(name,NF_NOWRITE,restid)
        if ( status.ne.NF_NOERR ) call handle_err (status, subnam)
        call restart_nc
        status = nf_close(restid)
        if (status.ne.NF_NOERR) call handle_err (status, subnam)
#else
*       Read unformatted restart dump
        open (11, file=name, form='unformatted', status='old')

c   time stamp
        read(11) tini

c   ocean pressure first
#  ifndef atmos_only
        read(11) po,pom
#  endif

c   initial condition for pa
#  ifndef ocean_only
        read(11) pa,pam
#  endif

c   sst initial condition
        read(11) sst,sstm

c   atmospheric surface temperature
        read(11) ast,astm

c   atmospheric mixed layer thickness
        read(11) hmixa,hmixam

        close (11)
#endif
      endif

*     Assume solution initially valid
      solnok = .true.

*     Derive run length and interval parameters
*     -----------------------------------------
      print *,' '
      print *,' Run length and interval parameters:'
      print *,' -----------------------------------'
      tend = tini + trun
      write(*,204) '  Start time tini     (years) = ',tini
      write(*,204) '  Run length trun     (years) = ',trun
      write(*,204) '  Final time tend     (years) = ',tend
      nsteps0 = nint(tini*secsyr/dta)
      nsteps = nint(tend*secsyr/dta)
      ntsrun = nsteps - nsteps0
      write(*,201) '  Start no. of (atmos)  steps = ',nsteps0
      write(*,201) '  Final no. of (atmos)  steps = ',nsteps
      write(*,201) '  Total no. of (atmos)  steps = ',ntsrun
*     Do all dumps when atmos. and ocean are in step
      nvalid = nint( valday*secday/dto )*nstr
      noutoc = nint( odiday*secday/dto )*nstr
      noutat = nint( adiday*secday/dto )*nstr
      write(*,203) '  Validity test int.   (days) = ',valday
      write(*,203) '  Ocean data dump int. (days) = ',odiday
      write(*,203) '  Atmos data dump int. (days) = ',adiday
      write(*,203) '  Diagnostics dump in. (days) = ',dgnday
      resdmp = resday.gt.0.0d0
      if ( resdmp ) then
        write(*,204) '  Restart dump interval (day) = ',resday
        noutre = nint( resday*secday/dto )*nstr
        iorest = 12
      endif
      nocmon = nint( dgnday*secday/dto )*nstr

*     Compute initial potential vorticities from pressures
*     ====================================================
*     Vorticity computations removed as not needed for testing xforc

!! Read in mean ocean temperature if we're in atmosphere only mode
#ifdef atmos_only
#  ifdef use_netcdf
      print *,' Mean ocean state for atmos_only case read from netCDF'
      status = nf_open('avges.nc',NF_NOWRITE,tempid)
      if ( status.ne.NF_NOERR ) call handle_err (status, subnam)
      status = nf_inq_varid(tempid, 'sst', varid)
      if (status.ne.NF_NOERR) call handle_err (status, subnam)
      status = nf_get_var_double(tempid, varid, sst)
      if (status.ne.NF_NOERR) call handle_err (status, subnam)
      status = nf_get_var_double(tempid, varid, sstm)
      if (status.ne.NF_NOERR) call handle_err (status, subnam)
      status = nf_close(tempid)
      if (status.ne.NF_NOERR) call handle_err (status, subnam)
#  else
      print *,' Mean forcing for atmos_only case read from binary'
      inquire (file=atforav, exist=fexist)
      if ( fexist ) then
        open (50, file=atforav, form='unformatted', status='old')
        read (50) sst
        do j=1,nyto
          do i=1,nxto
            sstm(i,j) = sst(i,j)
          enddo
        enddo
        close(50)
        print *,' Atmos mean forcing filename = ',atforav
       else
        print *,' Mean atmos forcing unformatted binary file not found'
        print *,' Filename = ',atforav
        print *,' Program terminates'
        stop
      endif
#  endif
#endif

!! Read in mean forcings if we're in ocean only mode
*     Need to supply fnetoc, tauxo and tauyo
#ifdef ocean_only
#  ifdef use_netcdf
      print *,' Mean forcing for ocean_only case read from netCDF'
      status = nf_open('avges.nc',NF_NOWRITE,tempid)
      if ( status.ne.NF_NOERR ) call handle_err (status, subnam)
      status = nf_inq_varid(tempid, 'fnetoc', varid)
      if (status.ne.NF_NOERR) call handle_err (status, subnam)
      status = nf_get_var_double(tempid, varid, fnetoc)
      if (status.ne.NF_NOERR) call handle_err (status, subnam)
      status = nf_inq_varid(tempid, 'tauxo', varid)
      if (status.ne.NF_NOERR) call handle_err (status, subnam)
      status = nf_get_var_double(tempid, varid, tauxo)
      if (status.ne.NF_NOERR) call handle_err (status, subnam)
      status = nf_inq_varid(tempid, 'tauyo', varid)
      if (status.ne.NF_NOERR) call handle_err (status, subnam)
      status = nf_get_var_double(tempid, varid, tauyo)
      if (status.ne.NF_NOERR) call handle_err (status, subnam)
      status = nf_close(tempid)
      if (status.ne.NF_NOERR) call handle_err (status, subnam)
#  else
      print *,' Mean forcing for ocean_only case read from binary'
      inquire (file=ocforav, exist=fexist)
      if ( fexist ) then
        open (50, file=ocforav, form='unformatted', status='old')
        read (50) fnetoc
        read (50) tauxo
        read (50) tauyo
        close(50)
        print *,' Ocean mean forcing filename = ',ocforav
       else
        print *,' Mean ocean forcing unformatted binary file not found'
        print *,' Filename = ',ocforav
        print *,' Program terminates'
        stop
      endif
#  endif
#endif

!! Next section is for the xcexp experiments
!! There are a few options which require adjustments to
!! the comments to make them current

!! Option 1: use astbar as determined by radiative equilibrium:
      do j=1,nyta
        do i=1,nxta
           xc1ast(i,j) = ( 1.0d0 - xcexp )*astbar(j)
        enddo
      enddo

*     Initialise forcing arrays
*     =========================
      call xforc

      nt = nsteps0
      tday = tini*daysyr
      tyrs = tini
      call prsamp (nt)

      stop
      end

c***********************************************************************
c
      SUBROUTINE ipbget (buffer, iounit)
*
*     Reads records from unit iounit until a valid one (i.e.
*     one not marked with a "comment" character) is found,
*     then returns this valid character string for processing.
*     The comment marker in Q-GCM is deemed to be
*     an exclamation mark "!" in the first column.

*     Modules
      USE parameters
*
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
      SUBROUTINE zeroin
*
*     Set initial state to zero pressure, and radiative
*     equilibrium with unperturbed (background) forcing,
*     i.e. mixed layer temperature anomalies are all zero.

*     Modules
      USE parameters
#ifndef ocean_only
      USE atstate
#endif
#ifndef atmos_only
      USE ocstate
#endif
      USE intrfac, ONLY : sst, sstm, ast, astm, hmixa, hmixam, hmat

      IMPLICIT NONE

*     Subroutine arguments
*
*     Local variables
      integer i,j,k

*     Initialise atmosphere fields (incl. mixed layer)
*     ------------------------------------------------
!$OMP PARALLEL DEFAULT (NONE)
!$OMP&         PRIVATE (i,j,k)
!$OMP&         SHARED  (ast,astm,hmixa,hmixam,hmat)
#ifndef ocean_only
!$OMP&         SHARED  (pa,pam)

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
#endif

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

*     Initialise ocean fields (incl. mixed layer)
*     -------------------------------------------
!$OMP PARALLEL DEFAULT (NONE)
!$OMP&         PRIVATE (i,j,k)
!$OMP&         SHARED  (sst,sstm)
#ifndef atmos_only
!$OMP&         SHARED  (po,pom)

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
#endif

*     Initialise oceanic mixed layer rel. temp.
!$OMP DO SCHEDULE (STATIC)
      do j=1,nyto
        do i=1,nxto
          sst(i,j) = 0.0d0
          sstm(i,j) = 0.0d0
        enddo
      enddo
!$OMP END DO NOWAIT

!$OMP END PARALLEL

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
#ifndef ocean_only
      USE atstate
#endif
#ifndef atmos_only
      USE ocstate
#endif
      USE intrfac, ONLY : sst, sstm, sstbar, ast, astm, astbar,
     &                    hmixa, hmixam, hmat
      USE radiate

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

#ifndef ocean_only
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
#endif

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

*     Initialise atmosphere fields (incl. mixed layer)
*     ------------------------------------------------
!$OMP PARALLEL DEFAULT (NONE)
!$OMP&         PRIVATE (i,j,k,play)
!$OMP&         SHARED  (ast,astm,astbar,hmixa,hmixam,hmat)
#ifndef ocean_only
!$OMP&         SHARED  (pa,pam,plfac,yparel)

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
#endif

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

*     Initialise ocean fields (incl. mixed layer)
*     -------------------------------------------
!$OMP PARALLEL DEFAULT (NONE)
!$OMP&         PRIVATE (i,j,k,play)
!$OMP&         SHARED  (sst,sstm,sstbar)
#ifndef atmos_only
!$OMP&         SHARED  (po,pom)

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
#endif

*     Initialise oceanic mixed layer rel. temp.
!$OMP DO SCHEDULE (STATIC)
      do j=1,nyto
        do i=1,nxto
          sst(i,j) = sstbar(j)
          sstm(i,j) = sstbar(j)
        enddo
      enddo
!$OMP END DO NOWAIT

!$OMP END PARALLEL

      END SUBROUTINE rbalin
c
c***********************************************************************
c
      SUBROUTINE prsamp (nt)
*
*     Print out a sample of the current state. Useful for
*     tracking progress, or for testing/debugging purposes

*     Modules
      USE parameters
#ifndef ocean_only
      USE atconst
      USE atstate
#endif
#ifndef atmos_only
      USE occonst
      USE ocstate
#endif
      USE intrfac, ONLY : sst, ast, hmixa
      USE timinfo, ONLY : tday, tyrs
      USE monitor, ONLY : pavgoc, qavgoc

      IMPLICIT NONE

*     Subroutine arguments
      integer, INTENT(IN) :: nt
*
*     Local parameters
      double precision extrem
      parameter ( extrem=1.0d30 )
*
*     Local variables
      integer nxco,nyco,nxca,nyca,i,j,k
      double precision astmin,astmax,hmxmin,hmxmax,sstmin,sstmax,
     &                 poref(nlo),paref(nla)

*     Centres of domains for sample printing
      nxco = (nxpo+1)/2
      nyco = (nypo+1)/2
      nxca = (nxpa+1)/2
      nyca = (nypa+1)/2

*     Define reference pressures on equatorward side of domain
*     --------------------------------------------------------
#ifndef ocean_only
      if ( fnot.gt.0.0d0 ) then
        do k=1,nla
          paref(k) = pa(1,1,k)
        enddo
       else if ( fnot.lt.0.0d0 ) then
        do k=1,nla
          paref(k) = pa(1,nypa,k)
        enddo
      endif
#endif
#ifndef atmos_only
      if ( fnot.gt.0.0d0 ) then
        do k=1,nlo
          poref(k) = po(1,1,k)
        enddo
       else if ( fnot.lt.0.0d0 ) then
        do k=1,nlo
          poref(k) = po(1,nypo,k)
        enddo
      endif
#endif

      print *,' '
      write(*,'(a,i12,f11.2,f11.4)')
     &      '  Sample output at nt, tday, tyrs = ',nt,tday,tyrs

*     Print some ocean spot values and extrema
*     ----------------------------------------
#ifndef atmos_only
      write(*,217) '  po(k) at centre = ',(po(nxco,nyco,k),k=1,nlo)
**    write(*,217) '  qo(k) at centre = ',(qo(nxco,nyco,k),k=1,nlo)
**    write(*,217) '  po(1) on S. bdy = ',(po(i,  1 ,1),i=nxco-1,nxco+1)
**    write(*,217) '  qo(1) on S. bdy = ',(qo(i,  1 ,1),i=nxco-1,nxco+1)
**    write(*,217) '  po(1) on N. bdy = ',(po(i,nypo,1),i=nxco-1,nxco+1)
**    write(*,217) '  qo(1) on N. bdy = ',(qo(i,nypo,1),i=nxco-1,nxco+1)
**    write(*,217) '   sst  on S. bdy = ',(sst(i,  1 ),i=nxco-1,nxco+1)
**    write(*,217) '   sst  on N. bdy = ',(sst(i,nyto),i=nxco-1,nxco+1)
      write(*,217) '  po(k)  average  = ',(pavgoc(k),k=1,nlo)
**    write(*,217) '  qo(k)  average  = ',(qavgoc(k),k=1,nlo)
**    write(*,217) '  Ocean kealoc(k) = ',(kealoc(k),k=1,nlo)
**    write(*,217) '  Ocean et2moc(k) = ',(et2moc(k),k=1,nlo-1)
**    write(*,217) '  osfmin(k)  (Sv) = ',(osfmin(k),k=1,nlo)
**    write(*,217) '  osfmax(k)  (Sv) = ',(osfmax(k),k=1,nlo)
**    write(*,217) '  occirc(k)  (Sv) = ',(occirc(k),k=1,nlo)
*     Meridional section at centre
**    write(*,*) ' Meridional section at ocean centre'
**    write(*,'(a,a)') '     j         po(1)             qo(1)       ',
**   &           '      wekpo           porel(1)           ddynoc'
**    do j=nypo,1,-1
**      write(*,240) j,po(nxco,j,1),qo(nxco,j,1),wekpo(nxco,j),
**   &               po(nxco,j,1)-poref(1),ddynoc(nxco,j)
**    enddo
*     Southern boundary values
**    write(*,*) ' Southern boundary values for ocean'
**    write(*,*) '    i        po(1)             qo(1)             ',
**   &           'wekpo           porel(1)'
**    do i=1,nxpo
**      write(*,240) i,po(i,1,1),qo(i,1,1),wekpo(i,1),
**   &               po(i,1,1)-poref(1)
**    enddo
      sstmin =  extrem
      sstmax = -extrem
!$OMP PARALLEL DEFAULT (NONE)
!$OMP&         PRIVATE (i,j)
!$OMP&         SHARED  (sstmin,sstmax,sst)

!$OMP DO SCHEDULE (STATIC)
!$OMP&   REDUCTION(min:sstmin) REDUCTION(max:sstmax)
      do j=1,nyto
        do i=1,nxto
          sstmin = min( sstmin, sst(i,j) )
          sstmax = max( sstmax, sst(i,j) )
        enddo
      enddo
!$OMP END DO NOWAIT

!$OMP END PARALLEL
      write(*,217) '  s.s.t: min, max = ',sstmin,sstmax
#endif

*     Print some atmos. spot values and extrema
*     -----------------------------------------
#ifndef ocean_only
      write(*,217) '  pa(k) at centre = ',(pa(nxca,nyca,k),k=1,nla)
**    write(*,217) '  qa(k) at centre = ',(qa(nxca,nyca,k),k=1,nla)
**    write(*,217) '  pa(1) on S. bdy = ',(pa(i,  1 ,1),i=nxca-1,nxca+1)
**    write(*,217) '  qa(1) on S. bdy = ',(qa(i,  1 ,1),i=nxca-1,nxca+1)
**    write(*,217) '  pa(1) on N. bdy = ',(pa(i,nypa,1),i=nxca-1,nxca+1)
**    write(*,217) '  qa(1) on N. bdy = ',(qa(i,nypa,1),i=nxca-1,nxca+1)
**    write(*,217) '   ast  on S. bdy = ',(ast(i,  1 ),i=nxca-1,nxca+1)
**    write(*,217) '   ast  on N. bdy = ',(ast(i,nyta),i=nxca-1,nxca+1)
**    write(*,217) '  Atmos kealat(k) = ',(kealat(k),k=1,nla)
**    write(*,217) '  Atmos et2mat(k) = ',(et2mat(k),k=1,nla-1)
*     Meridional section at centre
**    write(*,*) ' Meridional section at atmos. centre'
**    write(*,'(a,a)') '     j         pa(1)             qa(1)       ',
**   &           '      wekpa           parel(1)           ddynat'
**    do j=nypa,1,-1
**      write(*,240) j,pa(nxca,j,1),qa(nxca,j,1),wekpa(nxca,j),
**   &               pa(nxca,j,1)-paref(1),ddynat(nxca,j)
**    enddo
      astmin =  extrem
      astmax = -extrem
      hmxmin =  extrem
      hmxmax = -extrem
!$OMP PARALLEL DEFAULT (NONE)
!$OMP&         PRIVATE (i,j)
!$OMP&         SHARED  (astmin,astmax,ast,hmxmin,hmxmax,hmixa)

!$OMP DO SCHEDULE (STATIC)
!$OMP&   REDUCTION(min:astmin) REDUCTION(max:astmax)
!$OMP&   REDUCTION(min:hmxmin) REDUCTION(max:hmxmax)
      do j=1,nyta
        do i=1,nxta
          astmin = min( astmin, ast(i,j) )
          astmax = max( astmax, ast(i,j) )
          hmxmin = min( hmxmin, hmixa(i,j) )
          hmxmax = max( hmxmax, hmixa(i,j) )
        enddo
      enddo
!$OMP END DO NOWAIT

!$OMP END PARALLEL
      write(*,217) '  a.s.t: min, max = ',astmin,astmax
      write(*,217) '  hmixa: min, max = ',hmxmin,hmxmax
#endif

  217 format(a,1p,9d15.7)
  240 format(i6,1p,5d18.10)

      END SUBROUTINE prsamp
c
c***********************************************************************
c
      SUBROUTINE diffts (nord, nl, coeff, ncoef, dx, rdef)

*     Computes diffusive decay timescales for circular eddies whose radii
*     are the baroclinic Rossby radii, and for two-gridpoint noise.
*     See section 8.6 of the Userguide for derivation of timescales.

*     Input arguments:
*     nord  : order of the diffusive term
*     nl    : no. of QG layers  (=> nl-1 baroclinic modes)
*     coeff : vector of diffusion coefficients (should be .GE. 0)
*     ncoef : length of coefficient vector
*     dx    : gridlength (m)
*     rdef  : vector of nl modal deformation radii (m)
*             (infinite value for barotropic mode replaced by 0.0)
*     (all the above are unchanged on exit)

*     Modules
      USE parameters

      IMPLICIT NONE

*     Subroutine arguments
      integer nord,nl,ncoef
      double precision coeff(ncoef),dx,rdef(nl)
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
      double precision tdamp(nlmax),sinfac

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
*     at the Rossby radius for each baroclinic mode
      do m=2,nl
        sinfac = 2.0d0*sin( PIBY2*dx/rdef(m) )/dx
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

*     Compute decay timescale for two-gridpoint noise
*     for each coefft, avoiding infinities if coefft = 0
      do k=1,ncoef
        if ( coeff(k).eq.0.0d0 ) then
          tdamp(k) = 0.0d0
         else
          tdamp(k) = (0.5d0*dx)**nord/coeff(k)/3600.0d0
        endif
      enddo
      write(*,205) '  Gridpoint timescale (hours) = ',
     &             (tdamp(k),k=1,ncoef)

  205 format(a,9f13.5)
  225 format(a,i2,a,9f13.5)

      END SUBROUTINE diffts
c
c***********************************************************************
