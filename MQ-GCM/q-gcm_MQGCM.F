c***********************************************************************
c     Q-GCM Version 1.5.0 : last modified 24/09/2013
c
c     Modified on 11/13/2017 by S. Kravtsov to include high-res. atmospheric
c     output over ocean
c
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

*     Modules
      USE parameters
      USE atconst
      USE occonst
      USE radiate
      USE timinfo
      USE monitor
#ifndef ocean_only
      USE athomog
      USE amlsubs
      USE qgasubs
      USE atisubs
#endif
#if defined tau_udiff || !defined ocean_only
      USE atstate
#endif
#ifndef atmos_only
      USE ochomog
      USE ocstate
      USE omlsubs
      USE qgosubs
      USE ocisubs
#endif
      USE intrfac
      USE conhoms
      USE vorsubs
      USE xfosubs, ONLY : xforc
      USE eigmode
      USE valsubs
      USE radsubs
      USE topsubs
#ifdef use_netcdf
      USE nc_subs
#  ifndef ocean_only
#    ifdef highres_output
      USE hra_out
#    endif /* highres_output */
#  endif /* not ocean_only */
#endif /* use_netcdf */
*     Diagnostics modules
      USE mondiag
#ifdef get_areav
      USE areasubs
#endif
#ifndef atmos_only
#  ifdef qoc_diag
      USE qocdiag
#  endif /* qoc_diag */
#endif /* atmos_only */
#ifdef get_covar
      USE covaria
#endif /* get_covar */
      USE timavge

c      USE xfosubs, ONLY : fsprim
c      USE occonst, ONLY : ytorel

      IMPLICIT NONE
      
c      double precision fsp

#ifdef use_netcdf
      INCLUDE 'netcdf.inc'
#endif /* use_netcdf */
*
*     Local parameters
      character (len=*), parameter :: subnam = 'q-gcm_main'
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
      double precision st2da(1),st4da(1),at2da(1),at4da(1),ahmda(1)
*
*     Arrays passed as arguments
*
*     Local variables
      integer i,j,k,i1,iorest,nsko,nska,nt,lenod,lename,oceanonly,
     &        atmosonly,getcovar,outfloc(7),outflat(9),
     &        cyclicoc,hflxsb,hflxnb,tauudiff,
     &         tempfdbck,partialco,moist1,drylatent
      double precision tspind
      character (len=72) :: name
      character (len=80) :: outdir, inpbuf
      logical resdmp,solnok,avgatm,avgocn
*     Extra variables for topography
      character (len=80) :: topocname,topatname

#ifdef use_netcdf
*     Extra variables for netCDF input/output
      integer status,tempid,varid
#  ifndef atmos_only
      double precision poc4dsz
#  endif /* not atmos_only */
#  ifndef ocean_only
      double precision pat4dsz
#  endif /* not ocean_only */
#else
#  ifdef atmos_only
      character (len=*), parameter :: atforav='./atforce.avg.binary'
#  endif /* atmos_only */
#  ifdef ocean_only
      character (len=*), parameter :: ocforav='./ocforce.avg.binary'
#  endif /* ocean_only */
      logical fexist
#endif

*     Extra variables for OpenMP
!$    integer nprocs,OMP_GET_NUM_PROCS,nthmax,OMP_GET_MAX_THREADS,
!$   &        numthr,OMP_GET_NUM_THREADS,thrnum,OMP_GET_THREAD_NUM
!$    logical dynadj,OMP_GET_DYNAMIC,nested,OMP_GET_NESTED

# ifdef moist
      double precision dt1,dt2,dtainv,xfa0m
# endif

      print *,' '
      write(*,*) ' Computational parameters:'
      write(*,*) ' ========================='
      write(*,*) ' Model version is Q-GCM v1.5.0'
*     Compute typical array sizes, and local arrays for stack allocation
      call memreq

*     Examine OpenMP environment
*     --------------------------
!$    nprocs = OMP_GET_NUM_PROCS()
!$    nthmax = OMP_GET_MAX_THREADS()
!$    dynadj = OMP_GET_DYNAMIC()
!$    nested = OMP_GET_NESTED()
!$    write(*,*) ' '
!$    write(*,*) ' OpenMP parallelism is activated'
!$    write(*,*) ' -------------------------------'
!$    write(*,'(a,i5)') '  No. of processors available = ',nprocs
!$    write(*,'(a,i3)') '  Max. no. of threads available = ',nthmax
!$    write(*,*) ' Dynamic adjustment = ',dynadj
!$    write(*,*) ' Nested parallelism = ',nested
*     Test OpenMP is functioning correctly
!$    write(*,*) ' Outside parallel section:'
!$    numthr = OMP_GET_NUM_THREADS()
!$    write(*,*) ' Number of threads = ',numthr
!$    thrnum = OMP_GET_THREAD_NUM()
!$    write(*,'(a,i4)') '  thrnum = ',thrnum
!$    write(*,*) ' Test trivial parallel loop:'
!$OMP PARALLEL DEFAULT (NONE)
!$OMP&         PRIVATE (j,thrnum)
!$OMP&         SHARED (numthr)
!$OMP SINGLE
!$    numthr = OMP_GET_NUM_THREADS()
!$    write(*,*) ' Number of threads = ',numthr
!$OMP END SINGLE
!$OMP DO SCHEDULE (STATIC)
!$    do j=1,numthr
!$      thrnum = OMP_GET_THREAD_NUM()
!$      write(*,'(a,i6,i4)') '  j, thrnum = ',j,thrnum
!$    enddo
!$OMP END DO NOWAIT
!$OMP END PARALLEL
*     The following loop is the most realistic test of OpenMP //ism
*     Enable (if necessary) by removing the first ! from each line
*     Comment out again once OpenMP is working correctly,
*     to avoid lots of spurious output
!!$    write(*,*) ' Test typical j-range loop (ocean p-points)'
!!$OMP PARALLEL DO DEFAULT (NONE)
!!$OMP&         PRIVATE (j,thrnum)
!!$OMP&         SCHEDULE (STATIC)
!!$    do j=1,nypo
!!$      thrnum = OMP_GET_THREAD_NUM()
!!$      write(*,'(a,i6,i4)') '  j, thrnum = ',j,thrnum
!!$    enddo
!!$OMP END PARALLEL DO

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
#endif /* atmos_only */
#ifdef get_covar
      getcovar = 1
#else
      getcovar = 0
#endif /* get_covar */
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
#ifdef temp_fdbck
      tempfdbck = 1
#else
      tempfdbck = 0
#endif
#ifdef partial_coup
      partialco = 1
#else
      partialco = 0
#endif
#ifdef moist
      moist1 = 1
#else
      moist1 = 0
#endif
#ifdef dry_latent
      drylatent=1
#else
      drylatent=0
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
#ifdef qoc_diag
      print *,' Model computes ocean dq/dt diagnostics '
#else
      print *,' Model does not compute ocean dq/dt diagnostics'
#endif
      if ( getcovar.eq.1 ) then
        print *,' Model accumulates covariance statistics'
        write(*,201) '  Sampling intervals nscvat, nscvoc = ',
     &               nscvat,nscvoc
       else
        print *,' Model does not accumulate covariance stats'
      endif
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
#ifdef get_covar
      if ( mod(nxta,nscvat).ne.0 .or. mod(nyta,nscvat).ne.0 ) then
        print *,' '
        print *,' nxta, nyta, nscvat = ',nxta,nyta,nscvat
        print *,' Need nxta, nyta both to be integer multiples of'
        print *,' nscvat for covariance subsampling to work correctly'
        print *,' Program terminates'
        stop
      endif
      if ( mod(nxto,nscvoc).ne.0 .or. mod(nyto,nscvoc).ne.0 ) then
        print *,' '
        print *,' nxto, nyto, nscvoc = ',nxto,nyto,nscvoc
        print *,' Need nxto, nyto both to be integer multiples of'
        print *,' nscvoc for covariance subsampling to work correctly'
        print *,' Program terminates'
        stop
      endif
#endif /* get_covar */
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
*
*     Topography section
*     ==================
      call topset (topocname, topatname, outdir, lenod, .true.)

*     Radiation section
*     =================
*     Compute mean state radiative balance and perturbation
*     radiation coefficients A, B, C and D. Also compute atmosphere
*     and ocean mixed layer temperatures that ensure equilibrium
      call radiat


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
      write(*,213) '  Sp. ht. cap. (J kg^-1 K^-1) = ',cpoc
      write(*,213) '  SST grad-2 diff  (m^2 s^-1) = ',st2d
      st2da(1)=st2d
      call diffts (2, nlo, st2da, 1, dxo, rdefoc, dto)
      write(*,213) '  SST grad-4 diff  (m^4 s^-1) = ',st4d
      st4da(1)=st4d
      call diffts (4, nlo, st4da, 1, dxo, rdefoc, dto)
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
      write(*,206) '  Rossby number   Beta*yla/f0 = ',beta*yla/abs(fnot)
      write(*,215) '  f range S -> N   (rad s^-1) = ',
     &             fnot+beta*yparel(1),fnot+beta*yparel(nypa)
      write(*,215) '  Midlatitude Coriolis param  = ',
     &             fnot+beta*0.5d0*( yparel(1) + yparel(nypa) )
      write(*,205) '  Timestep dta      (minutes) = ',dta/60.0d0
      write(*,201) '  No. of timesteps per day    = ',nint(secday/dta)
      write(*,203) '  Mixed layer thickness   (m) = ',hmat
      write(*,203) '  Min. mixed layer thick. (m) = ',hmamin
      write(*,213) '  Quad. drag coefft Cd (nond) = ',cdat
      write(*,213) '  Sp. ht. cap. (J kg^-1 K^-1) = ',cpat
      write(*,213) '  AST grad-2 diff  (m^2 s^-1) = ',at2d
      at2da(1)=at2d
      call diffts (2, nla, at2da, 1, dxa, rdefat, dta)
      write(*,213) '  AST grad-4 diff  (m^4 s^-1) = ',at4d
      at4da(1)=at4d
      call diffts (4, nla, at4da, 1, dxa, rdefat, dta)
      write(*,213) '  hmix diffusivity (m^2 s^-1) = ',ahmd
      ahmda(1)=ahmd
      call diffts (2, nla, ahmda, 1, dxa, rdefat, dta)
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
* 214 format(a,1p,9d13.4)
  215 format(a,1p,9d13.5)
* 225 format(a,i2,a,9f13.5)
* 226 format(a,i2,a,9f13.6)


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
        status = nf_open(name, NF_NOWRITE, restid)
        if ( status.ne.NF_NOERR ) call handle_err (status, subnam)
        call restart_nc
        status = nf_close(restid)
        if ( status.ne.NF_NOERR ) call handle_err (status, subnam)


c******************************************************************************
c
c FOR THE INITIAL RUN (S. Kravtsov, 10/27/2020)
c
c

# ifndef moist
c
c        do j=1,nyta
c          do i=1,nxta
c            hmixa(i,j) = hmat
c            hmixam(i,j)= hmat
c          enddo
c        enddo
c
# else
c
        dt1=tat(2)-tat(1)
        dt2=tat(3)-tat(2)
        dtainv=1.d0/(dt1+f2*dt2)
        xfa0m=0.d0
        do j=1,nyta
          do i=1,nxta
c            hmixa(i,j) = qatm
c            hmixam(i,j)= qatm
c            hum(i,j,1)= qat(1)
c            humm(i,j,1)=qat(1)
c            hum(i,j,2)=qat(2)
c            humm(i,j,2)=qat(2)
c            hum(i,j,3)=qat(3)
c            humm(i,j,3)=qat(3)
c             xfa0(i,j)=0.d0
             xfa0(i,j)=Lv*(rho(1)*precip_i(i,j,1)
     &                    +rho(2)*precip_i(i,j,2)
     &                    +rho(3)*precip_i(i,j,3))*rrcpat*dtainv
             xfa0m=xfa0m+xfa0(i,j)
          enddo
        enddo
c        print*,xfa0m
c        pause
        do j=1,nyta
          do i=1,nxta
             xfa0(i,j)=xfa0(i,j)-xfa0m*atnorm
          enddo
        enddo
c
c
# endif
c
c
c*******************************************************************************
c
c
c
#else
*       Read unformatted restart dump
        open (11, file=name, form='unformatted', status='old')

c   time stamp
        read(11) tini

c   ocean pressure first
#  ifndef atmos_only
        read(11) po,pom
#  endif /* atmos_only */

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


c       pa(68,20,3)=100.d0
c       pam(68,20,3)=100.d0
c
c      do j=1,nyto
c        fsp = fsprim( ytorel(j) )
c        do i=1,nxto
c          sst(i,j)=-26.d0*(fsp/fspco)
c          sstm(i,j)=sst(i,j)
c        enddo
c        print*,sst(1,j)
c       enddo


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
      write(*,204) '  Validity test int.   (days) = ',valday
      write(*,204) '  Ocean data dump int. (days) = ',odiday
      write(*,204) '  Atmos data dump int. (days) = ',adiday
      write(*,204) '  Diagnostics dump in. (days) = ',dgnday
      write(*,204) '  Printout dump inter. (days) = ',prtday
      resdmp = resday.gt.0.0d0
      if ( resdmp ) then
        write(*,204) '  Restart dump interval (day) = ',resday
        noutre = nint( resday*secday/dto )*nstr
        iorest = 12
      endif
      nocmon = nint( dgnday*secday/dto )*nstr
      nocprt = nint( prtday*secday/dto )*nstr
*     Control intervals for time-averaging of forcings
      avgatm = dtavat.gt.0.0d0 .and. oceanonly.eq.0
      if ( avgatm ) then
        write(*,204) '  Atmos. averaging int. (day) = ',dtavat
        ntavat = nint( dtavat*secday/dto )*nstr
        write(*,201) '  No. of steps in avg. inter. = ',ntavat
        if ( mod(ntsrun,ntavat).ne.0 .or. mod(ntavat,2).ne.0 ) then
          print *,' Unsuitable choice of dtavat; program stops'
          stop
        endif
        nmidat = ntavat/2
      endif
      avgocn = dtavoc.gt.0.0d0 .and. atmosonly.eq.0
      if ( avgocn ) then
        write(*,204) '  Ocean  averaging int. (day) = ',dtavoc
        ntavoc = nint( dtavoc*secday/dto )*nstr
        write(*,201) '  No. of steps in avg. inter. = ',ntavoc
        if ( mod(ntsrun,ntavoc).ne.0 .or. mod(ntavoc,2).ne.0 ) then
          print *,' Unsuitable choice of dtavoc; program stops'
          stop
        endif
        nmidoc = ntavoc/2
      endif
#ifdef get_covar
*     Control intervals for taking covariances
      ntcovat = nint( dtcovat*secday/dto )*nstr
      ntcovoc = nint( dtcovoc*secday/dto )*nstr
#endif /* get_covar */

*     Write Matlab-readable copy of input and derived parameters
*     ----------------------------------------------------------
      open (10, file=outdir(1:lenod)//'/input_parameters.m',
     &      status='unknown')
      INCLUDE './out_param.f'
      close(10)

*     Compute initial mass and momentum integrals for constraints
*     -----------------------------------------------------------
      call constr

*     Compute initial potential vorticities from pressures
*     ====================================================
#ifndef atmos_only
*     Ocean
*     -----
*     Internal points
      call qcomp (qo, po, amatoc, yporel, dxom2,
     &            nxpo, nypo, nlo, ddynoc, nlo)
      call qcomp (qom,pom,amatoc, yporel, dxom2,
     &            nxpo, nypo, nlo, ddynoc, nlo)
*     Zonal boundaries + meridional if not cyclic (mixed condition)
      call ocqbdy (qo, po )
      call ocqbdy (qom,pom)
#  ifdef cyclic_ocean
*     Meridional boundaries (periodic)
      call merqcy (qo, po,  amatoc, yporel, dxom2,
     &             nxpo, nypo, nlo, ddynoc, nlo)
      call merqcy (qom,pom, amatoc, yporel, dxom2,
     &             nxpo, nypo, nlo, ddynoc, nlo)
#  endif
#endif /* atmos_only */
#ifndef ocean_only
*     Atmosphere
*     ----------
*     Internal points
      call qcomp (qa, pa, amatat, yparel, dxam2,
     &            nxpa, nypa, nla, ddynat, 1)
      call qcomp (qam,pam,amatat, yparel, dxam2,
     &            nxpa, nypa, nla, ddynat, 1)
*     Zonal boundaries (mixed condition)
      call atqzbd (qa, pa )
      call atqzbd (qam,pam)
*     Meridional boundaries (periodic)
      call merqcy (qa, pa,  amatat, yparel, dxam2,
     &             nxpa, nypa, nla, ddynat, 1)
      call merqcy (qam,pam, amatat, yparel, dxam2,
     &             nxpa, nypa, nla, ddynat, 1)
#endif

!! Read in mean ocean temperature if we're in atmosphere only mode
#ifdef atmos_only
#  ifdef use_netcdf
      print *,' Mean ocean state for atmos_only case read from netCDF'
      status = nf_open('avges.nc', NF_NOWRITE, tempid)
      if ( status.ne.NF_NOERR ) call handle_err (status, subnam)
      status = nf_inq_varid(tempid, 'sst', varid)
      if ( status.ne.NF_NOERR ) call handle_err (status, subnam)
      status = nf_get_var_double(tempid, varid, sst)
      if ( status.ne.NF_NOERR ) call handle_err (status, subnam)
      status = nf_get_var_double(tempid, varid, sstm)
      if ( status.ne.NF_NOERR ) call handle_err (status, subnam)
      status = nf_close(tempid)
      if ( status.ne.NF_NOERR ) call handle_err (status, subnam)
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
#endif /* atmos_only */

!! Read in mean forcings if we're in ocean only mode
*     Need to supply fnetoc, tauxo and tauyo
#ifdef ocean_only
#  ifdef use_netcdf
      print *,' Mean forcing for ocean_only case read from netCDF'
      status = nf_open('avges.nc', NF_NOWRITE, tempid)
      if ( status.ne.NF_NOERR ) call handle_err (status, subnam)
      status = nf_inq_varid(tempid, 'fnetoc', varid)
      if ( status.ne.NF_NOERR ) call handle_err (status, subnam)
      status = nf_get_var_double(tempid, varid, fnetoc)
      if ( status.ne.NF_NOERR ) call handle_err (status, subnam)
      status = nf_inq_varid(tempid, 'tauxo', varid)
      if ( status.ne.NF_NOERR ) call handle_err (status, subnam)
      status = nf_get_var_double(tempid, varid, tauxo)
      if ( status.ne.NF_NOERR ) call handle_err (status, subnam)
      status = nf_inq_varid(tempid, 'tauyo', varid)
      if ( status.ne.NF_NOERR ) call handle_err (status, subnam)
      status = nf_get_var_double(tempid, varid, tauyo)
      if ( status.ne.NF_NOERR ) call handle_err (status, subnam)
      status = nf_close(tempid)
      if ( status.ne.NF_NOERR ) call handle_err (status, subnam)
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

!! Options 2 & 3: both require reading in a file from the
!! local directory which holds averaged ast information:
!      print *,' Reading in mean AST for xcexp experiment'
!      status = nf_open('avges.nc', NF_NOWRITE, tempid)
!      if ( status.ne.NF_NOERR ) call handle_err (status, subnam)
!      status = nf_inq_varid(tempid, 'ast', varid)
!      if ( status.ne.NF_NOERR ) call handle_err (status, subnam)
!      status = nf_get_var_double(tempid, varid, xc1ast)
!      if ( status.ne.NF_NOERR ) call handle_err (status, subnam)

!! Option 2: Directly use the mean AST fields:
!      do j=1,nyta
!        do i=1,nxta
!           xc1ast(i,j) = ( 1.0d0 - xcexp )*xc1ast(i,j)
!        enddo
!      enddo

!! Option 3: Use the AST from data, but average
!! it in the x-direction, so that there are no 'bumps':
!      do j=1,nyta
!        do i=2,nxta
!           xc1ast(1,j) = xc1ast(1,j) + xc1ast(i,j)
!        enddo
!        xc1ast(1,j) = ( 1.0d0 - xcexp )*xc1ast(1,j)/dble(nxta)
!        do i=2,nxta
!           xc1ast(i,j) = xc1ast(1,j)
!        enddo
!      enddo

*     Initialise forcing arrays
*     =========================
      call xforc

#ifndef atmos_only
*     Ocean forcings
*     --------------
*     Set those quantities which have not already been either
*     set by xforc, or in the ocean only case, read from a file
!$OMP PARALLEL DO DEFAULT (NONE)
!$OMP&         PRIVATE (i,j)
!$OMP&         SHARED  (entoc)
!$OMP&         SCHEDULE (STATIC)
      do j=1,nypo
        do i=1,nxpo
          entoc(i,j) = 0.0d0
        enddo
      enddo
!$OMP END PARALLEL DO
      do k=1,nlo-1
        xon(k) = 0.0d0
      enddo
#  ifdef cyclic_ocean
      do k=1,nlo-1
        enisoc(k) = 0.0d0
        eninoc(k) = 0.0d0
      enddo
*     Also zero ocean boundary Jacobian arrays
*     These will be correctly recomputed in the cyclic
*     case, and are not used in the finite box case
      do k=1,nlo
        ajisoc(k) = 0.0d0
        ajinoc(k) = 0.0d0
      enddo
#  endif
#endif /* atmos_only */

#ifndef ocean_only
*     Atmosphere forcings
*     -------------------
      do j=1,nypa
        do i=1,nxpa
          entat(i,j) = 0.0d0
          entat1(i,j)=0.0d0
        enddo
      enddo
      do k=1,nla-1
        enisat(k) = 0.0d0
        eninat(k) = 0.0d0
        xan(k) = 0.0d0
      enddo
#endif

*     Compute tridiagonal coefficients for Helmholtz equations
*     ========================================================
*     Compute Del-sqd part once and for all; add def. rad. later.
*     In all cases the tridiagonal system is diagonally dominant and
*     thus numerically well-behaved. The Del-sqd operator to be
*     inverted is the finite-difference form in both x- and y-directions.
*     Ordering of wavenumber components in bd2oc,
*     bd2at is appropriate to FFTPACK.

#ifndef atmos_only
*     Oceanic
*     -------
      aoc = 1.0d0/( dyo*dyo )
#  ifdef cyclic_ocean
*     Zonally cyclic ocean
      do i=2,nxto/2
        i1 = 2*i - 1
        bd2oc(i1-1) = -2.0d0*aoc
     &               + 2.0d0*dxom2*( cos( (i-1)*TWOPI/nxto ) - 1.0d0 )
        bd2oc( i1 ) = bd2oc(i1-1)
      enddo
*     When using FFTPACK, special values go to 1 and nxto
      bd2oc(  1 ) = -2.0d0*aoc
      bd2oc(nxto) = -2.0d0*aoc - 4.0d0*dxom2
*     Initialise oceanic FFT coefficient array
      call drffti (nxto, oftwrk)
#  else
*     Finite box ocean
      do i=2,nxto
        bd2oc(i-1) = -2.0d0*aoc
     &              + 2.0d0*dxom2*( cos( (i-1)*PI/nxto ) - 1.0d0 )
      enddo
      bd2oc(nxto) = 0.0d0
*     Initialise oceanic FFT coefficient array
      call dsinti (nxto-1, oftwrk)
#  endif
#endif /* atmos_only */

#ifndef ocean_only
*     Atmospheric
*     -----------
      aat = 1.0d0/( dya*dya )
      do i=2,nxta/2
        i1 = 2*i - 1
        bd2at(i1-1) = -2.0d0*aat
     &               + 2.0d0*dxam2*( cos( (i-1)*TWOPI/nxta ) - 1.0d0 )
        bd2at( i1 ) = bd2at(i1-1)
      enddo
*     When using FFTPACK, special values go to 1 and nxta
      bd2at(  1 ) = -2.0d0*aat
      bd2at(nxta) = -2.0d0*aat - 4.0d0*dxam2
*     Initialise atmospheric FFT coefficient array
      call drffti (nxta, aftwrk)
#endif

*     Compute homogeneous solutions for ocean and atmosphere
      call homsol

      noutstepoc = ntsrun/noutoc + 1
      noutstepat = ntsrun/noutat + 1
      numoutsteps = ntsrun/nocmon + 1
      nt = nsteps0
      ntdone = 0
      tday = tini*daysyr
      tyrs = tini
      call monnc_comp

*     Check validity of initial state, before
*     proceeding with opening of output files
*     ---------------------------------------
      solnok = .true.
      call valids (solnok)
      if ( .not.solnok ) then
        print *,' *** valids has detected invalid values ***'
        write(*,'(a,i12,f12.2,f11.4)')
     &        '  problem occurs at nt, tday, tyrs = ',nt,tday,tyrs
        print *,' program will terminate after diagnostic printout'
        call prsamp (nt)
        call cfltry (nt)
        print *,' '
        print *,' program terminates'
        stop
      endif

#ifdef use_netcdf
*     netCDF output initialisation
*     ============================
*     N.B. if the error message "One or more variable sizes
*     violate format constraints" occurs when initialising
*     a netCDF file, try either coarsening spatial
*     and/or temporal sampling, or create that file with
*     "64-bit offset" format; see netCDF manuals for details
      print *,' '

#  ifndef atmos_only
*     Work out the size (in Gbytes) of a 4-D oceanic p-grid variable
*     (believed to be the maximal size variable in the netCDF dumps)
*     This assumes single precision (4 byte) dumps; if using
*     double precision the initial 4 should be changed to 8.
      poc4dsz = 4.0d0*(1 + nxto/nsko)*(1 + nyto/nsko)*nlo*noutstepoc
     &         /dble(1024**3)
      write(*,'(a,a,f7.3)') '  netCDF3: 4-dim oceanic p-grid',
     &                      ' record size (Gbytes) = ',poc4dsz
      if ( poc4dsz.ge.4.0d0 ) then
        print *,' netCDF3 ERROR: netCDF3 limit = 4.0 Gbytes'
        print *,' netCDF3 ERROR: program terminates'
        STOP
      endif
*     Ocean output files
*     ------------------
      status = nf_create (outdir(1:lenod)//'/ocpo.nc',
***  &                    NF_CLOBBER, ocpid)
*     Alternative creation call with "Large File Support"
     &                 IOR(NF_CLOBBER,NF_64BIT_OFFSET), ocpid)
      if ( status.ne.NF_NOERR ) then
        print *,' problem creating file ocpo.nc'
        call handle_err (status, subnam)
       else
        print *,' ocpo.nc file created'
      endif
      status = nf_create (outdir(1:lenod)//'/ocsst.nc',
     &                    NF_CLOBBER, octid)
*     Alternative creation call with "Large File Support"
***  &                    IOR(NF_CLOBBER,NF_64BIT_OFFSET), octid)
      if ( status.ne.NF_NOERR ) then
        print *,' problem creating file ocsst.nc'
        call handle_err (status, subnam)
       else
        print *,' ocsst.nc file created'
      endif
      call ocnc_init (nsko, outfloc)
      print *,' Ocean netCDF files initialised'
      call ocnc_out (nsko, outfloc)
#    ifdef qoc_diag
*     Optional ocean vorticity diagnostics
      status = nf_create (outdir(1:lenod)//'/qocdiag.nc',
***  &                    NF_CLOBBER, qocncid)
*     Alternative creation call with "Large File Support"
     &                    IOR(NF_CLOBBER,NF_64BIT_OFFSET), qocncid)
      if ( status.ne.NF_NOERR ) then
        print *,' problem creating file qocdiag.nc'
        call handle_err (status, subnam)
       else
        print *,' qocdiag.nc file created'
      endif
      call qocdiag_init (nsko)
      print *,' qocdiag netCDF file initialised'
#    endif
#  endif /* atmos_only */

#  ifndef ocean_only
*     Work out the size (in Gbytes) of a 4-D atmospheric p-grid variable
*     (believed to be the maximal size variable in the netCDF dumps)
*     This assumes single precision (4 byte) dumps; if using
*     double precision the initial 4 should be changed to 8.
      pat4dsz = 4.0d0*(1 + nxta/nska)*(1 + nyta/nska)*nla*noutstepat
     &         /dble(1024**3)
      write(*,'(a,a,f7.3)') '  netCDF3: 4-dim atmospheric p-grid',
     &                      ' record size (Gbytes) = ',pat4dsz
      if ( pat4dsz.ge.4.0d0 ) then
        print *,' netCDF3 ERROR: netCDF3 limit = 4.0 Gbytes'
        print *,' netCDF3 ERROR: program terminates'
        STOP
      endif
*     Atmosphere output files
*     -----------------------
      status = nf_create (outdir(1:lenod)//'/atpa.nc',
     &                    NF_CLOBBER, atpid)
*     Alternative creation call with "Large File Support"
***  &                    IOR(NF_CLOBBER,NF_64BIT_OFFSET), atpid)
      if ( status.ne.NF_NOERR ) then
        print *,' problem creating file atpa.nc'
        call handle_err (status, subnam)
       else
        print *,' atpa.nc file created'
      endif
      status = nf_create (outdir(1:lenod)//'/atast.nc',
     &                    NF_CLOBBER, attid)
*     Alternative creation call with "Large File Support"
***  &                    IOR(NF_CLOBBER,NF_64BIT_OFFSET), attid)
      if ( status.ne.NF_NOERR ) then
        print *,' problem creating file atast.nc'
        call handle_err (status, subnam)
       else
        print *,' atast.nc file created'
      endif
# ifdef moist
      status = nf_create (outdir(1:lenod)//'/athum.nc',
     &                    NF_CLOBBER, athid)
*     Alternative creation call with "Large File Support"
***  &                    IOR(NF_CLOBBER,NF_64BIT_OFFSET), athid)
      if ( status.ne.NF_NOERR ) then
        print *,' problem creating file athum.nc'
        call handle_err (status, subnam)
       else
        print *,' athum.nc file created'
      endif
# endif

      call atnc_init (nska, outflat)
      print *,' Atmos. netCDF files initialised'
      call atnc_out (nska, outflat)

#    ifdef highres_output
*     Work out the size (in Gbytes) of a 4-D atmospheric p-grid variable
*     (believed to be the maximal size variable in the netCDF dumps)
*     This assumes single precision (4 byte) dumps; if using
*     double precision the initial 4 should be changed to 8.
      pat4dsz = 4.0d0*(1 + nxtout/nska1)*(1 + nytout/nska1)
     &                                  *nla*noutstepat
     &         /dble(1024**3)
      write(*,'(a,a,f7.3)') ' netCDF3: 4-dim HR atmospheric p-grid',
     &                      ' record size (Gbytes) = ',pat4dsz
      if ( pat4dsz.ge.4.0d0 ) then
        print *,' netCDF3 ERROR: netCDF3 limit = 4.0 Gbytes'
        print *,' netCDF3 ERROR: program terminates'
        STOP
      endif

*     Atmosphere high-resolution output files (over ocean)
*     -----------------------
      status = nf_create (outdir(1:lenod)//'/atpahr.nc',
***     &                    NF_CLOBBER, atphrid)
*     Alternative creation call with "Large File Support"
     &          IOR(NF_CLOBBER,NF_64BIT_OFFSET), atphrid)
      if ( status.ne.NF_NOERR ) then
        print *,' problem creating file atpahr.nc'
        call handle_err (status, subnam)
       else
        print *,' atpahr.nc file created'
      endif
      status = nf_create (outdir(1:lenod)//'/atasthr.nc',
***     &                    NF_CLOBBER, atthrid)
*     Alternative creation call with "Large File Support"
     &          IOR(NF_CLOBBER,NF_64BIT_OFFSET), atthrid)
      if ( status.ne.NF_NOERR ) then
        print *,' problem creating file atasthr.nc'
        call handle_err (status, subnam)
       else
        print *,' atasthr.nc file created'
      endif
      call atnchr_init
      print *,' Atmos. high-res. netCDF files initialised'
      call atnchr_out
#    endif /* highres_output */

#  endif  /* not ocean_only */

*     Additional monitoring files
*     ---------------------------
      status = nf_create (outdir(1:lenod)//'/monit.nc',
     &                    NF_CLOBBER, monncid)
      if ( status.ne.NF_NOERR ) then
        print *,' problem creating file monit.nc'
        call handle_err (status, subnam)
       else
        print *,' monit.nc file created'
      endif
*     Ensure that monitoring and area-averaging are done
*     synchronously as post-processing routines assume this.
      call monnc_init
      print *,' Diagnostics netCDF file initialised'
      call monnc_out
#  ifdef get_areav
*     Area average output file
*     ------------------------
      status = nf_create (outdir(1:lenod)//'/areas.nc',
     &                    NF_CLOBBER, arencid)
*     Alternative creation call with "Large File Support"
***  &                    IOR(NF_CLOBBER,NF_64BIT_OFFSET), arencid)
      if ( status.ne.NF_NOERR ) then
        print *,' problem creating file areas.nc'
        call handle_err (status, subnam)
       else
        print *,' areas.nc file created'
      endif
#  endif
#endif

*     Print sample of initial state
*     -----------------------------
*     (validity of initial state already tested)
      call prsamp (nt)
#ifdef get_areav
      call areavg
#endif

*     Time averaging of forcing fields
      if ( avgatm .or. avgocn ) call tavini

#ifdef get_covar
      call covini
#endif /* get_covar */


****************    START OF MAIN TIMESTEPPING LOOP     ****************

      do 1000 nt=nsteps0+1,nsteps

        if ( mod(nt,nstr).eq.1 ) then

*         Compute forcing fields, if on ocean step
#ifndef ocean_only
          call xforc

#endif

#ifndef atmos_only
*         Step ocean mixed layer
          call oml
 
#    ifdef qoc_diag
*         Compute dq/dt diagnostics here, after entoc has
*         been computed, but before po and qo are updated
          if ( mod(ntdone,noutoc).eq.0 ) then
            call qocdiag_out (nsko)
          endif
#    endif

*         Step ocean qg model
          call qgostep

*         Invert pv into oceanic pressures
          call ocinvq

*         Compute pv on ocean boundaries (mixed condition)
          call ocqbdy (qo, po)

#endif /* not atmos_only */

        endif

#ifndef ocean_only
*       Step atmospheric mixed layer
        call aml

*       Step atmospheric qg channel model
        call qgastep

*       Invert pv into atmospheric pressures
        call atinvq

*       Compute pv on atmosphere zonal boundary (mixed condition)
        call atqzbd (qa, pa)
#endif

*       Timestep done; do checking and diagnostics as necessary
        ntdone = nt - nsteps0
        tday = nt*dta/secday
        tyrs = nt*dta/secsyr

*       Periodically check validity of solution
*       ---------------------------------------
        if ( mod(ntdone,nvalid).eq.0 ) then
          call valids (solnok)
          if ( .not.solnok ) then
            print *,' *** valids has detected invalid values ***'
            write(*,'(a,i12,f12.2,f11.4)')
     &            '  problem occurs at nt, tday, tyrs = ',nt,tday,tyrs
            print *,' program will terminate after diagnostic printout'
            call monnc_comp
            call prsamp (nt)
            call cfltry (nt)
#ifdef use_netcdf
            print *,' make netCDF dumps of fields before exit'

            call monnc_out
            status = nf_close(monncid)
            if ( status.ne.NF_NOERR ) call handle_err (status, subnam)
#  ifdef get_areav
            status = nf_close(arencid)
            if ( status.ne.NF_NOERR ) call handle_err (status, subnam)
#  endif

#  ifndef atmos_only
            call ocnc_out (nsko, outfloc)
            status = nf_close(ocpid)
            if ( status.ne.NF_NOERR ) call handle_err (status, subnam)
            status = nf_close(octid)
            if ( status.ne.NF_NOERR ) call handle_err (status, subnam)
#    ifdef qoc_diag
            call qocdiag_out (nsko)
            status = nf_close(qocncid)
            if ( status.ne.NF_NOERR ) call handle_err (status, subnam)
#    endif
#  endif /* atmos_only */

#  ifndef ocean_only
            call atnc_out (nska, outflat)
            status = nf_close(atpid)
            if ( status.ne.NF_NOERR ) call handle_err (status, subnam)
            status = nf_close(attid)
            if ( status.ne.NF_NOERR ) call handle_err (status, subnam)
# ifdef moist
            status = nf_close(athid)
            if ( status.ne.NF_NOERR ) call handle_err (status, subnam)
# endif
#    ifdef highres_output
            call atnchr_out
            status = nf_close(atphrid)
            if ( status.ne.NF_NOERR ) call handle_err (status, subnam)
            status = nf_close(atthrid)
            if ( status.ne.NF_NOERR ) call handle_err (status, subnam)
#    endif /* highres_output */

#  endif
#endif
            print *,' program terminates'
            stop
          endif
        endif

*       Suppress computational mode of leapfrog scheme
*       ==============================================
#ifndef atmos_only
        if ( mod(nt-1,25*nstr).eq.0 ) then
*         Average oceanic time levels
*         ---------------------------
!$OMP PARALLEL DEFAULT (NONE)
!$OMP&         PRIVATE (i,j,k)
!$OMP&         SHARED  (qo,qom,po,pom,sst,sstm)

          do k=1,nlo
!$OMP       DO SCHEDULE (STATIC)
            do j=1,nypo
              do i=1,nxpo
                qo(i,j,k) = 0.5d0*( qo(i,j,k)+qom(i,j,k) )
                po(i,j,k) = 0.5d0*( po(i,j,k)+pom(i,j,k) )
              enddo
            enddo
!$OMP       END DO
          enddo

!$OMP     DO SCHEDULE (STATIC)
          do j=1,nyto
            do i=1,nxto
              sst(i,j) = 0.5d0*( sst(i,j)+sstm(i,j) )
            enddo
          enddo
!$OMP     END DO NOWAIT

!$OMP     END PARALLEL

*         Also average constraint variables
          do k=1,nlo-1
            dpioc(k) = 0.5d0*( dpioc(k) + dpiocp(k) )
          enddo
#  ifdef cyclic_ocean
          do k=1,nlo
            ocncs(k) = 0.5d0*( ocncs(k) + ocncsp(k) )
            ocncn(k) = 0.5d0*( ocncn(k) + ocncnp(k) )
          enddo
#  endif
        endif
#endif /* atmos_only */

#ifndef ocean_only
        if ( mod(nt-1,100).eq.0 ) then
*         Average atmospheric time levels
*         -------------------------------
!$OMP PARALLEL DEFAULT (NONE)
!$OMP&         PRIVATE (i,j,k)
!$OMP&         SHARED  (qa,qam,pa,pam,ast,astm,hmixa,hmixam)
# ifdef moist
!$OMP&         SHARED  (hum,humm)
# endif

          do k=1,nla
!$OMP       DO SCHEDULE (STATIC)
            do j=1,nypa
              do i=1,nxpa
                qa(i,j,k) = 0.5d0*( qa(i,j,k)+qam(i,j,k) )
                pa(i,j,k) = 0.5d0*( pa(i,j,k)+pam(i,j,k) )
              enddo
            enddo
!$OMP       END DO

# ifdef moist

!$OMP       DO SCHEDULE (STATIC)
            do j=1,nypa
              do i=1,nxpa
               hum(i,j,k) = 0.5d0*( hum(i,j,k)+humm(i,j,k) )
              enddo
            enddo
!$OMP       END DO

# endif
          enddo

!$OMP     DO SCHEDULE (STATIC)
          do j=1,nyta
            do i=1,nxta
              ast(i,j) = 0.5d0*( ast(i,j)+astm(i,j) )
              hmixa(i,j) = 0.5d0*( hmixa(i,j)+hmixam(i,j) )
            enddo
          enddo
!$OMP     END DO NOWAIT

!$OMP     END PARALLEL

*         Also average constraint variables
          do k=1,nla-1
            dpiat(k) = 0.5d0*( dpiat(k) + dpiatp(k) )
          enddo
          do k=1,nla
            atmcs(k) = 0.5d0*( atmcs(k) + atmcsp(k) )
            atmcn(k) = 0.5d0*( atmcn(k) + atmcnp(k) )
          enddo
        endif
#endif

*       Occasionally dump restart file
*       ------------------------------
        if ( resdmp ) then
          if ( mod(ntdone,noutre).eq.0 ) then
*           First check the solution is valid
            call valids (solnok)
            if ( solnok ) then
              print *,' '
              print *,' Writing restart dump at nt, tyrs = ',nt,tyrs
*             Open and close restart file each time because
*             flush and rewind not yet working properly under Linux.
#ifdef use_netcdf
              status = nf_create (outdir(1:lenod)//'/restart.nc',
     &                            NF_CLOBBER, restid)
              if ( status.ne.NF_NOERR ) call handle_err (status, subnam)
              call resave_nc
              status = nf_close(restid)
              if ( status.ne.NF_NOERR ) call handle_err (status, subnam)
#else
              open (iorest, file=outdir(1:lenod)//'/restart',
     &              form='unformatted', status='unknown')
              call resave (iorest)
              close (iorest)
#endif
            endif
          endif
        endif

*       Monitoring of current state if required
*       ---------------------------------------
*       Ensure that monitoring and area-averaging are done
*       synchronously as post-processing routines assume this.
        if ( mod(ntdone,nocmon).eq.0 ) then
          call monnc_comp
#ifdef use_netcdf
          call monnc_out
#endif
#ifdef get_areav
          call areavg
#endif
        endif
        if ( mod(ntdone,nocprt).eq.0 ) then
          call prsamp (nt)
        endif

#ifdef use_netcdf
*       Write full fields to netCDF files
*       ---------------------------------
#  ifndef atmos_only
        if ( mod(ntdone,noutoc).eq.0 ) then
          call ocnc_out (nsko, outfloc)
        endif
#  endif /* atmos_only */
#  ifndef ocean_only
        if ( mod(ntdone,noutat).eq.0 ) then
          call atnc_out (nska, outflat)
#    ifdef highres_output
          call atnchr_out
#    endif /* highres_output */
        endif
#  endif
#endif

*       Add contributions to running means of forcing fields
*       ----------------------------------------------------
        if ( avgatm ) then
          if ( mod(ntdone,ntavat).eq.nmidat ) call tavatm
        endif
        if ( avgocn ) then
          if ( mod(ntdone,ntavoc).eq.nmidoc ) call tavocn
        endif

#ifdef get_covar
*       Accumulate contributions to covariance matrices
*       -----------------------------------------------
        if ( mod(ntdone,ntcovat).eq.0 ) call covatm
        if ( mod(ntdone,ntcovoc).eq.0 ) call covocn
#endif /* get_covar */

 1000 continue

****************     END OF MAIN TIMESTEPPING LOOP      ****************

#ifdef use_netcdf
*     Close monitoring, netCDF and optional restart files
      status = nf_close(monncid)
      if ( status.ne.NF_NOERR ) call handle_err (status, subnam)
#  ifdef get_areav
      status = nf_close(arencid)
      if ( status.ne.NF_NOERR ) call handle_err (status, subnam)
#  endif

#  ifndef atmos_only
      status = nf_close(ocpid)
      if ( status.ne.NF_NOERR ) call handle_err (status, subnam)
      status = nf_close(octid)
      if ( status.ne.NF_NOERR ) call handle_err (status, subnam)
#    ifdef qoc_diag
*     Compute qoc diagnostics on final state. Ideally oml would also be
*     called here to update entoc, but this would wrongly change sst and sstm
      call qocdiag_out (nsko)
      status = nf_close(qocncid)
      if ( status.ne.NF_NOERR ) call handle_err (status, subnam)
#    endif
#  endif /* atmos_only */

#  ifndef ocean_only
      status = nf_close(atpid)
      if ( status.ne.NF_NOERR ) call handle_err (status, subnam)
      status = nf_close(attid)
      if ( status.ne.NF_NOERR ) call handle_err (status, subnam)
# ifdef moist
      status = nf_close(athid)
      if ( status.ne.NF_NOERR ) call handle_err (status, subnam)
# endif
#  endif
#endif

*     dump restart file (either netCDF or unformatted)
#ifdef use_netcdf
      status = nf_create (outdir(1:lenod)//'/lastday.nc',
     &                    NF_CLOBBER, restid)
      if ( status.ne.NF_NOERR ) call handle_err (status, subnam)
      call resave_nc
      status = nf_close(restid)
      if ( status.ne.NF_NOERR ) call handle_err (status, subnam)
#else
      open (14, file=outdir(1:lenod)//'/last.day', form='unformatted',
     &      status='unknown')
      call resave (14)
      close(14)
#endif

*     Compute and output time-averaged fields
#ifdef use_netcdf
      status = nf_create (outdir(1:lenod)//'/avges.nc',
     &                    NF_CLOBBER, tavncid)
      if ( status.ne.NF_NOERR ) call handle_err (status, subnam)
#endif
      call tavout
#ifdef use_netcdf
      status = nf_close(tavncid)
      if ( status.ne.NF_NOERR ) call handle_err (status, subnam)
#endif

#ifdef get_covar
*     Write out covariance matrix
#  ifdef use_netcdf
      status = nf_create (outdir(1:lenod)//'/covar.nc',
     &                    NF_CLOBBER, covncid)
      if ( status.ne.NF_NOERR ) call handle_err (status, subnam)
#  endif
      call covout
#  ifdef use_netcdf
      status = nf_close(covncid)
      if ( status.ne.NF_NOERR ) call handle_err (status, subnam)
#  endif
#endif /* get_covar */

*     Precautionary call to prsamp, if not already done
      if ( mod(ntdone,nocprt).ne.0 ) then
        call prsamp (nt)
      endif

*     Compute final mass and momentum integrals for constraints
*     ---------------------------------------------------------
      call constr

      print *,' '
      write(*,'(a,i12,f11.2,f11.4)')
     &        '  End of run at nt, tday, tyrs = ',nt,tday,tyrs

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
*
      IMPLICIT NONE
*
*     Subroutine arguments
      character (len=80) :: buffer
      integer, INTENT(IN) :: iounit
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
      USE atconst
#endif
#ifndef atmos_only
      USE ocstate
#endif
      USE intrfac, ONLY : sst, sstm, ast, astm, hmixa, hmixam, 
     &                   hmat, xfa0

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
# ifdef moist
!$OMP&         SHARED  (hum,humm,qatm,qat,xfa0)
# endif

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


*     Initialise atmos. mixed layer rel. temp. and humidity
*
!$OMP DO SCHEDULE (STATIC)
      do j=1,nyta
        do i=1,nxta
          ast(i,j) = 0.0d0
          astm(i,j) = 0.0d0
# ifndef moist
            hmixa(i,j) = hmat
            hmixam(i,j)= hmat
# else

            hmixa(i,j) = qatm
            hmixam(i,j)= qatm
            hum(i,j,1)= qat(1)
            humm(i,j,1)=qat(1)
            hum(i,j,2)=qat(2)
            humm(i,j,2)=qat(2)
            hum(i,j,3)=qat(3)
            humm(i,j,3)=qat(3)
            xfa0(i,j)=0.d0

# endif
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
#endif /* atmos_only */

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
#ifndef ocean_only
      USE atstate
#endif /* not ocean_only */
#ifndef atmos_only
      USE ocstate
#endif /* not atmos_only */
      USE intrfac, ONLY : sst, sstm, sstbar, ast, astm, astbar,
     &                    hmixa, hmixam, hmat, xfa0
      USE xfosubs, ONLY : fsprim
      USE radiate

      IMPLICIT NONE

*     Subroutine arguments
*
*     Local parameters
      logical prtval
      parameter ( prtval = .false. )
*
*     Local variables
      integer i,j,k,nxco,nyco,nxca,nyca
      double precision plfac(nla),play
**    double precision sumc,sumcp,p1off

*     Centres of domains if needed for sample printing
      nxco = (nxpo+1)/2
      nyco = (nypo+1)/2
      nxca = (nxpa+1)/2
      nyca = (nypa+1)/2

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
        write (*,*) ' Initial relative ast:'
        do j=nyta,1,-1
          if ( astbar(j).gt.tat(1) ) then
            write (*,'(i6,f16.8,a)') j,astbar(j),'  convect'
           else
            write (*,'(i6,f16.8,a)') j,astbar(j)
          endif
        enddo
      endif

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

*     Initialise atmosphere fields (incl. mixed layer)
*     ------------------------------------------------
!$OMP PARALLEL DEFAULT (NONE)
!$OMP&         PRIVATE (i,j,k,play)
!$OMP&         SHARED  (ast,astm,astbar,hmixa,hmixam,hmat)
#ifndef ocean_only
!$OMP&         SHARED  (pa,pam,plfac,yparel)
# ifdef moist
!$OMP&         SHARED  (hum,humm,qatm,qat,xfa0)
# endif

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
#endif /* not ocean_only */

*     Initialise atmos. mixed layer rel. temp. and humidity
*     The a.m.l is initialised even in the ocean_only
*     case, to provide an upper boundary condition

!$OMP DO SCHEDULE (STATIC)
      do j=1,nyta
        do i=1,nxta
          ast(i,j) = astbar(j)
          astm(i,j) = astbar(j)
# ifndef moist
            hmixa(i,j) = hmat
            hmixam(i,j)= hmat
# else

            hmixa(i,j) = qatm
            hmixam(i,j)= qatm
            hum(i,j,1)= qat(1)
            humm(i,j,1)=qat(1)
            hum(i,j,2)=qat(2)
            humm(i,j,2)=qat(2)
            hum(i,j,3)=qat(3)
            humm(i,j,3)=qat(3)
            xfa0(i,j)=0.d0
# endif
        enddo
      enddo
!$OMP END DO NOWAIT

!$OMP END PARALLEL

#ifndef ocean_only
      if ( prtval ) then
        print *,' '
*       Print meridional sections at the centre of the domain,
*       since the initial values are a function of y only
        write(*,*) ' Meridional sections of initial',
     &             ' values across atmos. centre'
        do k=1,nla
          print *,' '
          write(*,'(a,i2,a)') '  Layer k = ',k,':'
          write(*,'(a,a)') '     j        pa(k)            pam(k)'
          do j=nypa,1,-1
            write(*,240) j,pa(nxca,j,k),pam(nxca,j,k)
          enddo
        enddo
      endif
#endif /* not ocean_only */

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
#endif /* not atmos_only */

*     Initialise oceanic mixed layer rel. temp.
*     The o.m.l is initialised even in the atmos_only
*     case, to provide a lower boundary condition
!$OMP DO SCHEDULE (STATIC)
      do j=1,nyto
        do i=1,nxto
          sst(i,j) = sstbar(j)
          sstm(i,j) = sstbar(j)
        enddo
      enddo
!$OMP END DO NOWAIT

!$OMP END PARALLEL

#ifndef atmos_only
      if ( prtval ) then
        print *,' '
*       Print meridional sections at the centre of the domain,
*       since the initial values are a function of y only
        write(*,*) ' Meridional sections of initial',
     &             ' values across atmos. centre'
        do k=1,nlo
          print *,' '
          write(*,'(a,i2,a)') '  Layer k = ',k,':'
          write(*,'(a,a)') '     j        po(k)            pom(k)'
          do j=nypo,1,-1
            write(*,240) j,po(nxco,j,k),pom(nxco,j,k)
          enddo
        enddo
      endif
#endif /* not atmos_only */

      print *,' '
      write (*,*) ' Rbal initialisation complete'

  240 format(i6,1p,5d18.10)

      END SUBROUTINE rbalin
c
c***********************************************************************
c
      SUBROUTINE prsamp (nt)
*
*     Print out a sample of the current state. Useful for
*     tracking progress, or for testing/debugging purposes

*     Modules
      USE parameters, ONLY : nxpo, nypo, nxto, nyto, nlo,
     &                       nxpa, nypa, nxta, nyta, nla, fnot
#ifndef ocean_only
      USE atconst, ONLY : ddynat
      USE atstate, ONLY : pa, qa, wekpa
#endif
#ifndef atmos_only
      USE occonst, ONLY : ddynoc
      USE ocstate, ONLY : po, qo, wekpo
#endif
      USE intrfac, ONLY : sst, ast, hmixa
      USE timinfo, ONLY : tday, tyrs
      USE monitor, ONLY : pavgoc, qavgoc, kealoc, et2moc, osfmin,
     &                    osfmax, occirc, kealat, et2mat

      IMPLICIT NONE

*     Subroutine arguments
      integer, INTENT(IN) :: nt
*
*     Local parameters
      double precision extrem
      parameter ( extrem=1.0d30 )
*
*     Local variables
      integer i,j,k
#ifndef atmos_only
      integer nxco,nyco
      double precision sstmin,sstmax,poref(nlo)
#endif /* not atmos_only */
#ifndef ocean_only
      integer nxca,nyca
      double precision astmin,astmax,hmxmin,hmxmax,paref(nla)
#endif /* not ocean_only */

*     Centres of domains for sample printing
#ifndef atmos_only
      nxco = (nxpo+1)/2
      nyco = (nypo+1)/2
#endif /* not atmos_only */
#ifndef ocean_only
      nxca = (nxpa+1)/2
      nyca = (nypa+1)/2
#endif /* not ocean_only */

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
#endif /* not ocean_only */
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
#endif /* not atmos_only */

      print *,' '
      write(*,'(a,i12,f11.2,f11.4)')
     &      '  Sample output at nt, tday, tyrs = ',nt,tday,tyrs

*     Print some ocean spot values and extrema
*     ----------------------------------------
#ifndef atmos_only
      write(*,217) '  po(k) at centre = ',(po(nxco,nyco,k),k=1,nlo)
      write(*,217) '  qo(k) at centre = ',(qo(nxco,nyco,k),k=1,nlo)
**    write(*,217) '  po(1) on S. bdy = ',(po(i,  1 ,1),i=nxco-1,nxco+1)
**    write(*,217) '  qo(1) on S. bdy = ',(qo(i,  1 ,1),i=nxco-1,nxco+1)
**    write(*,217) '  po(1) on N. bdy = ',(po(i,nypo,1),i=nxco-1,nxco+1)
**    write(*,217) '  qo(1) on N. bdy = ',(qo(i,nypo,1),i=nxco-1,nxco+1)
**    write(*,217) '   sst  on S. bdy = ',(sst(i,  1 ),i=nxco-1,nxco+1)
**    write(*,217) '   sst  on N. bdy = ',(sst(i,nyto),i=nxco-1,nxco+1)
      write(*,217) '  po(k)  average  = ',(pavgoc(k),k=1,nlo)
      write(*,217) '  qo(k)  average  = ',(qavgoc(k),k=1,nlo)
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
#endif /* atmos_only */

*     Print some atmos. spot values and extrema
*     -----------------------------------------
#ifndef ocean_only
      write(*,217) '  pa(k) at centre = ',(pa(nxca,nyca,k),k=1,nla)
      write(*,217) '  qa(k) at centre = ',(qa(nxca,nyca,k),k=1,nla)
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
#endif /* not ocean_only */

  217 format(a,1p,9d15.7)
  240 format(i6,1p,5d18.10)

      END SUBROUTINE prsamp
c
c***********************************************************************
c
      SUBROUTINE cfltry (nt)
*
*     Check how well the current state satisfies the CFL criterion.
*     From version 1.5.0, additionally checks atmos. mixed layer.

*     Modules
      USE parameters, ONLY : nxpo, nypo, nlo, nxpa, nypa, nla, fnot
#ifndef ocean_only
      USE atconst, ONLY : dta, dxa, rdxaf0
      USE atstate, ONLY : pa
      USE intrfac, ONLY : uekat, vekat
#endif /* not ocean_only */
#ifndef atmos_only
      USE occonst, ONLY : dto, dxo, rdxof0
      USE ocstate, ONLY : po
      USE intrfac, ONLY : tauxo, tauyo, hmoc
#endif /* not atmos_only */
      USE timinfo, ONLY : tday, tyrs

      IMPLICIT NONE

*     Subroutine arguments
      integer, INTENT(IN) :: nt
*
*     Local parameters
      double precision extrem,cflcrit
      logical cflocn,cflatm,cflmix
      parameter ( extrem=1.0d30, cflocn=.true., cflatm=.true.,
     &            cflmix= .true., cflcrit=0.8d0 )
*
*     Local variables
      integer i,j,k
      double precision uabs,uabmax,vabs,vabmax
#ifndef atmos_only
      double precision uamaxo(nlo),vamaxo(nlo),rhf0hm
#endif /* not atmos_only */
#ifndef ocean_only
      double precision uamaxa(nla),vamaxa(nla)
#endif /* not ocean_only */

      print *,' '
      write(*,'(a,i12,f11.2,f11.4)')
     &      '  CFL testing   at nt, tday, tyrs = ',nt,tday,tyrs

#ifndef atmos_only

*     Check extreme ocean velocity components and derive Courant numbers
*     ==================================================================
      if ( cflocn ) then

        if ( cflmix ) then
*         Optionally also infer um, vm including the
*         ageostrophic components (Ekman transport)
          rhf0hm = 0.5d0/(fnot*hmoc)
          uabmax = -extrem
          vabmax = -extrem
!$OMP     PARALLEL DEFAULT (NONE)
!$OMP&             PRIVATE (i,j,uabs,vabs)
!$OMP&             SHARED  (rdxof0,rhf0hm,po,tauxo,tauyo,uabmax,vabmax)

!$OMP     DO SCHEDULE (STATIC)
!$OMP&       REDUCTION(max:uabmax)
          do j=1,nypo-1
            do i=1,nxpo
              uabs = abs( -rdxof0*(  po(i,j+1,1) -  po(i,j,1) )
     &                    +rhf0hm*( tauyo(i,j+1) + tauyo(i,j) ) )
              uabmax = max( uabmax, uabs )
            enddo
          enddo
!$OMP     END DO NOWAIT
!$OMP     DO SCHEDULE (STATIC)
!$OMP&       REDUCTION(max:vabmax)
          do j=1,nypo
            do i=1,nxpo-1
              vabs = abs(  rdxof0*(  po(i+1,j,1) -  po(i,j,1) )
     &                    -rhf0hm*( tauxo(i+1,j) + tauxo(i,j) ) )
              vabmax = max( vabmax, vabs )
            enddo
          enddo
!$OMP     END DO NOWAIT

!$OMP     END PARALLEL

          write(*,217) '  Ocn. m.l. |u| max, CFL = ',
     &                 uabmax,uabmax*dto/dxo
          write(*,217) '  Ocn. m.l. |v| max, CFL = ',
     &                 vabmax,vabmax*dto/dxo
*         If have bad CFL values, scan for and print locations
          if ( uabmax*dto/dxo.ge.cflcrit ) then
            do j=1,nypo-1
              do i=1,nxpo
                uabs = (dto/dxo)*
     &                 abs( -rdxof0*(  po(i,j+1,1) -  po(i,j,1) )
     &                      +rhf0hm*( tauyo(i,j+1) + tauyo(i,j) ) )
                if ( uabs.ge.cflcrit ) then
                  write(*,250) '  Bad m.l. |u|; CFL, i, j = ',uabs,i,j
                endif
              enddo
            enddo
          endif
          if ( vabmax*dto/dxo.ge.cflcrit ) then
            do j=1,nypo
              do i=1,nxpo-1
                vabs = (dto/dxo)*
     &                 abs(  rdxof0*(  po(i+1,j,1) -  po(i,j,1) )
     &                      -rhf0hm*( tauxo(i+1,j) + tauxo(i,j) ) )
                if ( vabs.ge.cflcrit ) then
                  write(*,250) '  Bad m.l. |v|; CFL, i, j = ',vabs,i,j
                endif
              enddo
            enddo
          endif
        endif

*       Infer u, v in Q-G layers geostrophically
*       ----------------------------------------
        do k=1,nlo
          uabmax = -extrem
          vabmax = -extrem
!$OMP     PARALLEL DEFAULT (NONE)
!$OMP&             PRIVATE (i,j,uabs,vabs)
!$OMP&             SHARED  (k,rdxof0,po,uabmax,vabmax)

!$OMP     DO SCHEDULE (STATIC)
!$OMP&       REDUCTION(max:uabmax)
          do j=1,nypo-1
            do i=1,nxpo
              uabs = abs( -rdxof0*( po(i,j+1,k) - po(i,j,k) ) )
              uabmax = max( uabmax, uabs )
            enddo
          enddo
!$OMP     END DO NOWAIT
!$OMP     DO SCHEDULE (STATIC)
!$OMP&       REDUCTION(max:vabmax)
          do j=1,nypo
            do i=1,nxpo-1
              vabs = abs( rdxof0*( po(i+1,j,k) - po(i,j,k) ) )
              vabmax = max( vabmax, vabs )
            enddo
          enddo
!$OMP     END DO NOWAIT

!$OMP     END PARALLEL
          uamaxo(k) = uabmax
          vamaxo(k) = vabmax
        enddo
        write(*,217) '  Ocn. max geos |u|(k) = ',(uamaxo(k),k=1,nlo)
        write(*,217) '  Ocn. CFL geos |u|(k) = ',
     &               (uamaxo(k)*dto/dxo,k=1,nlo)
        write(*,217) '  Ocn. max geos |v|(k) = ',(vamaxo(k),k=1,nlo)
        write(*,217) '  Ocn. CFL geos |v|(k) = ',
     &               (vamaxo(k)*dto/dxo,k=1,nlo)
*       If have bad CFL values, scan for and print locations
        do k=1,nlo
          if ( uamaxo(k)*dto/dxo.ge.cflcrit ) then
            do j=1,nypo-1
              do i=1,nxpo
                uabs = (dto/dxo)*
     &                 abs( -rdxof0*( po(i,j+1,k) - po(i,j,k) ) )
                if ( uabs.ge.cflcrit ) then
                  write(*,250) '  Bad |u|; CFL, i, j, k = ',uabs,i,j,k
                endif
              enddo
            enddo
          endif
          if ( vamaxo(k)*dto/dxo.ge.cflcrit ) then
            do j=1,nypo
              do i=1,nxpo-1
                vabs = (dto/dxo)*
     &                 abs( rdxof0*( po(i+1,j,k) - po(i,j,k) ) )
                if ( vabs.ge.cflcrit ) then
                  write(*,250) '  Bad |v|; CFL, i, j, k = ',vabs,i,j,k
                endif
              enddo
            enddo
          endif
        enddo
      endif

#endif

#ifndef ocean_only

*     Check extreme atmos. velocity components and derive Courant numbers
*     ===================================================================
      if ( cflatm ) then
        if ( cflmix ) then
*         Optionally also infer um, vm including the ageostrophic
*         components (Ekman transport), which as uekat,
*         vekat, have already been computed in xfosubs.F
          uabmax = -extrem
          vabmax = -extrem
!$OMP     PARALLEL DEFAULT (NONE)
!$OMP&             PRIVATE (i,j,uabs,vabs)
!$OMP&             SHARED  (rdxaf0,pa,uekat,vekat,uabmax,vabmax)

!$OMP     DO SCHEDULE (STATIC)
!$OMP&       REDUCTION(max:uabmax)
          do j=1,nypa-1
            do i=1,nxpa
              uabs = abs( -rdxaf0*( pa(i,j+1,1) -  pa(i,j,1) )
     &                   + uekat(i,j) )
              uabmax = max( uabmax, uabs )
            enddo
          enddo
!$OMP     END DO NOWAIT
!$OMP     DO SCHEDULE (STATIC)
!$OMP&       REDUCTION(max:vabmax)
          do j=1,nypa
            do i=1,nxpa-1
              vabs = abs(  rdxaf0*( pa(i+1,j,1) -  pa(i,j,1) )
     &                   + vekat(i,j) )
              vabmax = max( vabmax, vabs )
            enddo
          enddo
!$OMP     END DO NOWAIT

!$OMP     END PARALLEL
          write(*,217) '  Atm. m.l. |u| max, CFL = ',
     &                 uabmax,uabmax*dta/dxa
          write(*,217) '  Atm. m.l. |v| max, CFL = ',
     &                 vabmax,vabmax*dta/dxa
*         If have bad CFL values, scan for and print locations
          if ( uabmax*dta/dxa.ge.cflcrit ) then
            do j=1,nypa-1
              do i=1,nxpa
                uabs = (dta/dxa)*
     &                 abs( -rdxaf0*( pa(i,j+1,1) -  pa(i,j,1) )
     &                     + uekat(i,j) )
                if ( uabs.ge.cflcrit ) then
                  write(*,250) '  Bad m.l. |u|; CFL, i, j = ',uabs,i,j
                endif
              enddo
            enddo
          endif
          if ( vabmax*dta/dxa.ge.cflcrit ) then
            do j=1,nypa
              do i=1,nxpa-1
                vabs = (dta/dxa)*
     &                 abs(  rdxaf0*( pa(i+1,j,1) -  pa(i,j,1) )
     &                     + vekat(i,j) )
                if ( vabs.ge.cflcrit ) then
                  write(*,250) '  Bad m.l. |v|; CFL, i, j = ',vabs,i,j
                endif
              enddo
            enddo
          endif
        endif

*       Infer u, v in Q-G layers geostrophically
*       ----------------------------------------
        do k=1,nla
          uabmax = -extrem
          vabmax = -extrem
!$OMP     PARALLEL DEFAULT (NONE)
!$OMP&             PRIVATE (i,j,uabs,vabs)
!$OMP&             SHARED  (k,rdxaf0,pa,uabmax,vabmax)

!$OMP     DO SCHEDULE (STATIC)
!$OMP&       REDUCTION(max:uabmax)
          do j=1,nypa-1
            do i=1,nxpa
              uabs = abs( -rdxaf0*( pa(i,j+1,k) - pa(i,j,k) ) )
              uabmax = max( uabmax, uabs )
            enddo
          enddo
!$OMP     END DO NOWAIT
!$OMP     DO SCHEDULE (STATIC)
!$OMP&       REDUCTION(max:vabmax)
          do j=1,nypa
            do i=1,nxpa-1
              vabs = abs( rdxaf0*( pa(i+1,j,k) - pa(i,j,k) ) )
              vabmax = max( vabmax, vabs )
            enddo
          enddo
!$OMP     END DO NOWAIT

!$OMP     END PARALLEL
          uamaxa(k) = uabmax
          vamaxa(k) = vabmax
        enddo
        write(*,217) '  Atm. max geos |u|(k) = ',(uamaxa(k),k=1,nla)
        write(*,217) '  Atm. CFL geos |u|(k) = ',
     &               (uamaxa(k)*dta/dxa,k=1,nla)
        write(*,217) '  Atm. max geos |v|(k) = ',(vamaxa(k),k=1,nla)
        write(*,217) '  Atm. CFL geos |v|(k) = ',
     &               (vamaxa(k)*dta/dxa,k=1,nla)
*       If have bad CFL values, scan for and print locations
        do k=1,nla
          if ( uamaxa(k)*dta/dxa.ge.cflcrit ) then
            do j=1,nypa-1
              do i=1,nxpa
                uabs = (dta/dxa)*
     &                 abs( -rdxaf0*( pa(i,j+1,k) - pa(i,j,k) ) )
                if ( uabs.ge.cflcrit ) then
                  write(*,250) '  Bad |u|; CFL, i, j, k = ',uabs,i,j,k
                endif
              enddo
            enddo
          endif
          if ( vamaxa(k)*dta/dxa.ge.cflcrit ) then
            do j=1,nypa
              do i=1,nxpa-1
                vabs = (dta/dxa)*
     &                 abs( rdxaf0*( pa(i+1,j,k) - pa(i,j,k) ) )
                if ( vabs.ge.cflcrit ) then
                  write(*,250) '  Bad |v|; CFL, i, j, k = ',vabs,i,j,k
                endif
              enddo
            enddo
          endif
        enddo
      endif

#endif /* not ocean_only */

  217 format(a,1p,9d15.7)
  250 format(a,1p,d15.7,2i8,i6)

      END SUBROUTINE cfltry
c
c***********************************************************************
c
      SUBROUTINE memreq

*     Compute size of ocean and atmosphere arrays, data module sizes,
*     and local memory allocation in subroutines if significant
*     and required (controlled by internal parameter stacks)

*     Modules
      USE parameters
#ifndef ocean_only
      USE atisubs
#endif
#ifndef atmos_only
      USE ocisubs
#endif
*     Use here those *.mod files of program modules
*     which define array-dimensioning parameters

      IMPLICIT NONE

*     Subroutine arguments
*
*     Local parameters
      double precision dim2mb
      logical stacks
      parameter ( dim2mb=1.0d0/131072.0d0, stacks = .false. )
*     dim2mb converts from no. of elements in a double
*     precision array to size in Megabytes (= 8/(1024**2) )
*
*     Local variables
      double precision sizarr,dynsiz,dgnsiz
      double precision sizmin

*     Compute sizes of main arrays
*     ============================
      print *,' '
      write(*,'(a)') '  Array sizes:'
      write(*,'(a)') '  ------------'
      sizarr = nxpo*nypo
      write(*,203) '  Ocean  2D p-array size (MB) = ',dim2mb*sizarr
      sizarr = nxpa*nypa
      write(*,203) '  Atmos. 2D p-array size (MB) = ',dim2mb*sizarr
*     Compute size of full 3D arrays
#ifndef atmos_only
      sizarr = nxpo*nypo*nlo
      write(*,203) '  Ocean  3D p-array size (MB) = ',dim2mb*sizarr
#endif
#ifndef ocean_only
      sizarr = nxpa*nypa*nla
      write(*,203) '  Atmos. 3D p-array size (MB) = ',dim2mb*sizarr
#endif
      sizarr = nxpaor*nypaor
      write(*,203) '  At higher (ocean) resolution,'
      write(*,203) '  atmos. 2D p-array size (MB) = ',dim2mb*sizarr
#ifdef get_covar
*     Compute size of covariance arrays
#  ifndef atmos_only
      sizarr = nmcvoc
      write(*,203) '  Ocean covar array size (MB) = ',dim2mb*sizarr
      sizarr = nvcvoc
      write(*,203) '  Ocean avge  array size (MB) = ',dim2mb*sizarr
#  endif
#  ifndef ocean_only
      sizarr = nmcvat
      write(*,203) '  Atmos covar array size (MB) = ',dim2mb*sizarr
      sizarr = nvcvat
      write(*,203) '  Atmos avge  array size (MB) = ',dim2mb*sizarr
#  endif
      write(*,201) ' '
#  ifndef atmos_only
      write(*,201) '  Ocean covar array dimension = ',nmcvoc
      write(*,201) '  Ocean avge  array dimension = ',nvcvoc
#  endif
#  ifndef ocean_only
      write(*,201) '  Atmos covar array dimension = ',nmcvat
      write(*,201) '  Atmoc avge  array dimension = ',nvcvat
#  endif
#endif /* get_covar */

*     Compute sizes of data storage in modules
*     ========================================
*     Accumulate the number of double precision variables
      dynsiz = 0.0d0
      dgnsiz = 0.0d0
      print *,' '
      write(*,'(a)') '  Data sizes of modules:'
      write(*,'(a)') '  ======================'
      write(*,'(a)') '  Dynamics:'
      write(*,'(a)') '  ---------'
      write(*,'(a)') '  Data modules:'
*     atconst data module
*     -------------------
      sizarr = 7 + nxpa + 2*nypa + nxta + 2*nyta + 2 + 3*nla + 4
     &        + 2*nla + nxta*nyta + 1 + nla*nla + 3*nla + 2*nla*nla
     &        + 2*nxpa*nypa + 1
      dynsiz = dynsiz + sizarr
      write(*,203) '  atconst mod total size (MB) = ',dim2mb*sizarr
*     athomog data module
*     -------------------
#ifndef ocean_only
      sizarr = 2*nypa*(nla-1) + 5*(nla-1) + nypa + 2 + 3*(nla-1)
     &        + 4*nla + 2*(nla-1) + 4*nla + 2
      dynsiz = dynsiz + sizarr
      write(*,203) '  athomog mod total size (MB) = ',dim2mb*sizarr
#endif
*     atstate data module
*     -------------------
      sizarr = 4*nxpa*nypa*nla + 2*nxpa*nypa + nxta*nyta
      dynsiz = dynsiz + sizarr
      write(*,203) '  atstate mod total size (MB) = ',dim2mb*sizarr
*     intrfac data module
*     -------------------
      sizarr = 2*nxto*nyto + nyto + 2*nxta*nyta + nyta + 2*nxpo*nypo
     &        + 2*nxpa*nypa + nxpa*nyta + nxta*nypa + nxpa*(nypa+1)
     &        + (nxpa+1)*nypa + nxto*nyto + 3*nxta*nyta + 13
      dynsiz = dynsiz + sizarr
      write(*,203) '  intrfac mod total size (MB) = ',dim2mb*sizarr
*     occonst data module
*     -------------------
      sizarr = 7 + nxpo + 2*nypo + nxto + 2*nyto + 2 + (nlo-1)
     &        + 2*nlo + 2 + 3*nlo + 5 + nlo*nlo + 3*nlo + 2*nlo*nlo
     &        + 2*nxpo*nypo + 1
      dynsiz = dynsiz + sizarr
      write(*,203) '  occonst mod total size (MB) = ',dim2mb*sizarr
*     ochomog data module
*     -------------------
#ifndef atmos_only
#  ifdef cyclic_ocean
*     Zonally cyclic ocean
      sizarr = 2*nypo*(nlo-1) + 4*(nlo-1) + nypo + (nlo-1) + 2
     &        + 4*nlo + 2*(nlo-1) + 6*nlo + 4
#  else
*     Finite box ocean
      sizarr = nxpo*nypo*(nlo-1) + (nlo-1) + nlo*(nlo-1)
     &        + 2*(nlo-1)*(nlo-1) + (nlo-1)/2
#  endif
      sizarr = sizarr + 3*(nlo-1)
      dynsiz = dynsiz + sizarr
      write(*,203) '  ochomog mod total size (MB) = ',dim2mb*sizarr
#endif
*     ocstate data module
*     -------------------
#ifndef atmos_only
      sizarr = 4*nxpo*nypo*nlo + 2*nxpo*nypo + nxto*nyto
      dynsiz = dynsiz + sizarr
      write(*,203) '  ocstate mod total size (MB) = ',dim2mb*sizarr
#endif
*     radiate data module
*     -------------------
      sizarr = 4 + nla + 4 + 2*nla*(nla-1) + 4 + 2*nla
     &        + 3 + nla + 2*(nla-1) + 3
      dynsiz = dynsiz + sizarr
      write(*,203) '  radiate mod total size (MB) = ',dim2mb*sizarr
*     timinfo data module
*     -------------------
      sizarr = 3 + 14 + (17+1)/2
#ifdef get_covar
      sizarr = sizarr + 1
#endif
      dynsiz = dynsiz + sizarr
      write(*,203) '  timinfo mod total size (MB) = ',dim2mb*sizarr
      write(*,'(a)') '  Program modules:'
*     atisubs prog module
*     -------------------
#ifndef ocean_only
      sizarr = lwftat + 1 + nxta
      dynsiz = dynsiz + sizarr
      write(*,203) '  atisubs mod total size (MB) = ',dim2mb*sizarr
#endif
*     ocisubs prog module
*     -------------------
#ifndef atmos_only
      sizarr = lwftoc + 1 + nxto
      dynsiz = dynsiz + sizarr
      write(*,203) '  ocisubs mod total size (MB) = ',dim2mb*sizarr
#endif
!!    [UP TO HERE]
*     xfosubs prog module
*     -------------------
      sizarr = 5*16*16 + 16*(ndxr+1)*(ndxr+1)
      dynsiz = dynsiz + sizarr
      write(*,203) '  xfosubs mod total size (MB) = ',dim2mb*sizarr

      if ( dim2mb*dynsiz.lt.10000.00 ) then
        write(*,'(a)') '                                  --------'
        write(*,203) '       Dynamics subtotal (MB) = ',dim2mb*dynsiz
        write(*,'(a)') '                                  --------'
       else
        write(*,'(a)') '                                 ---------'
        write(*,203) '       Dynamics subtotal (MB) = ',dim2mb*dynsiz
        write(*,'(a)') '                                 ---------'
      endif
      write(*,'(a)') '  Diagnostics:'
#ifdef get_covar
*     covaria prog module
*     -------------------
      sizarr = 0.5d0
#  ifndef ocean_only
      sizarr = sizarr + 2*( nmcvat + nvcvat + 1 ) + 1
#  endif
#  ifndef atmos_only
      sizarr = sizarr + 2*( nmcvoc + nvcvoc + 1 ) + 1
#  endif
      dgnsiz = dgnsiz + sizarr
      write(*,203) '  covaria mod total size (MB) = ',dim2mb*sizarr
#endif /* get_covar */
*     monitor prog module
*     -------------------
      sizarr = 9 + (nla-1) + 1 + (nla-1) + (nlo-1) + (nla-1) + 17
     &        + 2*(nla-1) + 1 + 2*(nlo-1) + 3*nlo + 3*nla + 2
     &        + 2*nlo + nla + (nlo-1)+ 2*(nla-1) + 5
     &        + (nla-1) + (nlo-1) + 3*nlo + nla + ( nla + nlo )/2
      dgnsiz = dgnsiz + sizarr
      write(*,203) '  monitor mod total size (MB) = ',dim2mb*sizarr
*     nc_subs prog module
*     -------------------
*     nc_subs: non-stack storage size is negligible
*     timavge prog module
*     -------------------
      sizarr = 0.5d0
#ifndef ocean_only
      sizarr =  sizarr + 2*nxpa*nypa + 3*nxta*nyta
     &        + 2*nxpa*nypa*nla + 3*nxpa*nyta + 3*nxta*nypa
#endif
#ifndef atmos_only
      sizarr =  sizarr + 3*nxpo*nypo + 3*nxto*nyto
     &        + 2*nxpo*nypo*nlo + 3*nxpo*nyto + 3*nxto*nypo + 2
#endif
      dgnsiz = dgnsiz + sizarr
      write(*,203) '  timavge mod total size (MB) = ',dim2mb*sizarr
      if ( dim2mb*dgnsiz.lt.10000.00 ) then
        write(*,'(a)') '                                  --------'
       else
        write(*,'(a)') '                                 ---------'
      endif
      write(*,203) '    Diagnostics subtotal (MB) = ',dim2mb*dgnsiz
      if ( dim2mb*(dynsiz+dgnsiz).lt.10000.00 ) then
        write(*,'(a)') '                                  --------'
        write(*,203) '  Total common blk  size (MB) = ',
     &               dim2mb*(dynsiz+dgnsiz)
        write(*,'(a)') '                                  ========'
       else
        write(*,'(a)') '                                 ---------'
        write(*,203) '  Total common blk  size (MB) = ',
     &               dim2mb*(dynsiz+dgnsiz)
        write(*,'(a)') '                                 ========='
      endif

      if ( stacks ) then
*     Compute size of local arrays/vectors
*     ====================================
*     Intended to be in order of first execution, to
*     investigate possible stack memory allocation problems
*     (N.B. former "small_local" option now removed)
      print *,' '
      print *,' Stack allocation of local arrays/vectors:'
      print *,' -----------------------------------------'
      print *,' Setup phase:'

*     constr
*     ------
      sizmin = 4
#  ifndef ocean_only
      sizmin = sizmin + nxpa*nypa*nla + 4*nla
#  endif
#  ifndef atmos_only
      sizmin = sizmin + 4*nlo
      sizarr = sizmin
      sizarr = sizarr + nxpo*nypo*nlo
#  endif
      write(*,203) '  constr     local memory   (MB) = ',dim2mb*sizarr
#  ifndef atmos_only
      write(*,203) '  constr     loc. with SAVE (MB) = ',dim2mb*sizmin
#  endif

*     xforc
*     -----
      sizmin = 8
#  ifndef ocean_only
      sizmin = sizmin + 5 + 34 + 1 + ndxr+1
#    ifdef tau_udiff
      sizmin = sizmin + 2
#    endif
      sizarr = sizmin + 2*nxpa*nypa + 4*nxpaor*nypaor + nxtaor*nytaor
     &        + nxto*nyto
#  endif
      write(*,203) '  xforc      local memory   (MB) = ',dim2mb*sizarr
#  ifndef ocean_only
      write(*,203) '  xforc      loc. with SAVE (MB) = ',dim2mb*sizmin
#  endif

*     bilint
*     ------
      sizarr = 7 + 3*nxpo
*     (scalars + 2 x integer + 2 x double)
      write(*,203) '  bilint     local memory   (MB) = ',dim2mb*sizarr

*     auvbli
*     ------
***   sizarr = 5 + 3*nxpaor
***   write(*,203) '  auvbli     local memory   (MB) = ',dim2mb*sizarr

*     auvbcu
*     ------
      sizarr = 5 + 4*16 + 6
      write(*,203) '  auvbcu     local memory   (MB) = ',dim2mb*sizarr

*     local common block bcudata
*     --------------------------
      sizarr = 5*16*16 + 16*(ndxr+1)*(ndxr+1)
      write(*,203) '  bcudata    local memory   (MB) = ',dim2mb*sizarr

*     bcuini
*     ------
      sizarr = 4*4*4*2*2 + 2 + 4
      write(*,203) '  bcuini     local memory   (MB) = ',dim2mb*sizarr

*     wts2bb
*     ------
      sizarr = 2*16*16 + 1 + 5
      write(*,203) '  wts2bb     local memory   (MB) = ',dim2mb*sizarr

*     homsol
*     ------
      sizmin = 7
#  ifndef ocean_only
      sizmin = sizmin + 2*nxpa*nypa + nxta
#  endif
#  ifndef atmos_only
      sizmin = sizmin + nxto
      sizarr = sizmin
#    ifdef cyclic_ocean
      sizarr = sizarr + 2*nxpo*nypo
#    endif
#  endif
      write(*,203) '  homsol     local memory   (MB) = ',dim2mb*sizarr
#  ifndef atmos_only
#    ifdef cyclic_ocean
      write(*,203) '  homsol     loc. with SAVE (MB) = ',dim2mb*sizmin
#    endif
#  endif

*     hsbxoc/hscyoc
*     -------------
#  ifndef atmos_only
#    ifdef cyclic_ocean
      sizarr = lwftoc + 2*(nypo-1) + 1
      write(*,203) '  hscyoc     local memory   (MB) = ',dim2mb*sizarr
#    else
      sizarr = lwftoc + 2*(nypo-1) + 1
      write(*,203) '  hsbxoc     local memory   (MB) = ',dim2mb*sizarr
#    endif
#  endif

*     monnc_comp
*     ----------
      sizmin = 21
#  ifndef ocean_only
      sizmin = sizmin + 3*nxta*nypa + 3*nxpa*nypa + nxpa*nypa
     &        + nxpa*nyta + nxta*nypa + nyta + nla
#  endif
      sizarr = sizmin
#  ifndef atmos_only
      sizmin = sizmin + nyto + nlo + 2
      sizarr = sizmin
      sizarr = sizarr + 4*nxto*nypo + 4*nxpo*nypo
     &        + nxpo*nypo + nxpo*nyto + nxto*nypo
#  endif
      write(*,203) '  monnc_comp local memory   (MB) = ',dim2mb*sizarr
#  ifndef atmos_only
      write(*,203) '  monnc_comp loc. with SAVE (MB) = ',dim2mb*sizmin
#  endif

#  ifdef use_netcdf
#    ifndef atmos_only
*     ocnc_init
      sizarr = nxpo + nypo + nlo
      write(*,203) '  ocnc_init  local arrays   (MB) = ',dim2mb*sizarr
*     ocnc_out
      sizarr = nxpo*nypo*nlo
      write(*,203) '  ocnc_out   local arrays   (MB) = ',dim2mb*sizarr
      sizarr = 0
      write(*,203) '  ocnc_out   loc. with SAVE (MB) = ',dim2mb*sizarr
#  ifdef qoc_diag
*     qocdiag_init
      sizarr = nxpo + nypo + nlo
      write(*,203) '  qocdiag_init loc. arrays (MB) = ',dim2mb*sizarr
*     qocdiag_out
      sizarr = 8*nxpo*nypo + nlo + 8
      write(*,203) '  qocdiag_out local arrays (MB) = ',dim2mb*sizarr
      sizarr = nlo + 8
      write(*,203) '  qocdiag_out   with SAVE  (MB) = ',dim2mb*sizarr
#  endif /* qoc_diag */
#    endif
#    ifndef ocean_only
*     atnc_init
      sizarr = nxpa + nypa + nla
      write(*,203) '  atnc_init local arrays   (MB) = ',dim2mb*sizarr
*     atnc_out
      sizarr = nxpa*nypa*nla
      write(*,203) '  atnc_out  local arrays   (MB) = ',dim2mb*sizarr
#    endif
#  endif
      print *,' Timestepping phase:'
#  ifndef atmos_only
*     oml
      sizarr = 2*nxto*nyto
      write(*,203) '  oml        local arrays   (MB) = ',dim2mb*sizarr
      sizarr = 0
      write(*,203) '  oml        loc. with SAVE (MB) = ',dim2mb*sizarr
*     omladf
      sizarr = (nxto+2)*nyto
      write(*,203) '  omladf     local arrays   (MB) = ',dim2mb*sizarr
      sizarr = 0
      write(*,203) '  omladf     loc. with SAVE (MB) = ',dim2mb*sizarr
*     qgostep
      sizarr = 4 + 2*nlo + nxpo*nypo + nxpo*nypo*nlo
      write(*,203) '  qgostep    local arrays   (MB) = ',dim2mb*sizarr
*     ocadif
      sizarr = 4 + nxpo*nypo
#  ifdef cyclic_ocean
      sizarr = sizarr + 8
#  endif
      write(*,203) '  ocadif     local arrays   (MB) = ',dim2mb*sizarr
*     ocinvq
      sizarr = 3*nlo + nxpo*nypo*nlo + nxto
#    ifdef cyclic_ocean
      sizarr = sizarr + 9*nlo + 2*(nlo-1)
#    else
      sizarr = sizarr + 3*nlo + 3*(nlo-1)
#    endif
      write(*,203) '  ocinvq     local arrays   (MB) = ',dim2mb*sizarr
      sizarr = sizarr - nxpo*nypo*nlo
      write(*,203) '  ocinvq     loc. with SAVE (MB) = ',dim2mb*sizarr
#  endif
#  ifndef ocean_only
*     aml
      sizarr = 2*(nla-1) + 3*nxta*nyta
      write(*,203) '  aml        local arrays   (MB) = ',dim2mb*sizarr
*     amladf
      sizarr = (nxta+2)*nyta
      write(*,203) '  amladf     local arrays   (MB) = ',dim2mb*sizarr
*     qgastep
      sizarr = 2*nla + nxpa*nypa + nxpa*nypa*nla
      write(*,203) '  qgastep    local arrays   (MB) = ',dim2mb*sizarr
*     atadif
      sizarr = 4*nxpa*nypa
      write(*,203) '  atadif     local arrays   (MB) = ',dim2mb*sizarr
*     atinvq
      sizarr = 14*nla + nxpa*nypa*nla + nxta + 2*(nla-1)
      write(*,203) '  atinvq     local arrays   (MB) = ',dim2mb*sizarr
*     tavatm
      sizarr = nxpa
      write(*,203) '  tavatm     local arrays   (MB) = ',dim2mb*sizarr
#  endif
#  ifndef atmos_only
      sizarr = 3*nxpo
      write(*,203) '  tavocn     local arrays   (MB) = ',dim2mb*sizarr
#  endif
      print *,' Finishing phase:'

*     tavout
*     ------
      sizmin = 0.0d0
#  ifndef ocean_only
      sizmin = sizmin + 1 + nxpa*nyta + nxta*nypa
#  endif
#  ifdef use_netcdf
#    ifndef ocean_only
      sizmin = sizmin + nxpa + nypa + nla
#    endif
#    ifndef atmos_only
      sizmin = sizmin + nxpo + nypo + nlo
#    endif
#  endif
      sizarr = sizmin
#  ifndef atmos_only
      sizmin = sizmin + 1
      sizarr = sizmin
      sizarr = sizarr + nxpo*nyto + nxto*nypo
#  endif
      write(*,203) '  tavout     local memory   (MB) = ',dim2mb*sizarr
#  ifndef atmos_only
      write(*,203) '  tavout     loc. with SAVE (MB) = ',dim2mb*sizmin
#  endif
      endif

  201 format(a,i12)
*     Integer format should be large enough to print max. I*4 number
  203 format(a,f10.3)

      END SUBROUTINE memreq
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
c
      SUBROUTINE resave (iounit)

*     Write an unformatted restart dump, either at end
*     of run, or periodically in case of program crash.
*     Dump contains time in years, and all progostic fields

*     Modules
      USE parameters
#ifndef ocean_only
      USE atstate
#endif
#ifndef atmos_only
      USE ocstate
#endif
      USE intrfac, ONLY : sst, sstm, ast, astm, hmixa, hmixam
      USE timinfo, ONLY : tyrs

      IMPLICIT NONE

*     Subroutine arguments
      integer iounit
*
*     Dump restart file on unit iounit
      write (iounit) tyrs
#ifndef atmos_only
      write (iounit) po,pom
#endif
#ifndef ocean_only
      write (iounit) pa,pam
#endif
      write (iounit) sst,sstm
      write (iounit) ast,astm
      write (iounit) hmixa,hmixam

      END SUBROUTINE resave
c
c***********************************************************************
