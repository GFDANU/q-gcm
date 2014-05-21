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
      MODULE topsubs

*     Contains subroutines topset and topout_nc, which set and
*     output in netCDF format the atmospheric and oceanic topographies.

*     Modules

      IMPLICIT NONE

      PRIVATE

      PUBLIC  :: topset
      PRIVATE :: topout_nc

      CONTAINS

c***********************************************************************
c
      SUBROUTINE topset (topocname, topatname, outdir, lenod, atmver)
*
*     Sets the topography in the atmosphere and ocean, as specified by
*     the input arguments topocname, topatname, which are character
*     strings either specifying a file to read existing topography from,
*     or determining the initialisation choice within this subroutine
*     or the treatment of an already existing field. Valid values are:
*     'flat', 'define', 'extant', or the name of a file to be read.
*     The input file will be in netCDF format (written by subroutine
*     topout_nc) if the use_netcdf option is activated, otherwise it is
*     assumed to be an unformatted binary dump of appropriate endian-ness.
*     outdir (with significant length lenod) is the name of the directory
*     in which the netCDF output file (called topog.nc) will be written.
*     The logical switch atmver controls whether consistency of atmosphere
*     and ocean topographies is verified (i.e. whether the program checks
*     that there is no atmospheric topography over ocean gridpoints.
*     This switch should be set to .true. when running the model normally,
*     but to .false. if necessary when preparing the topography dataset(s).

*     Modules
      USE parameters
      USE atconst
      USE occonst
      USE intsubs, ONLY : xintp
      USE nc_subs, ONLY : handle_err

      IMPLICIT NONE

#ifdef use_netcdf
      INCLUDE 'netcdf.inc'
#endif /* use_netcdf */
*
*     Subroutine arguments
      character (len=*) :: topocname,topatname,outdir
      integer lenod
      logical atmver
*
*     Local parameters
      character (len=*), parameter :: subnam = 'topset'
*
*     Local variables
      integer i,j
      double precision dxlo,dxhi,dcent,dhwid,xrel,davg,ddymin,ddymax
      double precision xacent
      integer itlo,ithi,jsec,jtlo,jthi
#ifdef use_netcdf
*     Extra variables for netCDF input/output
      integer ncstat,tempid,varid,dim,idim,jdim
#else
      logical fexist
#endif /* use_netcdf */

      print *,' '
      print *,' Topography:'
      print *,' -----------'

*     Initialise ocean topography
*     ===========================
      if ( trim(adjustl(topocname)).eq.'flat' ) then
*       Flat-bottom case
*       ----------------
*       Set topography to zero
        do j=1,nypo
          do i=1,nxpo
            dtopoc(i,j) = 0.0d0
          enddo
        enddo
        print *,' Ocean topography set to zero (flat)'
       else if ( trim(adjustl(topocname)).eq.'define' ) then
*       Define topography within this program
*       -------------------------------------
*       Default topography to zero
        do j=1,nypo
          do i=1,nxpo
            dtopoc(i,j) = 0.0d0
          enddo
        enddo
*       Optionally add topographic features
*       Define w.r.t. ocean basin coordinates
!!      Put in a mid-Atlantic ridge for testing:
        dxlo = 2000.0d3
        dxhi = 2600.0d3
        itlo = 1 + nint( dxlo/dxo )
        ithi = 1 + nint( dxhi/dxo )
        dcent = 0.5d0*( dxlo + dxhi )
        dhwid = 0.5d0*( dxhi - dxlo )
        do j=1,nypo
          do i=1,nxpo
            xrel = xpo(i) - xpo(1)
            dtopoc(i,j) = 1000.0d0*( 1.0d0 - abs(xrel-dcent)/dhwid )
            dtopoc(i,j) = max( 0.0d0, dtopoc(i,j) )
          enddo
        enddo
        print *,' Ocean topography defined in program'
*       Print zonal section of topography
        jsec = ( 1 + nypo )/2
        jtlo = max(  1 ,jsec-1)
        jthi = min(nxpo,jsec+1)
        write(*,'(a,i4)') '  Zonal ocean section around j = ',jsec
        write(*,'(a,i9,4i11)') '     i    xrel(km)',(j,j=jtlo,jthi)
        do i=max(1,itlo-3),min(nxpo,ithi+3)
          xrel = 1.0d-3*( xpo(i) - xpo(1) )
          write(*,'(i6,f11.2,5f11.3)') i,xrel,(dtopoc(i,j),j=jtlo,jthi)
        enddo
       else if ( trim(adjustl(topocname)).eq.'extant' ) then
*       Case of pre-existing topography
*       -------------------------------
        print *,' Ocean topography already extant'
*       In either case, print zonal section of topography
        jsec = ( 1 + nypo )/2
        jtlo = max(  1 ,jsec-3)
        jthi = min(nxpo,jsec+3)
        write(*,'(a,i4)') '  Zonal ocean section around j = ',jsec
        write(*,'(a,i8,6i10)') '      i   xrel(km)',(j,j=jtlo,jthi)
        do i=1,min(nxpo,25)
          xrel = 1.0d-3*( xpo(i) - xpo(1) )
          write(*,'(i7,f10.2,7f10.2)') i,xrel,(dtopoc(i,j),j=jtlo,jthi)
        enddo
        print *,' '
        do i=max(26,nxpo-24),nxpo
          xrel = 1.0d-3*( xpo(i) - xpo(1) )
          write(*,'(i7,f10.2,7f10.2)') i,xrel,(dtopoc(i,j),j=jtlo,jthi)
        enddo
       else
*       Case of topography to be read from a file
*       -----------------------------------------
#ifdef use_netcdf
*       Read ocean topography from an existing netCDF file
        ncstat = nf_open(trim(adjustl(topocname)), NF_NOWRITE, tempid)
        if ( ncstat.ne.NF_NOERR ) call handle_err (ncstat, subnam)
*       First check that the dimensions are correct
        ncstat = nf_inq_dimid(tempid, 'xpo', dim)
        if ( ncstat.ne.NF_NOERR ) call handle_err (ncstat, subnam)
        ncstat = nf_inq_dimlen(tempid, dim, idim)
        if ( ncstat.ne.NF_NOERR ) call handle_err (ncstat, subnam)
        if ( idim.ne.nxpo ) then
          print *,' ocean topography netCDF file error'
          print *,' xpo dimension not equal to nxpo:'
          print *,' dimensions are: ',idim,nxpo
          print *,' Program terminates in topset'
          stop
        endif
        ncstat = nf_inq_dimid(tempid, 'ypo', dim)
        if ( ncstat.ne.NF_NOERR ) call handle_err (ncstat, subnam)
        ncstat = nf_inq_dimlen(tempid, dim, jdim)
        if ( ncstat.ne.NF_NOERR ) call handle_err (ncstat, subnam)
        if ( jdim.ne.nypo ) then
          print *,' ocean topography netCDF file error'
          print *,' ypo dimension not equal to nypo:'
          print *,' dimensions are: ',jdim,nypo
          print *,' Program terminates in topset'
          stop
        endif
*       Now read ocean topography data
        ncstat = nf_inq_varid(tempid, 'dtopoc', varid)
        if ( ncstat.ne.NF_NOERR ) call handle_err (ncstat, subnam)
        ncstat = nf_get_var_double(tempid, varid, dtopoc)
        if ( ncstat.ne.NF_NOERR ) call handle_err (ncstat, subnam)
        ncstat = nf_close(tempid)
        if ( ncstat.ne.NF_NOERR ) call handle_err (ncstat, subnam)
        print *,' Ocean topography read in from netCDF file'
        write(*,'(a,a)') '  Ocean topography filename = ',
     &                   trim(adjustl(topocname))
#else
*       Read ocean topography from an existing binary file
        inquire (file=trim(adjustl(topocname)), exist=fexist)
        if ( fexist ) then
          open (50, file=trim(adjustl(topocname)),
     &          form='unformatted', status='old')
          read (50) dtopoc
          close(50)
          print *,' Ocean topography read in from unformatted binary',
     &            ' file '
          write(*,'(a,a)') '  Ocean topography filename = ',
     &                     trim(adjustl(topocname))
         else
          print *,' Ocean topography unformatted binary file not found'
          write(*,'(a,a)') '  Filename = ',trim(adjustl(topocname))
          print *,' Program terminates in topset'
          stop
        endif
#endif /* use_netcdf */
      endif

*     If cyclic and not flat-bottom, check cyclicity of topography
*     ------------------------------------------------------------
#ifdef cyclic_ocean
      if ( trim(adjustl(topocname)).ne.'flat' ) then
        do j=1,nypo
          if ( dtopoc(1,j).ne.dtopoc(nxpo,j) ) then
            print *,' *** WARNING *** problem with specified topography'
            print *,' Ocean topography not exactly cyclic for j = ',j
            print *,' dtopoc values are: ',dtopoc(1,j),dtopoc(nxpo,j)
          endif
        enddo
      endif
#endif

*     Compute average of ocean topography
*     -----------------------------------
      call xintp (davg, dtopoc, nxpo, nypo)
      davgoc = davg*ocnorm
      write(*,204) '  Avge of ocean  topogr.  (m) = ',davgoc

*     Initialise atmosphere topography
*     ================================
      print *,' '
      if ( trim(adjustl(topatname)).eq.'flat' ) then
*       Flat-bottom case
*       ----------------
*       Set topography to zero
        do j=1,nypa
          do i=1,nxpa
            dtopat(i,j) = 0.0d0
          enddo
        enddo
        print *,' Atmosphere topography set to zero (flat)'
       else if ( trim(adjustl(topatname)).eq.'define' ) then
*       Define topography within this program
*       -------------------------------------
*       Default topography to zero
        do j=1,nypa
          do i=1,nxpa
            dtopat(i,j) = 0.0d0
          enddo
        enddo
*       Optionally add topographic features
!!      Put in Rocky Mountains for testing:
!!      Original parameters from Hogg et al 2006 paper
!!      in JClimate, for 15360 km length atmosphere.
        dxlo = 1800.0d3
        dxhi = 5400.0d3
        dcent = 0.5d0*( dxlo + dxhi )
        dhwid = 0.5d0*( dxhi - dxlo )
!!      Revised parameters for extended 30720 km length domain
!!      1st attempt: same width & height as in
!!      Hogg et al 2006, but centred at x = 9000.0
        dxlo =  7200.0d3
        dxhi = 10800.0d3
        dcent = 0.5d0*( dxlo + dxhi )
        dhwid = 0.5d0*( dxhi - dxlo )
!!      2nd attempt: centred at x = 8800.0 km,
!!      slightly narrowed, and with the central axis
!!      sloping by +/- 2000.0 km across the domain
        dcent = 8800.0d3
        dhwid = 1440.0d3
        dxlo = dcent - dhwid
        dxhi = dcent + dhwid
        itlo = 1 + nint( dxlo/dxa )
        ithi = 1 + nint( dxhi/dxa )
        do j=1,nypa
!!        xacent = dcent
          xacent = dcent - 2000.0d3*yparel(j)/(0.5d0*yla)
          do i=1,nxpa
            dtopat(i,j) = 1000.0d0*( 1.0d0 - abs(xpa(i)-xacent)/dhwid )
            dtopat(i,j) = max( 0.0d0, dtopat(i,j) )
          enddo
        enddo
        print *,' Atmosphere topography defined in program'
*       Print zonal section of topography
*       Format statements currently permit j-range of up to jsec +/- 2
        jsec = ( 1 + nypa )/2
        jtlo = max(  1 ,jsec-2)
        jthi = min(nxpa,jsec+2)
        write(*,'(a,i4)') '  Zonal atmosphere section around j = ',jsec
        write(*,'(a,i10,4i11)') '     i    xpa(km)',(j,j=jtlo,jthi)
        do i=max(1,itlo-3),min(nxpa,ithi+3)
          xrel = 1.0d-3*xpa(i)
          write(*,'(i6,f11.2,5f11.3)') i,xrel,(dtopat(i,j),j=jtlo,jthi)
        enddo
       else if ( trim(adjustl(topatname)).eq.'extant' ) then
*       Case of pre-existing topography
*       -------------------------------
        print *,' Atmosphere topography already extant'
*       In either case, print zonal section of topography
        jsec = ( 1 + nypa )/2
        jtlo = max(  1 ,jsec-3)
        jthi = min(nxpa,jsec+3)
        write(*,'(a,i4)') '  Zonal atmosphere section around j = ',jsec
        write(*,'(a,i8,6i10)') '      i    xpa(km)',(j,j=jtlo,jthi)
        do i=1,min(nxpa,25)
          xrel = 1.0d-3*xpa(i)
          write(*,'(i7,f10.2,7f10.2)') i,xrel,(dtopat(i,j),j=jtlo,jthi)
        enddo
        print *,' '
        do i=max(26,nxpa-24),nxpa
          xrel = 1.0d-3*xpa(i)
          write(*,'(i7,f10.2,7f10.2)') i,xrel,(dtopat(i,j),j=jtlo,jthi)
        enddo
       else
*       Case of topography to be read from a file
*       -----------------------------------------
#ifdef use_netcdf
*       Read atmosphere topography from an existing netCDF file
        ncstat = nf_open(trim(adjustl(topatname)), NF_NOWRITE, tempid)
        if ( ncstat.ne.NF_NOERR ) call handle_err (ncstat, subnam)
*       First check that the dimensions are correct
        ncstat = nf_inq_dimid(tempid, 'xpa', dim)
        if ( ncstat.ne.NF_NOERR ) call handle_err (ncstat, subnam)
        ncstat = nf_inq_dimlen(tempid, dim, idim)
        if ( ncstat.ne.NF_NOERR ) call handle_err (ncstat, subnam)
        if ( idim.ne.nxpa ) then
          print *,' atmos. topography netCDF file error'
          print *,' xpa dimension not equal to nxpa:'
          print *,' dimensions are: ',idim,nxpa
          print *,' Program terminates in topset'
          stop
        endif
        ncstat = nf_inq_dimid(tempid, 'ypa', dim)
        if ( ncstat.ne.NF_NOERR ) call handle_err (ncstat, subnam)
        ncstat = nf_inq_dimlen(tempid, dim, jdim)
        if ( ncstat.ne.NF_NOERR ) call handle_err (ncstat, subnam)
        if ( jdim.ne.nypa ) then
          print *,' atmos. topography netCDF file error'
          print *,' ypa dimension not equal to nypa:'
          print *,' dimensions are: ',jdim,nypa
          print *,' Program terminates in topset'
          stop
        endif
*       Now read atmosphere topography data
        ncstat = nf_inq_varid(tempid, 'dtopat', varid)
        if ( ncstat.ne.NF_NOERR ) call handle_err (ncstat, subnam)
        ncstat = nf_get_var_double(tempid, varid, dtopat)
        if ( ncstat.ne.NF_NOERR ) call handle_err (ncstat, subnam)
        ncstat = nf_close(tempid)
        if ( ncstat.ne.NF_NOERR ) call handle_err (ncstat, subnam)
        print *,' Atmos. topography read in from netCDF file'
        write(*,'(a,a)') '  Atmos. topography filename = ',
     &                   trim(adjustl(topatname))
#else
*       Read atmosphere topography from an existing binary file
        inquire (file=trim(adjustl(topatname)), exist=fexist)
        if ( fexist ) then
          open (50, file=trim(adjustl(topatname)),
     &          form='unformatted', status='old')
          read (50) dtopat
          close(50)
          print *,' Atmos. topography read in from unformatted binary',
     &            ' file'
          write(*,'(a,a)') '  Atmos. topography filename = ',
     &                     trim(adjustl(topatname))
         else
          print *,' Atmos. topography unformatted binary file not found'
          write(*,'(a,a)') '  Filename = ',trim(adjustl(topatname))
          print *,' Program terminates in topset'
          stop
        endif
#endif /* use_netcdf */
      endif

*     Check there is no atmosphere topography over
*     the ocean (including ocean boundary points)
*     --------------------------------------------
*     This check may need to be switched off during
*     the phase of producing suitable topography,
*     but should be turned on during normal running
      if ( atmver ) then
        do j=ny1,ny1+nyaooc
          do i=nx1,nx1+nxaooc
            if ( dtopat(i,j).ne.0.0d0 ) then
              print *,' Nonzero atmosphere topography over ocean'
              print *,' Problem occurs at i, j, xpa(km), ypa(km) = ',
     &                i,j,1.0d-3*xpa(i),1.0d-3*ypa(j)
              write(*,'(a,2f11.2)') '  Ocean limits in x (km) are: ',
     &                1.0d-3*xpa(nx1),1.0d-3*xpa(nx1+nxaooc)
              write(*,'(a,2f11.2)') '  Ocean limits in y (km) are: ',
     &                1.0d-3*ypa(ny1),1.0d-3*ypa(ny1+nyaooc)
              print *,' Program terminates in topset'
              stop
            endif
          enddo
        enddo
      endif

*     If not flat-bottom, check cyclicity of topography
*     -------------------------------------------------
      if ( trim(adjustl(topatname)).ne.'flat' ) then
        do j=1,nypa
          if ( dtopat(1,j).ne.dtopat(nxpa,j) ) then
            print *,' *** WARNING *** problem with specified topography'
            print *,' Atmos. topography not exactly cyclic for j = ',j
            print *,' dtopat values are: ',dtopat(1,j),dtopat(nxpa,j)
          endif
        enddo
      endif

*     Compute average of atmosphere topography
*     ----------------------------------------
      call xintp (davg, dtopat, nxpa, nypa)
      davgat = davg*atnorm
      write(*,204) '  Avge of atmos. topogr.  (m) = ',davgat

*     Write out physical topography
*     =============================
#ifdef use_netcdf
      ncstat = nf_create (outdir(1:lenod)//'/topog.nc',
     &                    NF_CLOBBER, tempid)
      if ( ncstat.ne.NF_NOERR ) call handle_err (ncstat, subnam)
      call topout_nc(tempid)
      ncstat = nf_close(tempid)
      if ( ncstat.ne.NF_NOERR ) call handle_err (ncstat, subnam)
      print *,' netCDF copy of topography written'
#endif /* use_netcdf */

*     Rescale physical topography to get dynamic topography
*     =====================================================
      print *,' '
*     Ocean
*     -----
      ddymin =  1.0d30
      ddymax = -1.0d30
      do j=1,nypo
        do i=1,nxpo
          ddynoc(i,j) = (fnot/hoc(nlo))*dtopoc(i,j)
          ddymin = min( ddymin, ddynoc(i,j) )
          ddymax = max( ddymax, ddynoc(i,j) )
        enddo
      enddo
      write(*,214) '  Min/max frac. ocean  topog. = ',
     &             ddymin/fnot,ddymax/fnot
*     Atmosphere
*     ----------
      ddymin =  1.0d30
      ddymax = -1.0d30
      do j=1,nypa
        do i=1,nxpa
          ddynat(i,j) = (fnot/hat(1))*dtopat(i,j)
          ddymin = min( ddymin, ddynat(i,j) )
          ddymax = max( ddymax, ddynat(i,j) )
        enddo
      enddo
      write(*,214) '  Min/max frac. atmos. topog. = ',
     &             ddymin/fnot,ddymax/fnot

  204 format(a,5f13.4)
  214 format(a,1p,5d13.4)

      END SUBROUTINE topset

c***********************************************************************

      SUBROUTINE topout_nc (tempid)

*     Modules
      USE parameters

#ifdef use_netcdf
      USE atconst
      USE occonst
#endif /* use_netcdf */
      USE nc_subs, ONLY : handle_err

      IMPLICIT NONE

#ifdef use_netcdf
      INCLUDE 'netcdf.inc'
#endif /* use_netcdf */
*
*     Subroutine arguments
      integer tempid

#ifdef use_netcdf

*     Local parameters
      character (len=*), parameter :: subnam = 'topout_nc'
*
*     Local variables
      integer i,j,ncstat,xopdim,yopdim,xapdim,yapdim,xdims
      integer xop_id,yop_id,xap_id,yap_id, dtopoc_id,dtopat_id
      integer dims(2),start(2),count(2)
      double precision xxo(nxpo),yyo(nypo),xxa(nxpa),yya(nypa)

!! Initialise netcdf file
      ncstat = nf_def_dim(tempid, 'xpo', nxpo, xopdim)
      if ( ncstat.ne.NF_NOERR ) call handle_err (ncstat, subnam)
      ncstat = nf_def_dim(tempid, 'ypo', nypo, yopdim)
      if ( ncstat.ne.NF_NOERR ) call handle_err (ncstat, subnam)

      ncstat = nf_def_dim(tempid, 'xpa', nxpa, xapdim)
      if ( ncstat.ne.NF_NOERR ) call handle_err (ncstat, subnam)
      ncstat = nf_def_dim(tempid, 'ypa', nypa, yapdim)
      if ( ncstat.ne.NF_NOERR ) call handle_err (ncstat, subnam)

!! Define a one-dimensional variable called 'xpo'
!! which stores x gridpoints (p-grid)
      xdims = xopdim
      ncstat = nf_def_var(tempid, 'xpo', NF_DOUBLE, 1, xdims, xop_id)
      if ( ncstat.ne.NF_NOERR ) call handle_err (ncstat, subnam)
      ncstat = nf_put_att_text(tempid,xop_id,'units',2,'km')
      if ( ncstat.ne.NF_NOERR ) call handle_err (ncstat, subnam)

!! Define a one-dimensional variable called 'ypo'
!! which stores y gridpoints (p-grid)
      xdims = yopdim
      ncstat = nf_def_var(tempid, 'ypo', NF_DOUBLE, 1, xdims, yop_id)
      if ( ncstat.ne.NF_NOERR ) call handle_err (ncstat, subnam)
      ncstat = nf_put_att_text(tempid,yop_id,'units',2,'km')
      if ( ncstat.ne.NF_NOERR ) call handle_err (ncstat, subnam)

!! Define a one-dimensional variable called 'xpa'
!! which stores x gridpoints (p-grid)
      xdims = xapdim
      ncstat = nf_def_var(tempid, 'xpa', NF_DOUBLE, 1, xdims, xap_id)
      if ( ncstat.ne.NF_NOERR ) call handle_err (ncstat, subnam)
      ncstat = nf_put_att_text(tempid,xap_id,'units',2,'km')
      if ( ncstat.ne.NF_NOERR ) call handle_err (ncstat, subnam)

!! Define a one-dimensional variable called 'ypa'
!! which stores y gridpoints (p-grid)
      xdims = yapdim
      ncstat = nf_def_var(tempid, 'ypa', NF_DOUBLE, 1, xdims, yap_id)
      if ( ncstat.ne.NF_NOERR ) call handle_err (ncstat, subnam)
      ncstat = nf_put_att_text(tempid,yap_id,'units',2,'km')
      if ( ncstat.ne.NF_NOERR ) call handle_err (ncstat, subnam)

!! Define a 2-dimensional variable called 'dtopoc'
!! to store ocean topography
      dims(1) = xopdim
      dims(2) = yopdim
      ncstat = nf_def_var(tempid, 'dtopoc', NF_DOUBLE,
     &                    2, dims, dtopoc_id)
      if ( ncstat.ne.NF_NOERR ) call handle_err (ncstat, subnam)
      ncstat = nf_put_att_text(tempid,dtopoc_id,'units',1,'m')
      if ( ncstat.ne.NF_NOERR ) call handle_err (ncstat, subnam)

!! Define a 2-dimensional variable called 'dtopat'
!! to store atmosphere topography
      dims(1) = xapdim
      dims(2) = yapdim
      ncstat = nf_def_var(tempid, 'dtopat', NF_DOUBLE,
     &                    2, dims, dtopat_id)
      if ( ncstat.ne.NF_NOERR ) call handle_err (ncstat, subnam)
      ncstat = nf_put_att_text(tempid,dtopat_id,'units',1,'m')
      if ( ncstat.ne.NF_NOERR ) call handle_err (ncstat, subnam)

!! Leave definition mode: entering data mode.
      ncstat = nf_enddef(tempid)
      if ( ncstat.ne.NF_NOERR ) call handle_err (ncstat, subnam)

!! OCEAN AXES
!! Calculate x gridpoints (in km) and store in 'x' array.
      do i=1,nxpo
        xxo(i) = 1.0d-3*( xpo(i) - xpo(1) )
      enddo
      ncstat = nf_put_vara_double(tempid, xop_id, 1, nxpo, xxo)
      if ( ncstat.ne.NF_NOERR ) call handle_err (ncstat, subnam)
!! Calculate y gridpoints (in km) and store in 'y' array.
      do j=1,nypo
        yyo(j) = 1.0d-3*( ypo(j) - ypo(1) )
      enddo
      ncstat = nf_put_vara_double(tempid, yop_id, 1, nypo, yyo)
      if ( ncstat.ne.NF_NOERR ) call handle_err (ncstat, subnam)

!! ATMOSPHERE AXES
!! Calculate x gridpoints (in km) and store in 'x' array.
      do i=1,nxpa
        xxa(i) = 1.0d-3*xpa(i)
      enddo
      ncstat = nf_put_vara_double(tempid, xap_id, 1, nxpa, xxa)
      if ( ncstat.ne.NF_NOERR ) call handle_err (ncstat, subnam)
!! Calculate y gridpoints (in km) and store in 'y' array.
      do j=1,nypa
        yya(j) = 1.0d-3*ypa(j)
      enddo
      ncstat = nf_put_vara_double(tempid, yap_id, 1, nypa, yya)
      if ( ncstat.ne.NF_NOERR ) call handle_err (ncstat, subnam)

!! Now put ocean topography
      start(1) = 1
      start(2) = 1
      count(1) = nxpo
      count(2) = nypo
      ncstat = nf_put_vara_double(tempid, dtopoc_id,
     &                            start, count, dtopoc)
      if ( ncstat.ne.NF_NOERR ) call handle_err (ncstat, subnam)

!! Now put atmos topography
      start(1) = 1
      start(2) = 1
      count(1) = nxpa
      count(2) = nypa
      ncstat = nf_put_vara_double(tempid, dtopat_id,
     &                            start, count, dtopat)
      if ( ncstat.ne.NF_NOERR ) call handle_err (ncstat, subnam)

!! Done
#endif /* use_netcdf */

      END SUBROUTINE topout_nc

c***********************************************************************
c
      END MODULE topsubs
c
c***********************************************************************
