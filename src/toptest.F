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
*
*     Simplified version of q-gcm.F for testing of topsubs.F

*     Modules
      USE parameters
      USE atconst
      USE occonst
      USE intrfac, ONLY : hmoc, hmat, hmamin
      USE eigmode
      USE topsubs

      IMPLICIT NONE

#ifdef use_netcdf
      INCLUDE 'netcdf.inc'
#endif /* use_netcdf */
*
*     Local parameters
      DOUBLE PRECISION PI,TWOPI,PIBY2
      PARAMETER ( PI=3.14159265358979324D0,
     &            TWOPI=6.28318530717958648D0,
     &            PIBY2=1.57079632679489662D0 )
*
*     Local variables
      integer i,j,k,lenod,oceanonly,atmosonly
      character (len=80) :: outdir
*     Extra variables for topography
      character (len=80) :: topocname,topatname

*     Specify input parameters in included file

***   outdir = './outdata/'   !! Output directory
      outdir = './'           !! Output directory
*     (Write test results into current working directory)

      dxo = 10.0D3            !! Ocean grid spacing (km)

*     Mixed layer parameters
*     ----------------------
*     Ocean
      hmoc = 100.0d0          !! Fixed ocean  mixed layer depth  (m)
*     Atmos.
      hmat = 1000.0d0         !! Fixed atmos. mixed layer depth  (m)
      hmamin = 100.0d0        !! Min. allowed atmos. m. l. depth (m)

*     Oceanic QG layer parameters
*     ---------------------------
*     Reduced gravities, temperatures, friction coeffs and thicknesses

      hoc(1) = 300.d0         !! Height of ocean layer 1 (m)

      gpoc(1) = 0.05d0        !! Reduced gravity for ocean 1/2 interface (m s^-2)

      hoc(2) = 1100.d0        !! Height of ocean layer 2 (m)

      gpoc(2) = gpoc(1)/2.0d0 !! Reduced gravity for ocean 2/3 interface (m s^-2)

      hoc(3) = 2600.d0        !! Height of ocean layer 3 (m)

*     Atmospheric QG layer parameters
*     -------------------------------
*     Reduced gravities, temperatures, friction coeffs and thicknesses

      hat(1) = 2000.d0        !! Height of atmos. layer 1 (m)

      gpat(1) = 1.2d0         !! Reduced gravity for atmos. 1/2 interface (m s^-2)

      hat(2) = 3000.d0        !! Height of atmos. layer 2 (m)

      gpat(2) = 0.4d0         !! Reduced gravity for atmos. 2/3 interface (m s^-2)

      hat(3) = 4000.d0        !! Height of atmos. layer 3 (m)

      topocname = 'define'    !! Ocean topography file. Possible values:
                              !! 'zero', 'define', or <restart_filename>
      topatname = 'define'    !! Atmos topography file. Possible values:
                              !! 'zero', 'define', or <restart_filename>

      print *,' '
      lenod = index(outdir, '   ') - 1
      print *,' outdir = ',outdir(1:lenod)

*     Specify configuration being used
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

      print *,' '
      print *,' Control parameters:'
      print *,' -------------------'
      if ( oceanonly.eq.1 .and. atmosonly.eq.1 ) then
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
*     Check some of the grid parameters
      if ( nla.lt.2 .or. nlo.lt.2 ) then
        print *,' nla, nlo = ',nla,nlo
        print *,' Inadequate nla or nlo, needs to be at least 2'
        print *,' Program terminates'
        stop
      endif
      if ( nxta.lt.nxaooc ) then
        print *,' nxta, nxaooc = ',nxta,nxaooc
        print *,' Inadequate nxta needs to be at least nxaooc'
        print *,' Program terminates'
        stop
      endif
      if ( nyta.lt.nyaooc ) then
        print *,' nyta, nyaooc = ',nyta,nyaooc
        print *,' Inadequate nyta needs to be at least nyaooc'
        print *,' Program terminates'
        stop
      endif

*     Derive atmos gridspace and ocean timestep
*     Derive larger from smaller to ensure integer ratio
      dxa = ndxr*dxo
      write(*,201) '  Atmos/ocean grid ratio ndxr = ',ndxr
      write(*,201) '  Atmos. gridcells over ocean = ',nxaooc,nyaooc
      write(*,201) '  Ocn start indices  nx1, ny1 = ',nx1,ny1
      write(*,214) '  Coriolis param.  (rad s^-1) = ',fnot
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
      write(*,203) '  Mixed layer thickness   (m) = ',hmoc
      write(*,203) '  Layer thicknesses hoc   (m) = ',(hoc(k),k=1,nlo)
      write(*,207) '  Reduced gravities  (m s^-2) = ',
     &             (gpoc(k),k=1,nlo-1)
      write(*,206) '  Baroclinic wavespeeds (m/s) = ',
     &             (cphsoc(k),k=2,nlo)
      write(*,204) '  Deformation radii      (km) = ',
     &             (1.0d-3*rdefoc(k),k=2,nlo)
      write(*,205) '  Gridlengths per radius      = ',
     &             (rdefoc(k)/dxo,k=2,nlo)

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
      write(*,205) '  Rossby number   Beta*yla/f0 = ',beta*yla/abs(fnot)
      write(*,214) '  f range S -> N   (rad s^-1) = ',
     &             fnot+beta*yparel(1),fnot+beta*yparel(nypa)
      write(*,214) '  Midlatitude Coriolis param  = ',
     &             fnot+beta*0.5d0*( yparel(1) + yparel(nypa) )
      write(*,203) '  Mixed layer thickness   (m) = ',hmat
      write(*,203) '  Min. mixed layer thick. (m) = ',hmamin
      write(*,203) '  Layer thicknesses hat   (m) = ',(hat(k),k=1,nla)
      write(*,207) '  Reduced gravities  (m s^-2) = ',
     &             (gpat(k),k=1,nla-1)
      write(*,206) '  Baroclinic wavespeeds (m/s) = ',
     &             (cphsat(k),k=2,nla)
      write(*,204) '  Deformation radii      (km) = ',
     &             (1.0d-3*rdefat(k),k=2,nla)
      write(*,205) '  Gridlengths per radius      = ',
     &             (rdefat(k)/dxa,k=2,nla)

  201 format(a,9i13)
  203 format(a,9f13.3)
  204 format(a,9f13.4)
  205 format(a,9f13.5)
  206 format(a,9f13.6)
  207 format(a,9f13.7)
  213 format(a,1p,9d13.3)
  214 format(a,1p,9d13.4)

*     Topography section
*     ==================
      call topset (topocname, topatname, outdir, lenod, .true.)

      stop
      end

c***********************************************************************
