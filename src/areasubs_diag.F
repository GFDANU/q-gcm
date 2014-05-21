c***********************************************************************
c     Q-GCM Version 1.5.0 : last modified 27/06/2013
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
      MODULE areasubs

*     Contains subprograms for computing and outputting (in netCDF
*     format) the averages over various specified areas within
*     the domains of various atmospheric and oceanic fields.
*     The areas are specified by a file (or link) called 'areas.limits'.
*     Output interval is nocmon atmospheric timesteps = dgnday days
*     (controlled by when areavg is called from the main program)

      IMPLICIT NONE

      PRIVATE

      integer, PUBLIC, SAVE :: arencid

#ifdef get_areav

      PUBLIC  :: areavg
      PRIVATE :: areint

      CONTAINS
#endif /* get_areav */

c***********************************************************************
c
#ifdef get_areav
*
      SUBROUTINE areavg
*
*     Defines several subareas of ocean and atmosphere, and
*     computes averages of various quantities over these areas
*     Number and limits of areas specified in file "areas.limits"

*     Modules
      USE parameters
      USE nc_subs, ONLY : handle_err
      USE timinfo, ONLY : ntdone, nocmon, tyrs, numoutsteps
#ifdef get_areav
      USE atconst
      USE occonst
#endif
      USE intrfac, ONLY : sst, ast

      IMPLICIT NONE

#ifdef use_netcdf
      INCLUDE 'netcdf.inc'
#endif /* use_netcdf */
*
*     Subroutine arguments
*
#ifdef get_areav
*
*     Local parameters
      character (len=*), parameter :: subnam = 'areavg'
      integer namxoc,namxat,arunit
      parameter ( namxoc=9, namxat=9, arunit=47 )
*
*     Local arrays and variables
      integer m,nareoc,nareat,ia
*
*     Extra variables for limited subareas.
      double precision
     &        xlooc(namxoc),xhioc(namxoc),ylooc(namxoc),yhioc(namxoc),
     &        fwtoc(namxoc),fetoc(namxoc),fstoc(namxoc),fntoc(namxoc),
     &        fwpoc(namxoc),fepoc(namxoc),fspoc(namxoc),fnpoc(namxoc),
     &        tavoc(namxoc),
     &        xloat(namxat),xhiat(namxat),yloat(namxat),yhiat(namxat),
     &        fwtat(namxat),fetat(namxat),fstat(namxat),fntat(namxat),
     &        fwpat(namxat),fepat(namxat),fspat(namxat),fnpat(namxat),
     &        tavat(namxat),
     &        xoclim,yoclim,rilo,rihi,rjlo,rjhi,frlo,frhi
      integer i1toc(namxoc),i2toc(namxoc),j1toc(namxoc),j2toc(namxoc),
     &        i1poc(namxoc),i2poc(namxoc),j1poc(namxoc),j2poc(namxoc),
     &        i1tat(namxat),i2tat(namxat),j1tat(namxat),j2tat(namxat),
     &        i1pat(namxat),i2pat(namxat),j1pat(namxat),j2pat(namxat)
      character (len=3) :: areaat(namxat),areaoc(namxoc)
      logical first
*
      SAVE xlooc,xhioc,ylooc,yhioc,areaoc,nareoc,
     &     i1toc,i2toc,j1toc,j2toc,i1poc,i2poc,j1poc,j2poc,
     &     fwtoc,fetoc,fstoc,fntoc,fwpoc,fepoc,fspoc,fnpoc,
     &     xloat,xhiat,yloat,yhiat,areaat,nareat,
     &     i1tat,i2tat,j1tat,j2tat,i1pat,i2pat,j1pat,j2pat,
     &     fwtat,fetat,fstat,fntat,fwpat,fepat,fspat,fnpat,
     &     first

#  ifdef use_netcdf
      integer status, timedim, nocdim, natdim, edims(2)
      integer time_id, noc_id, nat_id, ocdat_id, atdat_id
      integer tomp(namxoc), tamp(namxat)
      integer startt, countt, start(2), count(2)

      SAVE time_id, ocdat_id, atdat_id
#  endif /* use_netcdf */
*
      DATA first/.true./
*
*     Work out subscript ranges and weights for the subareas
*     ======================================================
      if ( first ) then

        write(6,'(a)') ' '
        write(6,'(a)') '  Area-averaged diagnostics:'
        write(6,'(a)') '  =========================='

*       Read file specifying numbers, limits and names of areas
*       -------------------------------------------------------
        open (arunit, file='areas.limits', status='old')
*       Oceanic
*       -------
        read (arunit, *, end=100) nareoc
**      print *,' nareat: ',nareat
        if ( nareoc.gt.namxoc ) then
          write(6,'(a)') '  ERROR: nareoc > namxoc'
          write(6,'(a,2i5)') '  nareoc, namxoc = ',nareoc,namxoc
          write(6,'(a)')
     &          '  Increase parameter namxoc in areasubs_diag.F'
          write(6,'(a)') '  or reduce value of nareoc in areas.limits'
          write(6,'(a)') '  Program terminates in areavg'
          stop
        endif
        if ( nareoc.gt.9 ) then
          write(6,'(a)') '  WARNING: nareoc > 9'
          write(6,'(a)')
     &          '  Format stmnts in areasubs_diag.F need changing'
        endif
        read (arunit, *, err=100) (xlooc(ia),ia=1,nareoc)
**      print *,' xlooc: ',(xlooc(ia),ia=1,nareoc)
        read (arunit, *, err=100) (xhioc(ia),ia=1,nareoc)
**      print *,' xhioc: ',(xhioc(ia),ia=1,nareoc)
        read (arunit, *, err=100) (ylooc(ia),ia=1,nareoc)
**      print *,' ylooc: ',(ylooc(ia),ia=1,nareoc)
        read (arunit, *, err=100) (yhioc(ia),ia=1,nareoc)
**      print *,' yhioc: ',(yhioc(ia),ia=1,nareoc)
        read (arunit, '(9(8x,a3))') (areaoc(ia),ia=1,nareoc)
**      write(*,'(a,9(2x,a))') ' areaoc: ',(areaoc(ia),ia=1,nareoc)
*       Atmospheric
*       -----------
        read (arunit, *, err=100) nareat
**      print *,' nareat: ',nareat
        if ( nareat.gt.namxat ) then
          write(6,'(a)') '  ERROR: nareat > namxat'
          write(6,'(a,2i5)') '  nareat, namxat = ',nareat,namxat
          write(6,'(a)')
     &          '  Increase parameter namxat in areasubs_diag.F'
          write(6,'(a)') '  or reduce value of nareat in areas.limits'
          write(6,'(a)') '  Program terminates in areavg'
          stop
        endif
        if ( nareat.gt.9 ) then
          write(6,'(a)') '  WARNING: nareat > 9'
          write(6,'(a)')
     &          '  Format stmnts in areasubs_diag.F need changing'
        endif
        read (arunit, *, err=100) (xloat(ia),ia=1,nareat)
**      print *,' xloat: ',(xloat(ia),ia=1,nareat)
        read (arunit, *, err=100) (xhiat(ia),ia=1,nareat)
**      print *,' xhiat: ',(xhiat(ia),ia=1,nareat)
        read (arunit, *, err=100) (yloat(ia),ia=1,nareat)
**      print *,' yloat: ',(yloat(ia),ia=1,nareat)
        read (arunit, *, err=100) (yhiat(ia),ia=1,nareat)
**      print *,' yhiat: ',(yhiat(ia),ia=1,nareat)
        read (arunit, '(9(8x,a3))') (areaat(ia),ia=1,nareat)
**      write(*,'(a,9(2x,a))') ' areaat: ',(areaat(ia),ia=1,nareat)
        goto 200
  100   print *,' Error while reading file areas.limits in areavg'
        write(6,'(a)') '  Program terminates in areavg'
        stop
  200   continue
        close(arunit)

        write(6,'(a)') ' '
        write(6,'(a)') '  Oceanic subareas:'
        write(6,'(a)') '  -----------------'
*       Measure w.r.t. ocean coordinate system
        xoclim = 0.0d0
        yoclim = 0.0d0
*       Measure w.r.t. atmos. coordinate system
***     xoclim = (nx1-1)*dxa
***     yoclim = ypo(1)
        do m=1,nareoc
*         T points
*         --------
          rilo = 1.0d0 + ( xlooc(m) - 0.5d0*dxo - xoclim )/dxo
          rihi = 1.0d0 + ( xhioc(m) - 0.5d0*dxo - xoclim )/dxo
          frlo = mod( rilo, 1.0d0 )
          i1toc(m) = int( rilo )
          if ( frlo.ge.0.5d0 ) then
            frlo = frlo - 1.0d0
            i1toc(m) = i1toc(m) + 1
          endif
          fwtoc(m) = 0.5d0 - frlo
          frhi = mod( rihi, 1.0d0 )
          i2toc(m) = int( rihi )
          if ( frhi.gt.0.5d0 ) then
            frhi = frhi - 1.0d0
            i2toc(m) = i2toc(m) + 1
          endif
          fetoc(m) = 0.5d0 + frhi
          rjlo = 1.0d0 + ( ylooc(m) - 0.5d0*dyo - yoclim )/dyo
          rjhi = 1.0d0 + ( yhioc(m) - 0.5d0*dyo - yoclim )/dyo
          frlo = mod( rjlo, 1.0d0 )
          j1toc(m) = int( rjlo )
          if ( frlo.ge.0.5d0 ) then
            frlo = frlo - 1.0d0
            j1toc(m) = j1toc(m) + 1
          endif
          fstoc(m) = 0.5d0 - frlo
          frhi = mod( rjhi, 1.0d0 )
          j2toc(m) = int( rjhi )
          if ( frhi.gt.0.5d0 ) then
            frhi = frhi - 1.0d0
            j2toc(m) = j2toc(m) + 1
          endif
          fntoc(m) = 0.5d0 + frhi
*         p points
*         --------
          rilo = 1.0d0 + ( xlooc(m) - xoclim )/dxo
          rihi = 1.0d0 + ( xhioc(m) - xoclim )/dxo
          frlo = mod( rilo, 1.0d0 )
          i1poc(m) = int( rilo )
          if ( frlo.ge.0.5d0 ) then
            frlo = frlo - 1.0d0
            i1poc(m) = i1poc(m) + 1
          endif
          fwpoc(m) = 0.5d0 - frlo
          frhi = mod( rihi, 1.0d0 )
          i2poc(m) = int( rihi )
          if ( frhi.gt.0.5d0 ) then
            frhi = frhi - 1.0d0
            i2poc(m) = i2poc(m) + 1
          endif
          fepoc(m) = 0.5d0 + frhi
          rjlo = 1.0d0 + ( ylooc(m) - yoclim )/dyo
          rjhi = 1.0d0 + ( yhioc(m) - yoclim )/dyo
          frlo = mod( rjlo, 1.0d0 )
          j1poc(m) = int( rjlo )
          if ( frlo.ge.0.5d0 ) then
            frlo = frlo - 1.0d0
            j1poc(m) = j1poc(m) + 1
          endif
          fspoc(m) = 0.5d0 - frlo
          frhi = mod( rjhi, 1.0d0 )
          j2poc(m) = int( rjhi )
          if ( frhi.gt.0.5d0 ) then
            frhi = frhi - 1.0d0
            j2poc(m) = j2poc(m) + 1
          endif
          fnpoc(m) = 0.5d0 + frhi
          write(6,'(a)') ' '
          write(*,'(a,i3,a6)') '  m, area name       = ',
     &                         m,areaoc(m)
          write(*,'(a,4f9.1)') '  xlo, xhi, ylo, yhi = ',
     &                         1.0d-3*xlooc(m),1.0d-3*xhioc(m),
     &                         1.0d-3*ylooc(m),1.0d-3*yhioc(m)
          write(*,'(a,4i9)')   '  i1t, i2t, j1t, j2t = ',
     &                         i1toc(m),i2toc(m),j1toc(m),j2toc(m)
          write(*,'(a,4f9.3)') '             weights = ',
     &                         fwtoc(m),fetoc(m),fstoc(m),fntoc(m)
          write(*,'(a,4i9)')   '  i1p, i2p, j1p, j2p = ',
     &                         i1poc(m),i2poc(m),j1poc(m),j2poc(m)
          write(*,'(a,4f9.3)') '             weights = ',
     &                         fwpoc(m),fepoc(m),fspoc(m),fnpoc(m)
*
*         Check subarea limits against domain
*         -----------------------------------
          if ( ( xlooc(m) - xoclim ).lt.0.0d0 .or.
     &         ( xlooc(m) - xoclim ).gt.xlo .or.
     &         ( xhioc(m) - xoclim ).lt.0.0d0 .or.
     &         ( xhioc(m) - xoclim ).gt.xlo ) then
            print *,' limits outside x range for ocean area m = ',m
            print *,' xlooc(m), xhioc(m) = ',xlooc(m),xhioc(m)
            print *,' limits are : ',xoclim,xoclim+xlo
            print *,' program terminates in areavg'
            stop
          endif
          if ( ( ylooc(m) - yoclim ).lt.0.0d0 .or.
     &         ( ylooc(m) - yoclim ).gt.ylo .or.
     &         ( yhioc(m) - yoclim ).lt.0.0d0 .or.
     &         ( yhioc(m) - yoclim ).gt.ylo ) then
            print *,' limits outside y range for ocean area m = ',m
            print *,' ylooc(m), yhioc(m) = ',ylooc(m),yhioc(m)
            print *,' limits are : ',yoclim,yoclim+ylo
            print *,' program terminates in areavg'
            stop
          endif
*         Check subscript values
*         ----------------------
          if ( i1toc(m).lt.1 .or. i1toc(m).gt.nxto .or.
     &         i2toc(m).lt.1 .or. i2toc(m).gt.nxto ) then
            print *,' invalid subscript for ocean area m = ',m
            print *,' i1toc(m), i2toc(m) = ',i1toc(m),i2toc(m)
            print *,' program terminates in areavg'
            stop
          endif
          if ( j1toc(m).lt.1 .or. j1toc(m).gt.nyto .or.
     &         j2toc(m).lt.1 .or. j2toc(m).gt.nyto ) then
            print *,' invalid subscript for ocean area m = ',m
            print *,' j1toc(m), j2toc(m) = ',j1toc(m),j2toc(m)
            print *,' program terminates in areavg'
            stop
          endif
          if ( i1poc(m).lt.1 .or. i1poc(m).gt.nxpo .or.
     &         i2poc(m).lt.1 .or. i2poc(m).gt.nxpo ) then
            print *,' invalid subscript for ocean area m = ',m
            print *,' i1poc(m), i2poc(m) = ',i1poc(m),i2poc(m)
            print *,' program terminates in areavg'
            stop
          endif
          if ( j1poc(m).lt.1 .or. j1poc(m).gt.nypo .or.
     &         j2poc(m).lt.1 .or. j2poc(m).gt.nypo ) then
            print *,' invalid subscript for ocean area m = ',m
            print *,' j1poc(m), j2poc(m) = ',j1poc(m),j2poc(m)
            print *,' program terminates in areavg'
            stop
          endif
        enddo

        write(6,'(a)') ' '
        write(6,'(a)') '  Atmospheric subareas:'
        write(6,'(a)') '  ---------------------'
        do m=1,nareat
*         T points
*         --------
          rilo = 1.0d0 + ( xloat(m) - 0.5d0*dxa )/dxa
          rihi = 1.0d0 + ( xhiat(m) - 0.5d0*dxa )/dxa
          frlo = mod( rilo, 1.0d0 )
          i1tat(m) = int( rilo )
          if ( frlo.ge.0.5d0 ) then
            frlo = frlo - 1.0d0
            i1tat(m) = i1tat(m) + 1
          endif
          fwtat(m) = 0.5d0 - frlo
          frhi = mod( rihi, 1.0d0 )
          i2tat(m) = int( rihi )
          if ( frhi.gt.0.5d0 ) then
            frhi = frhi - 1.0d0
            i2tat(m) = i2tat(m) + 1
          endif
          fetat(m) = 0.5d0 + frhi
          rjlo = 1.0d0 + ( yloat(m) - 0.5d0*dya )/dya
          rjhi = 1.0d0 + ( yhiat(m) - 0.5d0*dya )/dya
          frlo = mod( rjlo, 1.0d0 )
          j1tat(m) = int( rjlo )
          if ( frlo.ge.0.5d0 ) then
            frlo = frlo - 1.0d0
            j1tat(m) = j1tat(m) + 1
          endif
          fstat(m) = 0.5d0 - frlo
          frhi = mod( rjhi, 1.0d0 )
          j2tat(m) = int( rjhi )
          if ( frhi.gt.0.5d0 ) then
            frhi = frhi - 1.0d0
            j2tat(m) = j2tat(m) + 1
          endif
          fntat(m) = 0.5d0 + frhi
*         p points
*         --------
          rilo = 1.0d0 + xloat(m)/dxa
          rihi = 1.0d0 + xhiat(m)/dxa
          frlo = mod( rilo, 1.0d0 )
          i1pat(m) = int( rilo )
          if ( frlo.ge.0.5d0 ) then
            frlo = frlo - 1.0d0
            i1pat(m) = i1pat(m) + 1
          endif
          fwpat(m) = 0.5d0 - frlo
          frhi = mod( rihi, 1.0d0 )
          i2pat(m) = int( rihi )
          if ( frhi.gt.0.5d0 ) then
            frhi = frhi - 1.0d0
            i2pat(m) = i2pat(m) + 1
          endif
          fepat(m) = 0.5d0 + frhi
          rjlo = 1.0d0 + yloat(m)/dya
          rjhi = 1.0d0 + yhiat(m)/dya
          frlo = mod( rjlo, 1.0d0 )
          j1pat(m) = int( rjlo )
          if ( frlo.ge.0.5d0 ) then
            frlo = frlo - 1.0d0
            j1pat(m) = j1pat(m) + 1
          endif
          fspat(m) = 0.5d0 - frlo
          frhi = mod( rjhi, 1.0d0 )
          j2pat(m) = int( rjhi )
          if ( frhi.gt.0.5d0 ) then
            frhi = frhi - 1.0d0
            j2pat(m) = j2pat(m) + 1
          endif
          fnpat(m) = 0.5d0 + frhi
          write(6,'(a)') ' '
          write(*,'(a,i3,a6)') '  m, area name       = ',
     &                         m,areaat(m)
          write(*,'(a,4f9.1)') '  xlo, xhi, ylo, yhi = ',
     &                         1.0d-3*xloat(m),1.0d-3*xhiat(m),
     &                         1.0d-3*yloat(m),1.0d-3*yhiat(m)
          write(*,'(a,4i9)')   '  i1t, i2t, j1t, j2t = ',
     &                         i1tat(m),i2tat(m),j1tat(m),j2tat(m)
          write(*,'(a,4f9.3)') '             weights = ',
     &                         fwtat(m),fetat(m),fstat(m),fntat(m)
          write(*,'(a,4i9)')   '  i1p, i2p, j1p, j2p = ',
     &                         i1pat(m),i2pat(m),j1pat(m),j2pat(m)
          write(*,'(a,4f9.3)') '             weights = ',
     &                         fwpat(m),fepat(m),fspat(m),fnpat(m)
*
*         Check subarea limits against domain
*         -----------------------------------
          if ( xloat(m).lt.0.0d0 .or. xloat(m).gt.xla .or.
     &         xhiat(m).lt.0.0d0 .or. xhiat(m).gt.xla ) then
            print *,' limits outside x range for atmos. area m = ',m
            print *,' xloat(m), xhiat(m) = ',xloat(m),xhiat(m)
            print *,' limits are : ',0.0d0,xla
            print *,' limits outside domain for atmos. area m = ',m
            print *,' program terminates in areavg'
            stop
          endif
          if ( yloat(m).lt.0.0d0 .or. yloat(m).gt.yla .or.
     &         yhiat(m).lt.0.0d0 .or. yhiat(m).gt.yla ) then
            print *,' limits outside y range for atmos. area m = ',m
            print *,' yloat(m), yhiat(m) = ',yloat(m),yhiat(m)
            print *,' limits are : ',0.0d0,yla
            print *,' program terminates in areavg'
            stop
          endif
        enddo

#  ifdef use_netcdf
!! Initialise netcdf file

!!  Define time dimension:
        status = nf_def_dim(arencid, 'time', numoutsteps, timedim)
        if (status.ne.NF_NOERR) call handle_err (status, subnam)

!!  Define a one-dimensional variable called 'time'
!!  which stores time data in years
        status = nf_def_var(arencid,'time',NF_FLOAT,1,timedim,time_id)
        if (status.ne.NF_NOERR) call handle_err (status, subnam)
        status = nf_put_att_text(arencid,time_id,'units',5,'years')
        if (status.ne.NF_NOERR) call handle_err (status, subnam)
        status = nf_put_att_text(arencid,time_id,'long_name',4,'Time')
        if (status.ne.NF_NOERR) call handle_err (status, subnam)

#    ifndef atmos_only
!!  Define nareoc dimension
        status = nf_def_dim(arencid, 'nareoc', nareoc, nocdim)
        if (status.ne.NF_NOERR) call handle_err (status, subnam)
!!  Define a 1-dimensional variable called 'nareoc'
        status=nf_def_var(arencid, 'nareoc', NF_INT, 1, nocdim, noc_id)
        if (status.ne.NF_NOERR) call handle_err (status, subnam)
        status = nf_put_att_text(arencid,noc_id,'units',2,' ')
        if (status.ne.NF_NOERR) call handle_err (status, subnam)
        status = nf_put_att_text(arencid,noc_id,'long_name',21,
     &       'Ocean averaging areas')
        if (status.ne.NF_NOERR) call handle_err (status, subnam)
!!  Define a two-dimensional variable called 'ocdata'
        edims(1) = nocdim
        edims(2) = timedim
        status = nf_def_var(arencid,'ocdata',NF_FLOAT,2,edims,ocdat_id)
        if (status.ne.NF_NOERR) call handle_err (status, subnam)
        status = nf_put_att_text(arencid,ocdat_id,'units',1,'K')
        if (status.ne.NF_NOERR) call handle_err (status, subnam)
        status = nf_put_att_text(arencid,ocdat_id,'long_name',31,
     &       'Area averaged ocean temperature')
        if (status.ne.NF_NOERR) call handle_err (status, subnam)
#    endif

#    ifndef ocean_only
!!  Define nareat dimension
        status = nf_def_dim(arencid, 'nareat', nareat, natdim)
        if (status.ne.NF_NOERR) call handle_err (status, subnam)
!!  Define a 1-dimensional variable called 'nareat'
        status=nf_def_var(arencid, 'nareat', NF_INT, 1, natdim, nat_id)
        if (status.ne.NF_NOERR) call handle_err (status, subnam)
        status = nf_put_att_text(arencid,nat_id,'units',2,' ')
        if (status.ne.NF_NOERR) call handle_err (status, subnam)
        status = nf_put_att_text(arencid,nat_id,'long_name',26,
     &       'Atmosphere averaging areas')
        if (status.ne.NF_NOERR) call handle_err (status, subnam)
!!  Define a two-dimensional variable called 'atdata'
        edims(1) = natdim
        edims(2) = timedim
        status = nf_def_var(arencid,'atdata',NF_FLOAT,2,edims,atdat_id)
        if (status.ne.NF_NOERR) call handle_err (status, subnam)
        status = nf_put_att_text(arencid,atdat_id,'units',1,'K')
        if (status.ne.NF_NOERR) call handle_err (status, subnam)
        status = nf_put_att_text(arencid,atdat_id,'long_name',36,
     &       'Area averaged atmosphere temperature')
        if (status.ne.NF_NOERR) call handle_err (status, subnam)
#    endif

!! Leave definition mode: entering data mode.
        status = nf_enddef(arencid)
        if (status.ne.NF_NOERR) call handle_err (status, subnam)

!! write ocean areas
#    ifndef atmos_only
        do m=1,nareoc
           tomp(m) = m
        enddo
        status = nf_put_var_int(arencid, noc_id, tomp)
        if (status.ne.NF_NOERR) call handle_err (status, subnam)
#    endif

!! write atmos areas
#    ifndef ocean_only
        do m=1,nareat
           tamp(m) = m
        enddo
        status = nf_put_var_int(arencid, nat_id, tamp)
        if (status.ne.NF_NOERR) call handle_err (status, subnam)
#    endif

!! NetCDF file all ready to go...
#  endif /* use_netcdf */

      endif
      first = .false.
*
*     Compute limited-area diagnostics
*     ================================
*     Oceanic
*     -------
#  ifndef atmos_only
      do m=1,nareoc
        call areint (sst, nxto, nyto, i1toc(m), i2toc(m), j1toc(m),
     &               j2toc(m), fwtoc(m), fetoc(m), fstoc(m), fntoc(m),
     &               tavoc(m))
      enddo

#  endif
*     Atmospheric
*     -----------
#  ifndef ocean_only
      do m=1,nareat
        call areint (ast, nxta, nyta, i1tat(m), i2tat(m), j1tat(m),
     &               j2tat(m), fwtat(m), fetat(m), fstat(m), fntat(m),
     &               tavat(m))
      enddo
#  endif

#  ifdef use_netcdf
!! Now write to NetCDF file
!! Store current time as part of 'time' vector
      startt = ntdone/nocmon + 1
      countt = 1
      status=nf_put_vara_double(arencid, time_id, startt, countt, tyrs)
      if (status.ne.NF_NOERR) call handle_err (status, subnam)

!! Store tavoc in 'ocdata'
#    ifndef atmos_only
      start(1) = 1
      start(2) = ntdone/nocmon+1
      count(1) = nareoc
      count(2) = 1
      status = nf_put_vara_double(arencid,ocdat_id,start,count,tavoc)
      if (status.ne.NF_NOERR) call handle_err (status, subnam)
#    endif

!! Store tavat in 'atdata'
#    ifndef ocean_only
      start(1) = 1
      start(2) = ntdone/nocmon+1
      count(1) = nareat
      count(2) = 1
      status = nf_put_vara_double(arencid,atdat_id,start,count,tavat)
      if (status.ne.NF_NOERR) call handle_err (status, subnam)
#    endif

#  endif /* use_netcdf */
#endif
*
      END SUBROUTINE areavg

#endif /* get_areav */
*
c***********************************************************************
*
#ifdef get_areav
*
      SUBROUTINE areint (array, nx, ny, ilo, ihi, jlo, jhi,
     &                   facw, face, facs, facn, answer)
*
*     Computes area integral of subsection of array(nx,ny) from
*     i=ilo to i=ihi and j=jlo to j=jhi. Weighting factors facw,
*     face, facs and facn are applied to the boundary points.
*     Modified version with reduced accumulator error
*
      IMPLICIT NONE
*
*     Subroutine arguments
      integer nx,ny,ilo,ihi,jlo,jhi
      double precision array(nx,ny),facw,face,facs,facn,answer
*
*     Local variables
      integer i,j
      double precision asums,wsums,asumi,wsumi,arsum,wtsum,asumn,wsumn
*
      arsum = 0.0d0
      wtsum = 0.0d0

*     Southern boundary
      asums = facw*array(ilo,jlo)
      wsums = facw
*     Northern boundary
      asumn = facw*array(ilo,jhi)
      wsumn = facw

*     Inner latitudes
!$OMP PARALLEL DEFAULT (NONE)
!$OMP&         PRIVATE (i,j,asumi,wsumi)
!$OMP&         SHARED  (jlo,jhi,ilo,ihi,facw,face,array,
!$OMP&                  arsum,wtsum,asums,wsums,asumn,wsumn)

!$OMP DO SCHEDULE (STATIC)
!$OMP&   REDUCTION(+:arsum,wtsum)
      do j=jlo+1,jhi-1
        asumi = facw*array(ilo,j)
        wsumi = facw
        do i=ilo+1,ihi-1
          asumi = asumi + array(i,j)
          wsumi = wsumi + 1.0d0
        enddo
        asumi = asumi + face*array(ihi,j)
        wsumi = wsumi + face
        arsum = arsum + asumi
        wtsum = wtsum + wsumi
      enddo
!$OMP END DO NOWAIT

!$OMP DO SCHEDULE (STATIC)
!$OMP&   REDUCTION(+:asums,wsums,asumn,wsumn)
      do i=ilo+1,ihi-1
*       Southern boundary
        asums = asums + array(i,jlo)
        wsums = wsums + 1.0d0
*       Northern boundary
        asumn = asumn + array(i,jhi)
        wsumn = wsumn + 1.0d0
      enddo
!$OMP END DO NOWAIT

!$OMP END PARALLEL

*     Southern boundary
      asums = asums + face*array(ihi,jlo)
      wsums = wsums + face
*     Northern boundary
      asumn = asumn + face*array(ihi,jhi)
      wsumn = wsumn + face

      answer = ( facs*asums + facn*asumn + arsum )
     &        /( facs*wsums + facn*wsumn + wtsum )

      END SUBROUTINE areint

#endif /* get_areav */
*
c***********************************************************************
*
      END MODULE areasubs
*
c***********************************************************************
