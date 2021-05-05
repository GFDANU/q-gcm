c***********************************************************************
c     Q-GCM Version 1.5.0 : module introduced 11/11/2017 by S. Kravtsov
c
c     modified 10/23/2020 to include additional output for the moist
c     version of the model
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
c
      MODULE hra_out

*     Contains subprograms for input, creation and output of high-
*     resolution atmospheric output in netCDF format.
*     For files created by the subroutines in this program file,
*     all real fields are single precision (float) format.
*
*     As all input/output from these routines is explicitly intended
*     to be in netCDF format, the corresponding #ifdef is used to
*     remove code to the widest possible extent, with configuration-
*     dependant #ifdefs being nested within the use_netcdf test.

      IMPLICIT NONE

      PRIVATE

      PUBLIC :: atnchr_init, atnchr_out

#ifdef use_netcdf

#  ifndef ocean_only

#    ifdef highres_output

*     Storage for identifiers for atmospheric netCDF dumps
      integer, PUBLIC, SAVE :: atphrid, atthrid
      integer, PRIVATE, SAVE :: ast_id, pa_id, eta_id,
     &  timap_id, timat_id

#    endif /* highres_output */

#  endif /* not ocean_only */

#endif /* use_netcdf */

      CONTAINS

c***********************************************************************


      SUBROUTINE atnchr_init

*     Modules
      USE parameters
#ifndef ocean_only
      USE atconst
#endif
      USE timinfo, ONLY : noutstepat
      USE nc_subs, ONLY : handle_err

      IMPLICIT NONE

#ifdef use_netcdf
      INCLUDE 'netcdf.inc'
#endif /* use_netcdf */
*
*     Subroutine arguments

#ifdef use_netcdf
#  ifndef ocean_only
#    ifdef highres_output

*     Local parameters
      character (len=*), parameter :: subnam = 'atnchr_init'

!!  netCDF variables used locally
      integer ncstat, timapdim, xapdim, yapdim, ladim, lamdim
      integer timatdim, xatdim, yatdim, axadims, adims(3), padims(4)
      integer xap_id, yap_id, xat_id, yat_id, la_id, lam_id

!!  Other variables used locally
      double precision xa(nxpa),ya(nypa),tmp(nla)
      integer i, iwk, mwk

!!  Define four dimensions: x, y, z, time
!!  Do pressure fields first
      ncstat = nf_def_dim(atphrid, 'time', noutstepat, timapdim)
      if ( ncstat.ne.NF_NOERR ) call handle_err (ncstat, subnam)

!!  Note that x and y are designed for subsampling
!!  Here are the x and y dimensions for the p-grid
      mwk = mod(nxpout, nska1)
      iwk = min(mwk,1) + (nxpout-mwk)/nska1
      ncstat = nf_def_dim(atphrid, 'xp', iwk, xapdim)
      if ( ncstat.ne.NF_NOERR ) call handle_err (ncstat, subnam)
      mwk = mod(nypout, nska1)
      iwk = min(mwk,1) + (nypout-mwk)/nska1
      ncstat = nf_def_dim(atphrid, 'yp', iwk, yapdim)
      if ( ncstat.ne.NF_NOERR ) call handle_err (ncstat, subnam)

!!  And here the z dimension
      ncstat = nf_def_dim(atphrid, 'z', nla, ladim)
      if ( ncstat.ne.NF_NOERR ) call handle_err (ncstat, subnam)
!!  Also define a zi dimension for the interfaces
      ncstat = nf_def_dim(atphrid, 'zi', nla-1, lamdim)
      if ( ncstat.ne.NF_NOERR ) call handle_err (ncstat, subnam)

!!  Now do it for the T-grid
      ncstat = nf_def_dim(atthrid, 'time', noutstepat, timatdim)
      if ( ncstat.ne.NF_NOERR ) call handle_err (ncstat, subnam)

!!  Here are the x and y dimensions for the T-grid
      mwk = mod(nxtout, nska1)
      iwk = min(mwk,1) + (nxtout-mwk)/nska1
      ncstat = nf_def_dim(atthrid, 'xt', iwk, xatdim)
      if ( ncstat.ne.NF_NOERR ) call handle_err (ncstat, subnam)
      mwk = mod(nytout, nska1)
      iwk = min(mwk,1) + (nytout-mwk)/nska1
      ncstat = nf_def_dim(atthrid, 'yt', iwk, yatdim)
      if ( ncstat.ne.NF_NOERR ) call handle_err (ncstat, subnam)

!!  Define a one-dimensional variable called
!!  'xp' which stores x gridpoints (p-grid)
      axadims = xapdim
      ncstat = nf_def_var(atphrid, 'xp', NF_FLOAT, 1, axadims, xap_id)
      if ( ncstat.ne.NF_NOERR ) call handle_err (ncstat, subnam)
      ncstat = nf_put_att_text(atphrid, xap_id, 'units', 2, 'km')
      if ( ncstat.ne.NF_NOERR ) call handle_err (ncstat, subnam)
      ncstat = nf_put_att_text(atphrid, xap_id, 'long_name',
     &            26, 'Atmosphere X axis (p-grid)')
      if ( ncstat.ne.NF_NOERR ) call handle_err (ncstat, subnam)

!!  Define a one-dimensional variable called
!!  'xt' which stores x gridpoints (T-grid)
      axadims = xatdim
      ncstat = nf_def_var(atthrid, 'xt', NF_FLOAT, 1, axadims, xat_id)
      if ( ncstat.ne.NF_NOERR ) call handle_err (ncstat, subnam)
      ncstat = nf_put_att_text(atthrid, xat_id, 'units', 2, 'km')
      if ( ncstat.ne.NF_NOERR ) call handle_err (ncstat, subnam)
      ncstat = nf_put_att_text(atthrid, xat_id, 'long_name',
     &            26, 'Atmosphere X axis (T-grid)')
      if ( ncstat.ne.NF_NOERR ) call handle_err (ncstat, subnam)

!!  Define a one-dimensional variable called
!!  'yp' which stores y gridpoints (p-grid)
      axadims = yapdim
      ncstat = nf_def_var(atphrid, 'yp', NF_FLOAT, 1, axadims, yap_id)
      if ( ncstat.ne.NF_NOERR ) call handle_err (ncstat, subnam)
      ncstat = nf_put_att_text(atphrid, yap_id, 'units', 2, 'km')
      if ( ncstat.ne.NF_NOERR ) call handle_err (ncstat, subnam)
      ncstat = nf_put_att_text(atphrid, yap_id, 'long_name',
     &            26, 'Atmosphere Y axis (p-grid)')
      if ( ncstat.ne.NF_NOERR ) call handle_err (ncstat, subnam)

!!  Define a one-dimensional variable called
!!  'yt' which stores y gridpoints (T-grid)
      axadims = yatdim
      ncstat = nf_def_var(atthrid, 'yt', NF_FLOAT, 1, axadims, yat_id)
      if ( ncstat.ne.NF_NOERR ) call handle_err (ncstat, subnam)
      ncstat = nf_put_att_text(atthrid, yat_id, 'units', 2, 'km')
      if ( ncstat.ne.NF_NOERR ) call handle_err (ncstat, subnam)
      ncstat = nf_put_att_text(atthrid, yat_id, 'long_name',
     &            26, 'Atmosphere Y axis (T-grid)')
      if ( ncstat.ne.NF_NOERR ) call handle_err (ncstat, subnam)

!!  Define a one-dimensional variable called
!!  'time' which stores time data in years
!!  p-grid:
      axadims = timapdim
      ncstat=nf_def_var(atphrid, 'time', NF_FLOAT, 1, axadims, timap_id)
      if ( ncstat.ne.NF_NOERR ) call handle_err (ncstat, subnam)
      ncstat = nf_put_att_text(atphrid, timap_id, 'units', 5, 'years')
      if ( ncstat.ne.NF_NOERR ) call handle_err (ncstat, subnam)
      ncstat=nf_put_att_text(atphrid,timap_id,'long_name',9,'Time axis')
      if ( ncstat.ne.NF_NOERR ) call handle_err (ncstat, subnam)
!!  T-grid:
      axadims = timatdim
      ncstat=nf_def_var(atthrid, 'time', NF_FLOAT, 1, axadims, timat_id)
      if ( ncstat.ne.NF_NOERR ) call handle_err (ncstat, subnam)
      ncstat = nf_put_att_text(atthrid, timat_id, 'units', 5, 'years')
      if ( ncstat.ne.NF_NOERR ) call handle_err (ncstat, subnam)
      ncstat=nf_put_att_text(atthrid,timat_id,'long_name',9,'Time axis')
      if ( ncstat.ne.NF_NOERR ) call handle_err (ncstat, subnam)

!!  Define a 1-dimensional variable called 'z'
!!  to store layer thicknesses; p-grid only
      axadims = ladim
      ncstat = nf_def_var(atphrid, 'z', NF_FLOAT, 1, axadims, la_id)
      if ( ncstat.ne.NF_NOERR ) call handle_err (ncstat, subnam)
      ncstat = nf_put_att_text(atphrid, la_id, 'units', 2, 'km')
      if ( ncstat.ne.NF_NOERR ) call handle_err (ncstat, subnam)
      ncstat = nf_put_att_text(atphrid, la_id, 'long_name',
     &            32, 'Atmosphere mid-layer height axis')
      if ( ncstat.ne.NF_NOERR ) call handle_err (ncstat, subnam)

!!  Define a 1-dimensional variable called 'zi'
!!  to store interface heights; p-grid only
      axadims = lamdim
      ncstat = nf_def_var(atphrid, 'zi', NF_FLOAT, 1, axadims, lam_id)
      if ( ncstat.ne.NF_NOERR ) call handle_err (ncstat, subnam)
      ncstat = nf_put_att_text(atphrid, lam_id, 'units', 2, 'km')
      if ( ncstat.ne.NF_NOERR ) call handle_err (ncstat, subnam)
      ncstat = nf_put_att_text(atphrid, lam_id, 'long_name',
     &            32, 'Atmosphere interface height axis')
      if ( ncstat.ne.NF_NOERR ) call handle_err (ncstat, subnam)

!!  Define a 3-dimensional variable called 'ast'
!!  to store atmosphere surface temperatures
        adims(1) = xatdim
        adims(2) = yatdim
        adims(3) = timatdim
        ncstat = nf_def_var(atthrid, 'ast', NF_FLOAT, 3, adims, ast_id)
        if ( ncstat.ne.NF_NOERR ) call handle_err (ncstat, subnam)
        ncstat = nf_put_att_text(atthrid, ast_id, 'units', 1, 'K')
        if ( ncstat.ne.NF_NOERR ) call handle_err (ncstat, subnam)
        ncstat = nf_put_att_text(atthrid, ast_id, 'long_name',
     &              30, 'Atmosphere surface temperature')
        if ( ncstat.ne.NF_NOERR ) call handle_err (ncstat, subnam)

!!  Define a 4-dimensional variable called 'p' to store pressure
        padims(1) = xapdim
        padims(2) = yapdim
        padims(3) = ladim
        padims(4) = timapdim
        ncstat = nf_def_var(atphrid, 'p', NF_FLOAT, 4, padims, pa_id)
        if ( ncstat.ne.NF_NOERR ) call handle_err (ncstat, subnam)
        ncstat = nf_put_att_text(atphrid, pa_id, 'units', 7, 'm^2/s^2')
        if ( ncstat.ne.NF_NOERR ) call handle_err (ncstat, subnam)
        ncstat = nf_put_att_text(atphrid, pa_id, 'long_name',
     &              27, 'Atmosphere dynamic pressure')
        if ( ncstat.ne.NF_NOERR ) call handle_err (ncstat, subnam)

!!  Define a 3-dimensional variable called 'hmixa'
!!  to store atmospheric specific humidity
        adims(1) = xatdim
        adims(2) = yatdim
        adims(3) = timatdim
        ncstat=nf_def_var(atthrid,'hmixa',NF_FLOAT,3,adims,eta_id)
        if ( ncstat.ne.NF_NOERR ) call handle_err (ncstat, subnam)
        ncstat = nf_put_att_text(atthrid,eta_id,'units',1,'kg kg^-1')
        if ( ncstat.ne.NF_NOERR ) call handle_err (ncstat, subnam)
        ncstat = nf_put_att_text(atthrid, eta_id, 'long_name',
     &              20, 'Mixed layer humidity')
        if ( ncstat.ne.NF_NOERR ) call handle_err (ncstat, subnam)

!!  Leave definition mode: entering data mode
      ncstat = nf_enddef(atphrid)
      if ( ncstat.ne.NF_NOERR ) then
        print *,' atnchr_init: enddef problem for atmos p-grid file'
        print *,' if problem is "One or more variable sizes'
        print *,' violate format constraints", either coarsen'
        print *,' spatial and/or temporal sampling of the data,'
        print *,' or create the file with "64-bit offset" format'
        call handle_err (ncstat, subnam)
      endif
      ncstat = nf_enddef(atthrid)
      if ( ncstat.ne.NF_NOERR ) then
        print *,' atnchr_init: enddef problem for atmos T-grid file'
        print *,' if problem is "One or more variable sizes'
        print *,' violate format constraints", either coarsen'
        print *,' spatial and/or temporal sampling of the data,'
        print *,' or create the file with "64-bit offset" format'
        call handle_err (ncstat, subnam)
      endif

!!  Calculate x gridpoints and store in 'x' arrays
!!  p-grid points
      mwk = mod(nxpout, nska1)
      iwk = min(mwk,1) + (nxpout-mwk)/nska1
      do i=1,iwk
        xa(i) = 1.0d-3*xpa(nx1-nxbd+(i-1)*nska1)
      enddo
      ncstat = nf_put_vara_double(atphrid, xap_id, 1, iwk, xa)
      if ( ncstat.ne.NF_NOERR ) call handle_err (ncstat, subnam)
!!  T-grid points
      mwk = mod(nxtout, nska1)
      iwk = min(mwk,1) + (nxtout-mwk)/nska1
      do i=1,iwk
        xa(i) = 1.0d-3*xta(nx1-nxbd+(i-1)*nska1)
      enddo
      ncstat = nf_put_vara_double(atthrid, xat_id, 1, iwk, xa)
      if ( ncstat.ne.NF_NOERR ) call handle_err (ncstat, subnam)

!!  Calculate y gridpoints and store in 'y' arrays
!!  p-grid points
      mwk = mod(nypout, nska1)
      iwk = min(mwk,1) + (nypout-mwk)/nska1
      do i=1,iwk
        ya(i) = 1.0d-3*ypa(ny1-nybd+(i-1)*nska1)
      enddo
      ncstat = nf_put_vara_double(atphrid, yap_id, 1, iwk, ya)
      if ( ncstat.ne.NF_NOERR ) call handle_err (ncstat, subnam)
!!  T-grid points
      mwk = mod(nytout, nska1)
      iwk = min(mwk,1) + (nytout-mwk)/nska1
      do i=1,iwk
        ya(i) = 1.0d-3*yta(ny1-nybd+(i-1)*nska1)
      enddo
      ncstat = nf_put_vara_double(atthrid, yat_id, 1, iwk, ya)
      if ( ncstat.ne.NF_NOERR ) call handle_err (ncstat, subnam)

!!  Convert mid-layer heights into km and store in 'z'
      tmp(1) = 0.5d-3*hat(1)
      do i=2,nla
        tmp(i) = tmp(i-1) + 0.5d-3*( hat(i-1) + hat(i) )
      enddo
      ncstat = nf_put_var_double(atphrid,la_id, tmp)
      if ( ncstat.ne.NF_NOERR ) call handle_err (ncstat, subnam)

!!  Convert interface heights into km and store in 'zi'
      tmp(1) = 1.0d-3*hat(1)
      do i=2,nla-1
        tmp(i) = tmp(i-1) + 1.0d-3*hat(i)
      enddo
      ncstat = nf_put_var_double(atphrid, lam_id, tmp)
      if ( ncstat.ne.NF_NOERR ) call handle_err (ncstat, subnam)

#    endif /* highres_output */
#  endif /* not ocean_only */
#endif /* use_netcdf */

      END SUBROUTINE atnchr_init

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


      SUBROUTINE atnchr_out 

*     Modules
      USE parameters
#ifdef use_netcdf
#  ifndef ocean_only
#    ifdef highres_output
      USE atconst
      USE atstate
      USE intrfac, ONLY : ast, hmixa
#    endif /* highres_output */
#  endif /* not ocean_only */
#endif /* use_netcdf */
      USE timinfo, ONLY : ntdone, noutat, tyrs
      USE nc_subs, ONLY : handle_err

      IMPLICIT NONE

#ifdef use_netcdf
      INCLUDE 'netcdf.inc'
#endif /* use_netcdf */
*
*     Subroutine arguments

#ifdef use_netcdf
#  ifndef ocean_only
#    ifdef highres_output

*     Local parameters
      character (len=*), parameter :: subnam = 'atnchr_out'

!!  netCDF variables used locally
      integer ncstat, startt, countt, start(3), count(3),
     &     startp(4), countp(4)

!!  Other variables used locally
      double precision wrk(nxpout*nypout*nla)
      integer i,j,k,ipwk,jpwk,itwk,jtwk,mwk
      integer nx2,ny2,nxbd2,nybd2,nska2

!!  Define local array-dimension parameters
      nx2=nx1
      ny2=ny1
      nxbd2=nxbd
      nybd2=nybd
      nska2=nska1

!!  Store current time as part of 'time' vector
      startt = ntdone/noutat + 1
      countt = 1
!!  Pressure file:
      ncstat=nf_put_vara_double(atphrid, timap_id, startt, countt, tyrs)
      if ( ncstat.ne.NF_NOERR ) call handle_err (ncstat, subnam)
!!  Temperature file:
      ncstat=nf_put_vara_double(atthrid, timat_id, startt, countt, tyrs)
      if ( ncstat.ne.NF_NOERR ) call handle_err (ncstat, subnam)

!!  Compute subsampling array indices
      mwk = mod(nxpout, nska1)
      ipwk = min(mwk,1) + (nxpout-mwk)/nska1
      mwk = mod(nxtout, nska1)
      itwk = min(mwk,1) + (nxtout-mwk)/nska1
      mwk = mod(nypout, nska1)
      jpwk = min(mwk,1) + (nypout-mwk)/nska1
      mwk = mod(nytout, nska1)
      jtwk = min(mwk,1) + (nytout-mwk)/nska1

!!  Store AST in 'ast'
         start(1) = 1
         start(2) = 1
         start(3) = ntdone/noutat + 1
         count(1) = itwk
         count(2) = jtwk
         count(3) = 1

!!  Looks ridiculous, but we now write out the sub-sampled
!!  ast array into wrk which is a very long vector.
!!  One needs to be careful about the order things are written
!!  out to get meaningful netCDF data here....
!$OMP    PARALLEL DO DEFAULT (NONE)
!$OMP&            PRIVATE (i,j)
!$OMP&            SHARED  (itwk,jtwk,wrk,ast,nska2,nx2,ny2,nxbd2,nybd2)
!$OMP&            SCHEDULE (STATIC)
         do j=1,jtwk
            do i=1,itwk
               wrk(i+itwk*(j-1)) = ast(nx2-nxbd2+(i-1)*nska2,
     &                                 ny2-nybd2+(j-1)*nska2)
            enddo
         enddo
!$OMP    END PARALLEL DO
         ncstat = nf_put_vara_double(atthrid, ast_id, start, count, wrk)
         if ( ncstat.ne.NF_NOERR ) call handle_err (ncstat, subnam)

!!  Interpolate pressure onto T-grid points by averaging,
!!  and then write to 'p'. Note that we once again use the
!!  long wrk vector here.
         startp(1) = 1
         startp(2) = 1
         startp(3) = 1
         startp(4) = ntdone/noutat + 1
         countp(1) = ipwk
         countp(2) = jpwk
         countp(3) = nla
         countp(4) = 1
         do k=1,nla
!$OMP       PARALLEL DO DEFAULT (NONE)
!$OMP&             PRIVATE (i,j)
!$OMP&             SHARED  (ipwk,jpwk,wrk,pa,nska2,nx2,ny2,nxbd2,nybd2,k)
!$OMP&             SCHEDULE (STATIC)
            do j=1,jpwk
               do i=1,ipwk
                  wrk(i+ipwk*(j-1)+ipwk*jpwk*(k-1))=
     &                 pa(nx2-nxbd2+(i-1)*nska2,ny2-nybd2+(j-1)*nska2,k)
               enddo
            enddo
!$OMP       END PARALLEL DO
         enddo
         ncstat=nf_put_vara_double(atphrid, pa_id, startp, countp, wrk)
         if ( ncstat.ne.NF_NOERR ) call handle_err (ncstat, subnam)

!!  hmixa in hmixa
         start(1) = 1
         start(2) = 1
         start(3) = ntdone/noutat + 1
         count(1) = itwk
         count(2) = jtwk
         count(3) = 1
!$OMP    PARALLEL DO DEFAULT (NONE)
!$OMP&            PRIVATE (i,j)
!$OMP&            SHARED  (itwk,jtwk,wrk,hmixa,nska2,nx2,ny2,nxbd2,nybd2)
!$OMP&            SCHEDULE (STATIC)
         do j=1,jtwk
            do i=1,itwk
               wrk(i+itwk*(j-1)) = 
     &               hmixa(nx2-nxbd2+(i-1)*nska2,ny2-nybd2+(j-1)*nska2)
            enddo
         enddo
!$OMP    END PARALLEL DO
         ncstat = nf_put_vara_double(atthrid, eta_id, start, count, wrk)
         if ( ncstat.ne.NF_NOERR ) call handle_err (ncstat, subnam)

#    endif /* highres_output */
#  endif /* not ocean_only */
#endif /* use_netcdf */

      END SUBROUTINE atnchr_out

c
c***********************************************************************
*
      END MODULE hra_out
*
c***********************************************************************
