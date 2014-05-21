c***********************************************************************
c     Q-GCM Version 1.5.0 : last modified 18/09/2013
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
      MODULE qocdiag

*     Contains subprograms for computing and outputting (in netCDF format)
*     the time-derivative of ocean vorticity qo, and its component terms.
*     Output interval is noutoc atmospheric timesteps = odiday days
*     (controlled by when qocdiag_out is called from the main program)

      IMPLICIT NONE

      PRIVATE

      PUBLIC :: qocdiag_init, qocdiag_out

      integer, PUBLIC, SAVE  :: qocncid
      integer, PRIVATE, SAVE :: dqdt_id, qotjac_id, qt2dif_id,
     &                          qt4dif_id, qotent_id, timop_id
      logical, PRIVATE, SAVE :: qt2dif_out
*     Controls whether the 2nd order diffusive contribution is written
*     out, since this is usually exactly zero by choice of ah2oc

      CONTAINS

c***********************************************************************

      SUBROUTINE qocdiag_init (nsko)

*     Modules
      USE parameters
#ifndef atmos_only
#  ifdef qoc_diag
#    ifdef use_netcdf
      USE occonst
#    endif /* use_netcdf */
#  endif
#endif
      USE timinfo, ONLY : noutstepoc
      USE nc_subs, ONLY : handle_err

      IMPLICIT NONE

#ifdef use_netcdf
      INCLUDE 'netcdf.inc'
#endif /* use_netcdf */
*
*     Subroutine arguments
      integer, INTENT(IN) :: nsko

#ifndef atmos_only
#  ifdef qoc_diag
#    ifdef use_netcdf

*     Local parameters
      character (len=*), parameter :: subnam = 'qocdiag_init'

!!    netCDF variables used locally
      integer ncstat, timopdim, xopdim, yopdim, lodim
      integer axodims, podims(4), xop_id, yop_id, lo_id

!!    Other variables used locally
      double precision xo(nxpo), yo(nypo), zo(nlo)
      integer i, j, k, iwk, mwk

!!    Definition section: define dimensions and variables
!!    ===================================================
!!    Define four dimensions: x, y, z, time

*     Decide if any variables need not be defined
*     We usually run with ah2oc = 0.0 at all levels
*     If the value at any level is nonzero, then
*     output the contribution of this term to dq/dt
      qt2dif_out = .false.
      do k=1,nlo
        if ( ah2oc(k).ne.0.0d0 ) qt2dif_out = .true.
      enddo

!!    Dimension definitions for p-grid output file
!!    --------------------------------------------
      ncstat = nf_def_dim(qocncid, 'time', noutstepoc, timopdim)
      if ( ncstat.ne.NF_NOERR ) call handle_err (ncstat, subnam)

!!    x and y dimensions for the (subsampled) p-grid
      mwk = mod(nxpo,nsko)
      iwk = min(mwk,1) + (nxpo-mwk)/nsko
      ncstat = nf_def_dim(qocncid, 'xp', iwk, xopdim)
      if ( ncstat.ne.NF_NOERR ) call handle_err (ncstat, subnam)
      mwk = mod(nypo,nsko)
      iwk = min(mwk,1) + (nypo-mwk)/nsko
      ncstat = nf_def_dim(qocncid, 'yp', iwk, yopdim)
      if ( ncstat.ne.NF_NOERR ) call handle_err (ncstat, subnam)

!!    And the z dimension
      ncstat = nf_def_dim(qocncid, 'z', nlo, lodim)
      if ( ncstat.ne.NF_NOERR ) call handle_err (ncstat, subnam)

!!    Grid variable definitions for p-grid files
!!    ------------------------------------------
!!    Define a one-dimensional variable called
!!    'xp' which stores x gridpoints (p-grid)
      axodims = xopdim
      ncstat = nf_def_var(qocncid, 'xp', NF_FLOAT, 1, axodims, xop_id)
      if ( ncstat.ne.NF_NOERR ) call handle_err (ncstat, subnam)
      ncstat = nf_put_att_text(qocncid, xop_id, 'units', 2, 'km')
      if ( ncstat.ne.NF_NOERR ) call handle_err (ncstat, subnam)
      ncstat = nf_put_att_text(qocncid, xop_id, 'long_name',
     &            21, 'Ocean X axis (p-grid)')
      if ( ncstat.ne.NF_NOERR ) call handle_err (ncstat, subnam)

!!    Define a one-dimensional variable called
!!    'yp' which stores y gridpoints (p-grid)
      axodims = yopdim
      ncstat = nf_def_var(qocncid, 'yp', NF_FLOAT, 1, axodims, yop_id)
      if ( ncstat.ne.NF_NOERR ) call handle_err (ncstat, subnam)
      ncstat = nf_put_att_text(qocncid, yop_id, 'units', 2, 'km')
      if ( ncstat.ne.NF_NOERR ) call handle_err (ncstat, subnam)
      ncstat = nf_put_att_text(qocncid, yop_id, 'long_name',
     &            21, 'Ocean Y axis (p-grid)')
      if ( ncstat.ne.NF_NOERR ) call handle_err (ncstat, subnam)

!!    Define a one-dimensional variable called
!!    'time' which stores time data in years
      axodims = timopdim
      ncstat=nf_def_var(qocncid, 'time', NF_FLOAT, 1, axodims, timop_id)
      if ( ncstat.ne.NF_NOERR ) call handle_err (ncstat, subnam)
      ncstat = nf_put_att_text(qocncid, timop_id, 'units', 5, 'years')
      if ( ncstat.ne.NF_NOERR ) call handle_err (ncstat, subnam)
      ncstat=nf_put_att_text(qocncid,timop_id,'long_name',9,'Time axis')
      if ( ncstat.ne.NF_NOERR ) call handle_err (ncstat, subnam)

!!    Define a 1-dimensional variable called 'z'
!!    to store mid-layer depths (in km)
      axodims = lodim
      ncstat = nf_def_var(qocncid, 'z', NF_FLOAT, 1, axodims, lo_id)
      if ( ncstat.ne.NF_NOERR ) call handle_err (ncstat, subnam)
      ncstat = nf_put_att_text(qocncid, lo_id, 'units', 2, 'km')
      if ( ncstat.ne.NF_NOERR ) call handle_err (ncstat, subnam)
      ncstat = nf_put_att_text(qocncid, lo_id, 'long_name',
     &            26, 'Ocean mid-layer depth axis')
      if ( ncstat.ne.NF_NOERR ) call handle_err (ncstat, subnam)

!!    Define a 4-dimensional variable called 'dqdt'
!!    to store the time derivative of qo
      podims(1) = xopdim
      podims(2) = yopdim
      podims(3) = lodim
      podims(4) = timopdim
      ncstat = nf_def_var(qocncid, 'dqdt', NF_FLOAT, 4, podims, dqdt_id)
      if ( ncstat.ne.NF_NOERR ) then
        print *,' problem defining dqdt in ocean q-diag file'
        call handle_err (ncstat, subnam)
      endif
      ncstat = nf_put_att_text(qocncid, dqdt_id, 'units', 4, 's^-2')
      if ( ncstat.ne.NF_NOERR ) call handle_err (ncstat, subnam)
      ncstat = nf_put_att_text(qocncid, dqdt_id, 'long_name',
     &            23, 'd/dt(oceanic vorticity)')
      if ( ncstat.ne.NF_NOERR ) call handle_err (ncstat, subnam)

!!    Define a 4-dimensional variable called 'qotjac'
!!    to store advective (Jacobian) component of d(qo)/dt
      podims(1) = xopdim
      podims(2) = yopdim
      podims(3) = lodim
      podims(4) = timopdim
      ncstat = nf_def_var(qocncid, 'qotjac', NF_FLOAT,
     &                    4, podims, qotjac_id)
      if ( ncstat.ne.NF_NOERR ) then
        print *,' problem defining qotjac in ocean q-diag file'
        call handle_err (ncstat, subnam)
      endif
      ncstat = nf_put_att_text(qocncid, qotjac_id, 'units', 4, 's^-2')
      if ( ncstat.ne.NF_NOERR ) call handle_err (ncstat, subnam)
      ncstat = nf_put_att_text(qocncid, qotjac_id, 'long_name',
     &            31, 'advective component of d(qo)/dt')
      if ( ncstat.ne.NF_NOERR ) call handle_err (ncstat, subnam)

!!    Define a 4-dimensional variable called 'qt2dif'
!!    to store dissipation (2nd order) component of d(qo)/dt
      if ( qt2dif_out ) then
        podims(1) = xopdim
        podims(2) = yopdim
        podims(3) = lodim
        podims(4) = timopdim
        ncstat = nf_def_var(qocncid, 'qt2dif', NF_FLOAT,
     &                      4, podims, qt2dif_id)
        if ( ncstat.ne.NF_NOERR ) then
          print *,' problem defining qt2dif in ocean q-diag file'
          call handle_err (ncstat, subnam)
        endif
        ncstat = nf_put_att_text(qocncid, qt2dif_id, 'units', 4, 's^-2')
        if ( ncstat.ne.NF_NOERR ) call handle_err (ncstat, subnam)
        ncstat = nf_put_att_text(qocncid, qt2dif_id, 'long_name',
     &              45, 'dissipative component of d(qo)/dt (2nd order)')
        if ( ncstat.ne.NF_NOERR ) call handle_err (ncstat, subnam)
      endif

!!    Define a 4-dimensional variable called 'qt4dif'
!!    to store dissipation (4th order) component of d(qo)/dt
      podims(1) = xopdim
      podims(2) = yopdim
      podims(3) = lodim
      podims(4) = timopdim
      ncstat = nf_def_var(qocncid, 'qt4dif', NF_FLOAT,
     &                    4, podims, qt4dif_id)
      if ( ncstat.ne.NF_NOERR ) then
        print *,' problem defining qt4dif in ocean q-diag file'
        call handle_err (ncstat, subnam)
      endif
      ncstat = nf_put_att_text(qocncid, qt4dif_id, 'units', 4, 's^-2')
      if ( ncstat.ne.NF_NOERR ) call handle_err (ncstat, subnam)
      ncstat = nf_put_att_text(qocncid, qt4dif_id, 'long_name',
     &            45, 'dissipative component of d(qo)/dt (4th order)')
      if ( ncstat.ne.NF_NOERR ) call handle_err (ncstat, subnam)

!!    Define a 4-dimensional variable called 'qotent'
!!    to store entrainment (forcing) component of d(qo)/dt
      podims(1) = xopdim
      podims(2) = yopdim
      podims(3) = lodim
      podims(4) = timopdim
      ncstat = nf_def_var(qocncid, 'qotent', NF_FLOAT,
     &                    4, podims, qotent_id)
      if ( ncstat.ne.NF_NOERR ) then
        print *,' problem defining qotent in ocean q-diag file'
        call handle_err (ncstat, subnam)
      endif
      ncstat = nf_put_att_text(qocncid, qotent_id, 'units', 4, 's^-2')
      if ( ncstat.ne.NF_NOERR ) call handle_err (ncstat, subnam)
      ncstat = nf_put_att_text(qocncid, qotent_id, 'long_name',
     &            43, 'entrainment (forcing) component of d(qo)/dt')
      if ( ncstat.ne.NF_NOERR ) call handle_err (ncstat, subnam)

!!    Leave definition mode: entering data mode
!!    -----------------------------------------
      ncstat = nf_enddef(qocncid)
      if ( ncstat.ne.NF_NOERR ) then
        print *,' qocdiag_init: enddef problem for ocean q-diag file'
        print *,' if problem is "One or more variable sizes'
        print *,' violate format constraints", either coarsen'
        print *,' spatial and/or temporal sampling of the data,'
        print *,' or create the file with "64-bit offset" format'
        call handle_err (ncstat, subnam)
      endif

!!    Calculate x gridpoints and store in 'x' arrays
!!    p-grid points
      mwk = mod(nxpo,nsko)
      iwk = min(mwk,1) + (nxpo-mwk)/nsko
      do i=1,iwk
        xo(i) = 1.0d-3*( xpo(1+(i-1)*nsko) - xpo(1) )
      enddo
      ncstat = nf_put_vara_double(qocncid, xop_id, 1, iwk, xo)
      if ( ncstat.ne.NF_NOERR ) call handle_err (ncstat, subnam)

!!    Calculate y gridpoints and store in 'y' arrays
!!    p-grid points
      mwk = mod(nypo,nsko)
      iwk = min(mwk,1) + (nypo-mwk)/nsko
      do j=1,iwk
        yo(j) = 1.0d-3*( ypo(1+(j-1)*nsko) - ypo(1) )
      enddo
      ncstat = nf_put_vara_double(qocncid, yop_id, 1, iwk, yo)
      if ( ncstat.ne.NF_NOERR ) call handle_err (ncstat, subnam)

!!    Convert mid-layer depths into km and store in 'z'
      zo(1) = 0.5d-3*hoc(1)
      do k=2,nlo
        zo(k) = zo(k-1) + 0.5d-3*( hoc(k-1) + hoc(k) )
      enddo
      ncstat = nf_put_var_double(qocncid, lo_id, zo)
      if ( ncstat.ne.NF_NOERR ) call handle_err (ncstat, subnam)
!!!!! put operation seems to lack some arguments!!!
!!!!! This may be because it's a put_var not put_vara!!!

#    endif /* use_netcdf */
#  endif /* qoc_diag */
#endif /* not atmos_only */

      END SUBROUTINE qocdiag_init

c***********************************************************************

      SUBROUTINE qocdiag_out (nsko)

*     Compute oceanic dq/dt and its component terms. The code to do this
*     is borrowed from subroutines qgostep and ocadif in qgosubs.F.
*     An alternative approach would be to be to save the terms in
*     those routines, but this would (a) complicate them, and (b)
*     entail declaring much more storage for transfer. Here the arrays are
*     invoked only when necessary, at the expense of some recomputation.

*     Modules
      USE parameters
#ifndef atmos_only
#  ifdef qoc_diag
#    ifdef use_netcdf
      USE occonst
      USE ocstate
#    endif /* use_netcdf */
#  endif
#endif
      USE timinfo, ONLY : ntdone, noutoc, tyrs
      USE nc_subs, ONLY : handle_err

      IMPLICIT NONE

#ifdef use_netcdf
      INCLUDE 'netcdf.inc'
#endif /* use_netcdf */
*
*     Subroutine arguments
      integer, INTENT(IN) :: nsko

#ifndef atmos_only
#  ifdef qoc_diag
#    ifdef use_netcdf

*     Local parameters
      character (len=*), parameter :: subnam = 'qocdiag_out'

!!    netCDF variables used locally
!!    start4, count4 are for 4-D variables
      integer ncstat,startt,countt,start4(4),count4(4)

!!    Other variables used locally
      integer i,j,k,ipwk,jpwk,mwk
      double precision adfaco,bcfaco,fohfac(nlo),bdrfac,ah2fac,
     &                 ah4fac,del2p(nxpo,nypo),del4p(nxpo,nypo),
     &                 qt2dif(nxpo,nypo),qt4dif(nxpo,nypo),d6p,
     &                 qotjac(nxpo,nypo),qotent(nxpo,nypo),
     &                 dqdt(nxpo,nypo),wrk(nxpo*nypo),rdto

!!    Store current time as part of 'time' vector
      startt = ntdone/noutoc + 1
      countt = 1
      ncstat=nf_put_vara_double(qocncid, timop_id, startt, countt, tyrs)
      if ( ncstat.ne.NF_NOERR ) call handle_err (ncstat, subnam)

!!    Compute subsampling array indices
      mwk = mod(nxpo, nsko)
      ipwk = min(mwk,1) + (nxpo-mwk)/nsko
      mwk = mod(nypo, nsko)
      jpwk = min(mwk,1) + (nypo-mwk)/nsko

*     Compute vorticity derivative and components with code from qgosubs.F
*     --------------------------------------------------------------------
      adfaco = 1.0d0/(12.0d0*dxo*dyo*fnot)
*     Version with nondimensional bccooc
      bcfaco = bccooc*dxom2/( 0.5d0*bccooc + 1.0d0 )
      do k=1,nlo
        fohfac(k) = fnot/hoc(k)
      enddo
      bdrfac = 0.5d0*sign(1.0d0,fnot)*delek/hoc(nlo)
      rdto = 1.0d0/dto

      do k=1,nlo

        ah2fac = ah2oc(k)/fnot
        ah4fac = ah4oc(k)/fnot
*       Set netCDF output counters once only, since they
*       are the same for all fields in any given layer
        start4(1) = 1
        start4(2) = 1
        start4(3) = k
        start4(4) = ntdone/noutoc + 1
        count4(1) = ipwk
        count4(2) = jpwk
        count4(3) = 1
        count4(4) = 1

!$OMP   PARALLEL DEFAULT (NONE)
!$OMP&           PRIVATE (i,j,d6p)
!$OMP&           SHARED  (k,del2p,bcfaco,pom,dxom2,del4p,qt2dif,qt4dif,
!$OMP&                    qotjac,qotent,dqdt,po,qo,fohfac,wekpo,entoc,
!$OMP&                    bdrfac,ah2fac,ah4fac,start4,count4,adfaco,
!$OMP&                    rdto,ipwk,jpwk,wrk,nsko,ncstat,qocncid,
!$OMP&                    dqdt_id,qotjac_id,qt2dif_id,qt4dif_id,
!$OMP&                    qotent_id,qom,qt2dif_out)

*       Compute Del-sqd(p) at previous time level for dissipation
*       ---------------------------------------------------------
*       N & S boundaries (including corners) - mixed BCs
*       No x-derivative term because p is constant along boundary
!$OMP   DO SCHEDULE (STATIC)
        do i=1,nxpo
          del2p(i,  1 ) = bcfaco*( pom(i,   2  ,k) - pom(i,  1 ,k) )
          del2p(i,nypo) = bcfaco*( pom(i,nypo-1,k) - pom(i,nypo,k) )
        enddo
!$OMP   END DO NOWAIT

!$OMP   DO SCHEDULE (STATIC)
        do j=2,nypo-1
*         Western boundary (i=1)
#  ifdef cyclic_ocean
*         Zonally cyclic ocean
          del2p(1,j) = (  pom(1,j-1,k) + pom(nxpo-1,j,k) + pom(2,j,k)
     &                  + pom(1,j+1,k) - 4.0d0*pom(1,j,k) )*dxom2
#  else
*         Finite box ocean - mixed BCs
          del2p( 1 ,j) = bcfaco*( pom(  2  ,j,k) - pom( 1 ,j,k) )
#  endif /* cyclic_ocean */
*         Inner points (standard case)
          do i=2,nxpo-1
            del2p(i,j) = (  pom(i,j-1,k) + pom(i-1,j,k) + pom(i+1,j,k)
     &                    + pom(i,j+1,k) - 4.0d0*pom(i,j,k) )*dxom2
          enddo
*         Eastern boundary (i=nxpo)
#  ifdef cyclic_ocean
*         Zonally cyclic ocean
          del2p(nxpo,j) = del2p(1,j)
#  else
*         Finite box ocean - mixed BCs
          del2p(nxpo,j) = bcfaco*( pom(nxpo-1,j,k) - pom(nxpo,j,k) )
#  endif /* cyclic_ocean */
        enddo
!$OMP   END DO

*       Compute Del-4th(p) at previous time level for dissipation
*       ---------------------------------------------------------
*       N & S boundaries (including corners) - mixed BCs
*       No x-derivative term because p is constant along boundary
!$OMP   DO SCHEDULE (STATIC)
        do i=1,nxpo
          del4p(i,  1 ) = bcfaco*( del2p(i,   2  ) - del2p(i,  1 ) )
          del4p(i,nypo) = bcfaco*( del2p(i,nypo-1) - del2p(i,nypo) )
        enddo
!$OMP   END DO NOWAIT

!$OMP   DO SCHEDULE (STATIC)
        do j=2,nypo-1
*         Western boundary (i=1)
#  ifdef cyclic_ocean
*         Zonally cyclic ocean
          del4p(1,j) = (  del2p(1,j-1) + del2p(nxpo-1,j) + del2p(2,j)
     &                  + del2p(1,j+1) - 4.0d0*del2p(1,j) )*dxom2
#  else
*         Finite box ocean - mixed BCs
          del4p( 1 ,j) = bcfaco*( del2p(  2  ,j) - del2p( 1 ,j) )
#  endif /* cyclic_ocean */
*         Inner points (standard case)
          do i=2,nxpo-1
            del4p(i,j) = (  del2p(i,j-1) + del2p(i-1,j) + del2p(i+1,j)
     &                    + del2p(i,j+1) - 4.0d0*del2p(i,j) )*dxom2
          enddo
*         Eastern boundary (i=nxpo)
#  ifdef cyclic_ocean
*         Zonally cyclic ocean
          del4p(nxpo,j) = del4p(1,j)
#  else
*         Finite box ocean - mixed BCs
          del4p(nxpo,j) = bcfaco*( del2p(nxpo-1,j) - del2p(nxpo,j) )
#  endif /* cyclic_ocean */
        enddo
!$OMP   END DO

*       Zero all the results arrays as a precaution
*       -------------------------------------------
*       This can be removed when the boundary terms are computed
!$OMP   DO SCHEDULE (STATIC)
        do j=1,nypo
          do i=1,nxpo
            qt2dif(i,j) = 0.0d0
            qt4dif(i,j) = 0.0d0
            qotjac(i,j) = 0.0d0
            qotent(i,j) = 0.0d0
            dqdt(i,j)   = 0.0d0
          enddo
        enddo
!$OMP   END DO

*       Compute advective (Jacobian) and diffusive contributions to dq/dt
*       -----------------------------------------------------------------
*       Compute Del-6th(p) except at N & S boundaries
*       No need to apply N & S boundary conditions
!$OMP   DO SCHEDULE (STATIC)
        do j=2,nypo-1
*         Western boundary (i=1)
#  ifdef cyclic_ocean
*         Zonally cyclic ocean
          d6p = dxom2*(  del4p(1,j-1) + del4p(nxpo-1,j) + del4p(2,j)
     &                 + del4p(1,j+1) - 4.0d0*del4p(1,j) )
          qt2dif(1,j) = ah2fac*del4p(1,j)
          qt4dif(1,j) = -ah4fac*d6p
          qotjac(1,j) = adfaco*
     & ( (qo( 2 , j ,k)-qo(nxpo-1, j ,k))*(po( 1 ,j+1,k)-po(  1 ,j-1,k))
     &  +(qo( 1 ,j-1,k)-qo(  1 ,j+1,k))*(po( 2 , j ,k)-po(nxpo-1, j ,k))
     &  + qo(   2  , j ,k)*( po(   2  ,j+1,k) - po(   2  ,j-1,k) )
     &  - qo(nxpo-1, j ,k)*( po(nxpo-1,j+1,k) - po(nxpo-1,j-1,k) )
     &  - qo(   1  ,j+1,k)*( po(   2  ,j+1,k) - po(nxpo-1,j+1,k) )
     &  + qo(   1  ,j-1,k)*( po(   2  ,j-1,k) - po(nxpo-1,j-1,k) )
     &  + po(   1  ,j+1,k)*( qo(   2  ,j+1,k) - qo(nxpo-1,j+1,k) )
     &  - po(   1  ,j-1,k)*( qo(   2  ,j-1,k) - qo(nxpo-1,j-1,k) )
     &  - po(   2  , j ,k)*( qo(   2  ,j+1,k) - qo(   2  ,j-1,k) )
     &  + po(nxpo-1, j ,k)*( qo(nxpo-1,j+1,k) - qo(nxpo-1,j-1,k) ) )
*         (layer-dependent) entrainment/forcing term
          if ( k.eq.1 ) then
            qotent(1,j) = fohfac(1)*( wekpo(1,j) - entoc(1,j) )
           else if ( k.eq.2 ) then
            qotent(1,j) = fohfac(2)*entoc(1,j)
           else
            qotent(1,j) = 0.0d0
          endif
*         Add bottom drag contribution
          if ( k.eq.nlo ) then
            qotent(1,j) = qotent(1,j) - bdrfac*del2p(1,j)
          endif
          dqdt(1,j) =  qotjac(1,j) + qt2dif(1,j)
     &               + qt4dif(1,j) + qotent(1,j)
#  else
*         Finite box ocean
          qt2dif(1,j) = 0.0d0
          qt4dif(1,j) = 0.0d0
          qotjac(1,j) = 0.0d0
          qotent(1,j) = 0.0d0
          dqdt(1,j)   = rdto*( qo(1,j,k) - qom(1,j,k) )
#  endif /* cyclic_ocean */
*         Inner points (standard case)
          do i=2,nxpo-1
            d6p = dxom2*(  del4p(i,j-1) + del4p(i-1,j) + del4p(i+1,j)
     &                   + del4p(i,j+1) - 4.0d0*del4p(i,j) )
            qt2dif(i,j) = ah2fac*del4p(i,j)
            qt4dif(i,j) = -ah4fac*d6p
            qotjac(i,j) = adfaco*
     &   ( (qo(i+1, j ,k)-qo(i-1, j ,k))*(po( i ,j+1,k)-po( i ,j-1,k))
     &    +(qo( i ,j-1,k)-qo( i ,j+1,k))*(po(i+1, j ,k)-po(i-1, j ,k))
     &    + qo(i+1, j ,k)*( po(i+1,j+1,k) - po(i+1,j-1,k) )
     &    - qo(i-1, j ,k)*( po(i-1,j+1,k) - po(i-1,j-1,k) )
     &    - qo( i ,j+1,k)*( po(i+1,j+1,k) - po(i-1,j+1,k) )
     &    + qo( i ,j-1,k)*( po(i+1,j-1,k) - po(i-1,j-1,k) )
     &    + po( i ,j+1,k)*( qo(i+1,j+1,k) - qo(i-1,j+1,k) )
     &    - po( i ,j-1,k)*( qo(i+1,j-1,k) - qo(i-1,j-1,k) )
     &    - po(i+1, j ,k)*( qo(i+1,j+1,k) - qo(i+1,j-1,k) )
     &    + po(i-1, j ,k)*( qo(i-1,j+1,k) - qo(i-1,j-1,k) ) )
*           (layer-dependent) entrainment/forcing term
            if ( k.eq.1 ) then
              qotent(i,j) = fohfac(1)*( wekpo(i,j) - entoc(i,j) )
             else if ( k.eq.2 ) then
              qotent(i,j) = fohfac(2)*entoc(i,j)
             else
              qotent(i,j) = 0.0d0
            endif
*           Add bottom drag contribution
            if ( k.eq.nlo ) then
              qotent(i,j) = qotent(i,j) - bdrfac*del2p(i,j)
            endif
            dqdt(i,j) =  qotjac(i,j) + qt2dif(i,j)
     &                 + qt4dif(i,j) + qotent(i,j)
          enddo
*         Eastern boundary (i=nxpo)
#  ifdef cyclic_ocean
*         Zonally cyclic ocean
          qt2dif(nxpo,j) = qt2dif(1,j)
          qt4dif(nxpo,j) = qt4dif(1,j)
          qotjac(nxpo,j) = qotjac(1,j)
          qotent(nxpo,j) = qotent(1,j)
          dqdt(nxpo,j)   = dqdt(1,j)
#  else
*         Finite box ocean
          qt2dif(nxpo,j) = 0.0d0
          qt4dif(nxpo,j) = 0.0d0
          qotjac(nxpo,j) = 0.0d0
          qotent(nxpo,j) = 0.0d0
          dqdt(nxpo,j)   = rdto*( qo(nxpo,j,k) - qom(nxpo,j,k) )
#  endif /* cyclic_ocean */
        enddo
!$OMP   END DO

*       Add zonal boundary contributions (to dqdt only), by differencing
*       qo in time. This is slightly displaced in time, but only by dto/2
!$OMP   DO SCHEDULE (STATIC)
        do i=1,nxpo
          dqdt(i,  1 ) = rdto*( qo(i,  1 ,k) - qom(i,  1 ,k) )
          dqdt(i,nypo) = rdto*( qo(i,nypo,k) - qom(i,nypo,k) )
        enddo
!$OMP   END DO

!!      Subsample the dqdt array for the current layer
!!      into contiguous values in the wrk vector
!$OMP   DO SCHEDULE (STATIC)
        do j=1,jpwk
           do i=1,ipwk
              wrk(i+ipwk*(j-1)) = dqdt(1+(i-1)*nsko,1+(j-1)*nsko)
           enddo
        enddo
!$OMP   END DO
!!      Output the subsampled values in wrk to the netCDF file
!$OMP   SINGLE
        ncstat = nf_put_vara_double(qocncid, dqdt_id,
     &                              start4, count4, wrk)
        if ( ncstat.ne.NF_NOERR ) call handle_err (ncstat, subnam)
!$OMP   END SINGLE

!!      Subsample the qotjac array for the current layer
!!      into contiguous values in the wrk vector
!$OMP   DO SCHEDULE (STATIC)
        do j=1,jpwk
           do i=1,ipwk
              wrk(i+ipwk*(j-1)) = qotjac(1+(i-1)*nsko,1+(j-1)*nsko)
           enddo
        enddo
!$OMP   END DO
!!      Output the subsampled values in wrk to the netCDF file
!$OMP   SINGLE
        ncstat = nf_put_vara_double(qocncid, qotjac_id,
     &                              start4, count4, wrk)
        if ( ncstat.ne.NF_NOERR ) call handle_err (ncstat, subnam)
!$OMP   END SINGLE

!!      Subsample the qt2dif array for the current layer
!!      into contiguous values in the wrk vector
        if ( qt2dif_out ) then
!$OMP     DO SCHEDULE (STATIC)
          do j=1,jpwk
             do i=1,ipwk
                wrk(i+ipwk*(j-1)) = qt2dif(1+(i-1)*nsko,1+(j-1)*nsko)
             enddo
          enddo
!$OMP     END DO
!!        Output the subsampled values in wrk to the netCDF file
!$OMP     SINGLE
          ncstat = nf_put_vara_double(qocncid, qt2dif_id,
     &                                start4, count4, wrk)
          if ( ncstat.ne.NF_NOERR ) call handle_err (ncstat, subnam)
!$OMP     END SINGLE
        endif

!!      Subsample the qt4dif array for the current layer
!!      into contiguous values in the wrk vector
!$OMP   DO SCHEDULE (STATIC)
        do j=1,jpwk
           do i=1,ipwk
              wrk(i+ipwk*(j-1)) = qt4dif(1+(i-1)*nsko,1+(j-1)*nsko)
           enddo
        enddo
!$OMP   END DO
!!      Output the subsampled values in wrk to the netCDF file
!$OMP   SINGLE
        ncstat = nf_put_vara_double(qocncid, qt4dif_id,
     &                              start4, count4, wrk)
        if ( ncstat.ne.NF_NOERR ) call handle_err (ncstat, subnam)
!$OMP   END SINGLE

!!      Subsample the qotent array for the current layer
!!      into contiguous values in the wrk vector
!$OMP   DO SCHEDULE (STATIC)
        do j=1,jpwk
           do i=1,ipwk
              wrk(i+ipwk*(j-1)) = qotent(1+(i-1)*nsko,1+(j-1)*nsko)
           enddo
        enddo
!$OMP   END DO
!!      Output the subsampled values in wrk to the netCDF file
!$OMP   SINGLE
        ncstat = nf_put_vara_double(qocncid, qotent_id,
     &                              start4, count4, wrk)
        if ( ncstat.ne.NF_NOERR ) call handle_err (ncstat, subnam)
!$OMP   END SINGLE NOWAIT

!$OMP   END PARALLEL

      enddo

#    endif /* use_netcdf */
#  endif /* qoc_diag */
#endif /* not atmos_only */

      END SUBROUTINE qocdiag_out

c***********************************************************************
*
      END MODULE qocdiag
*
c***********************************************************************
