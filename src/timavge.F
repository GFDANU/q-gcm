c***********************************************************************
c     Q-GCM Version 1.5.0 : last modified 05/12/2013
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
      MODULE timavge

*     Contains subprograms for computing and outputting
*     (in netCDF format) various quantities (mainly forcing
*      fields) time averaged over the course of a run.

*     Modules
      USE parameters

      IMPLICIT NONE

      PRIVATE

      PUBLIC :: tavini, tavatm, tavocn, tavout

      integer, PUBLIC, SAVE :: tavncid

#ifndef ocean_only
*     Storage for time averages of atmospheric fields (mainly
*     mixed layer and forcing fields) during the course of a run

      integer, PRIVATE, SAVE :: nsumat
      double precision, PRIVATE :: txatav(nxpa,nypa),tyatav(nxpa,nypa),
     &                 wtatav(nxta,nyta),
     &                 fmatav(nxta,nyta),astav(nxta,nyta),
     &                 patav(nxpa,nypa,nla),qatav(nxpa,nypa,nla),
     &                 uufa(nxpa,nyta),tufa(nxpa,nyta),utufa(nxpa,nyta),
     &                 vvfa(nxta,nypa),tvfa(nxta,nypa),vtvfa(nxta,nypa)

*     Quantities in the arrays during the run are the
*     accumulated sums of the quantities to be averaged,
*     which are transferred in by subroutine tavatm.
*     The arrays are initialised by subroutine tavini, and the
*     values averaged and the results output by subroutine tavout
*
*     txatav, tyatav are dynamic stress components          (m^2 s^-2)
*     wtatav is the Ekman velocity at atmospheric T points  (m s^-1)
*     wtatav is the Ekman velocity at atmospheric T points  (m s^-1)
*     fmatav is the net atmos. mixed layer forcing at T pts (W m^-2)
*     astav  is the atmos. mix. layer temperature (anomaly) (K)
*     patav  is the dynamic pressure in each layer          (m^2 s^-1)
*     qatav  is the vorticity in each layer                 (s^-1)
*     nsumat is the number of contributions
*     to the running totals for the atmosphere
*
*     (replaces the former timavat.cmn)
#endif /* not ocean_only */

#ifndef atmos_only
*     Storage for time averages of oceanic fields (mainly
*     mixed layer and forcing fields) during the course of a run.

      integer, PRIVATE, SAVE :: nsumoc
      double precision, PRIVATE :: txocav(nxpo,nypo),tyocav(nxpo,nypo),
     &                 wpocav(nxpo,nypo),wtocav(nxto,nyto),
     &                 fmocav(nxto,nyto),sstav(nxto,nyto),
     &                 pocav(nxpo,nypo,nlo),qocav(nxpo,nypo,nlo),
     &                 uufo(nxpo,nyto),tufo(nxpo,nyto),utufo(nxpo,nyto),
     &                 vvfo(nxto,nypo),tvfo(nxto,nypo),vtvfo(nxto,nypo)
*
*     Quantities in the arrays during the run are the
*     accumulated sums of the quantities to be averaged,
*     which are transferred in by subroutine tavocn.
*     The arrays are initialised by subroutine tavini, and the
*     values averaged and the results output by subroutine tavout
*
*     txocav, tyocav are dynamic stress components          (m^2 s^-2)
*     woocav is the Ekman velocity at oceanic     p points  (m s^-1)
*     wtocav is the Ekman velocity at oceanic     T points  (m s^-1)
*     fmocav is the net ocean  mixed layer forcing at T pts (W m^-2)
*     sstav  is the ocean  mix. layer temperature (anomaly) (K)
*     pocav  is the dynamic pressure in each layer          (m^2 s^-1)
*     qocav  is the vorticity in each layer                 (s^-1)
*     nsumoc is the number of contributions
*     to the running totals for the ocean
*
*     (replaces the former timavoc.cmn)
#endif /* not atmos_only */

      CONTAINS

c***********************************************************************
c
      SUBROUTINE tavini

*     Initialises the storage used for computing time averages
*     of various quantities (mainly forcing fields) during a run.
*     As we use a simple sum of quantities, initialise everything
*     to zero, roughly in the order it occurs in the declaration.

*     Modules
      USE parameters

      IMPLICIT NONE

*     Subroutine arguments
*
*     Local variables
      integer i,j,k

#ifndef ocean_only
*     Atmosphere fields
*     =================
!$OMP PARALLEL DEFAULT (NONE)
!$OMP&         PRIVATE (i,j,k)
!$OMP&         SHARED  (txatav,tyatav,wtatav,fmatav,astav,patav,qatav,
!$OMP&                  uufa,tufa,utufa,vvfa,tvfa,vtvfa)

*     Windstress fields (p-grid fields)
*     ---------------------------------
!$OMP DO SCHEDULE (STATIC)
      do j=1,nypa
        do i=1,nxpa
          txatav(i,j) = 0.0d0
          tyatav(i,j) = 0.0d0
        enddo
      enddo
!$OMP END DO NOWAIT

*     Ekman velocity, forcing fields & mixed
*     layer temperatures (T-grid fields)
*     --------------------------------------
!$OMP DO SCHEDULE (STATIC)
      do j=1,nyta
        do i=1,nxta
          wtatav(i,j) = 0.0d0
          fmatav(i,j) = 0.0d0
          astav(i,j) = 0.0d0
        enddo
      enddo
!$OMP END DO NOWAIT

*     Dynamic pressure and vorticity in QG layers
*     -------------------------------------------
      do k=1,nla
!$OMP   DO SCHEDULE (STATIC)
        do j=1,nypa
          do i=1,nxpa
            patav(i,j,k) = 0.0d0
            qatav(i,j,k) = 0.0d0
          enddo
        enddo
!$OMP   END DO NOWAIT
      enddo

*     Fields for mixed layer temperature advection calculations
*     ---------------------------------------------------------
!$OMP DO SCHEDULE (STATIC)
      do j=1,nyta
        do i=1,nxpa
          uufa (i,j) = 0.0d0
          tufa (i,j) = 0.0d0
          utufa(i,j) = 0.0d0
        enddo
      enddo
!$OMP END DO NOWAIT
!$OMP DO SCHEDULE (STATIC)
      do j=1,nypa
        do i=1,nxta
          vvfa (i,j) = 0.0d0
          tvfa (i,j) = 0.0d0
          vtvfa(i,j) = 0.0d0
        enddo
      enddo
!$OMP END DO NOWAIT

!$OMP END PARALLEL
#endif /* not ocean_only */

#ifndef atmos_only
*     Ocean fields
*     ============
!$OMP PARALLEL DEFAULT (NONE)
!$OMP&         PRIVATE (i,j,k)
!$OMP&         SHARED  (txocav,tyocav,wtocav,fmocav,sstav,pocav,qocav,
!$OMP&                  wpocav,uufo,tufo,utufo,vvfo,tvfo,vtvfo)

*     Windstress fields & wekpo (p-grid fields)
*     -----------------------------------------
!$OMP DO SCHEDULE (STATIC)
      do j=1,nypo
        do i=1,nxpo
          txocav(i,j) = 0.0d0
          tyocav(i,j) = 0.0d0
          wpocav(i,j) = 0.0d0
        enddo
      enddo
!$OMP END DO NOWAIT

*     Ekman velocity, forcing fields & mixed
*     layer temperatures (T-grid fields)
*     --------------------------------------
!$OMP DO SCHEDULE (STATIC)
      do j=1,nyto
        do i=1,nxto
          wtocav(i,j) = 0.0d0
          fmocav(i,j) = 0.0d0
          sstav(i,j) = 0.0d0
        enddo
      enddo
!$OMP END DO NOWAIT

*     Dynamic pressure and vorticity in QG layers
*     -------------------------------------------
      do k=1,nlo
!$OMP   DO SCHEDULE (STATIC)
        do j=1,nypo
          do i=1,nxpo
            pocav(i,j,k) = 0.0d0
            qocav(i,j,k) = 0.0d0
          enddo
        enddo
!$OMP   END DO NOWAIT
      enddo

*     Fields for mixed layer temperature advection calculations
*     ---------------------------------------------------------
!$OMP DO SCHEDULE (STATIC)
      do j=1,nyto
        do i=1,nxpo
          uufo (i,j) = 0.0d0
          tufo (i,j) = 0.0d0
          utufo(i,j) = 0.0d0
        enddo
      enddo
!$OMP END DO NOWAIT
!$OMP DO SCHEDULE (STATIC)
      do j=1,nypo
        do i=1,nxto
          vvfo (i,j) = 0.0d0
          tvfo (i,j) = 0.0d0
          vtvfo(i,j) = 0.0d0
        enddo
      enddo
!$OMP END DO NOWAIT

!$OMP END PARALLEL
#endif /* not atmos_only */

*     Counters
*     --------
#ifndef ocean_only
      nsumat = 0
#endif
#ifndef atmos_only
      nsumoc = 0
#endif

      END SUBROUTINE tavini
c
c***********************************************************************
c
      SUBROUTINE tavatm

*     Accumulates various atmospheric quantities (mainly forcing
*     fields) during a run, for later computation of time averages.

*     Modules
      USE parameters
#ifndef ocean_only
      USE atconst
      USE atstate
#endif
      USE intrfac, ONLY : ast, tauxa, tauya, fnetat, hmat

      IMPLICIT NONE

*     Subroutine arguments
*
#ifndef ocean_only
*
*     Local variables
      integer i,j,k
      double precision rhf0hm,uuf,tuf(nxpa),vvf,tvf,vtvf

      rhf0hm = 0.5d0/(fnot*hmat)

!$OMP PARALLEL DEFAULT (NONE)
!$OMP&         PRIVATE (i,j,k,uuf,tuf,vvf,tvf,vtvf)
!$OMP&         SHARED  (txatav,tauxa,tyatav,tauya,wtatav,wekta,
!$OMP&                  fmatav,fnetat,astav,ast,rdxaf0,pa,qa,rhf0hm,
!$OMP&                  uufa,tufa,utufa,vvfa,tvfa,vtvfa,patav,qatav)

*     Windstress fields
*     -----------------
!$OMP DO SCHEDULE (STATIC)
      do j=1,nypa
        do i=1,nxpa
          txatav(i,j) = txatav(i,j) + tauxa(i,j)
          tyatav(i,j) = tyatav(i,j) + tauya(i,j)
        enddo
      enddo
!$OMP END DO NOWAIT

*     Ekman velocity, forcing fields & mixed
*     layer temperatures (T-grid fields)
*     --------------------------------------
!$OMP DO SCHEDULE (STATIC)
      do j=1,nyta
        do i=1,nxta
          wtatav(i,j) = wtatav(i,j) + wekta(i,j)
          fmatav(i,j) = fmatav(i,j) + fnetat(i,j)
          astav(i,j) = astav(i,j) + ast(i,j)
        enddo
      enddo
!$OMP END DO NOWAIT

*     Fields for mixed layer temperature advection calculations
*     ---------------------------------------------------------
*
*     C-grid advection, second order accurate,
*     borrowed from subroutine amladv in amlsubs.f

*     Zonal advection
*     ---------------
!$OMP DO SCHEDULE (STATIC)
      do j=1,nyta
*       Western boundary (periodic)
        tuf( 1 ) = 0.5d0*( ast(1,j) + ast(nxta,j) )
*       Inner points + Eastern boundary
        do i=2,nxpa-1
          tuf(i) = 0.5d0*( ast(i,j) + ast(i-1,j) )
        enddo
        tuf(nxpa) = 0.5d0*( ast(1,j) + ast(nxta,j) )
        do i=1,nxpa
          uuf = -rdxaf0*( pa(i,j+1,1)-pa(i,j,1) )
     &          -rhf0hm*(tauya(i,j+1)+tauya(i,j))
          uufa (i,j) = uufa (i,j) + uuf
          tufa (i,j) = tufa (i,j) + tuf(i)
          utufa(i,j) = utufa(i,j) + uuf*tuf(i)
        enddo
      enddo
!$OMP END DO NOWAIT

*     Meridional advection
*     --------------------
*     Inner points + meridional boundaries
!$OMP DO SCHEDULE (STATIC)
      do j=2,nyta
        do i=1,nxta
          vvf  =  rdxaf0*( pa(i+1,j,1)-pa(i,j,1) )
     &           +rhf0hm*(tauxa(i+1,j)+tauxa(i,j))
          tvf  = 0.5d0*( ast(i,j) + ast(i,j-1) )
          vtvf = vvf*tvf
          vvfa (i,j) = vvfa (i,j) + vvf
          tvfa (i,j) = tvfa (i,j) + tvf
          vtvfa(i,j) = vtvfa(i,j) + vtvf
        enddo
      enddo
!$OMP END DO NOWAIT

*     Zonal boundaries, including corners
*     -----------------------------------
*     N.B. P constant on boundaries, so Px vanishes.
*     No normal flux of heat, so deem v to be zero.
*     Deem boundary temperature to be that of internal point.
!$OMP DO SCHEDULE (STATIC)
      do i=1,nxta
*       Southern boundary
        vvf  = 0.0d0
        tvf  = ast(i, 1 )
        vtvf = 0.0d0
        vvfa (i, 1 ) = vvfa (i, 1 ) + vvf
        tvfa (i, 1 ) = tvfa (i, 1 ) + tvf
        vtvfa(i, 1 ) = vtvfa(i, 1 ) + vtvf
*       Northern boundary
        vvf  = 0.0d0
        tvf  = ast(i,nypa-1)
        vtvf = 0.0d0
        vvfa (i,nypa) = vvfa (i,nypa) + vvf
        tvfa (i,nypa) = tvfa (i,nypa) + tvf
        vtvfa(i,nypa) = vtvfa(i,nypa) + vtvf
      enddo
!$OMP END DO NOWAIT

*     Dynamic pressure and vorticity in QG layers
*     -------------------------------------------
      do k=1,nla
!$OMP   DO SCHEDULE (STATIC)
        do j=1,nypa
          do i=1,nxpa
            patav(i,j,k) = patav(i,j,k) + pa(i,j,k)
            qatav(i,j,k) = qatav(i,j,k) + qa(i,j,k)
          enddo
        enddo
!$OMP   END DO NOWAIT
      enddo

!$OMP END PARALLEL

*     Update counter
*     --------------
      nsumat = nsumat + 1
#endif

      END SUBROUTINE tavatm
c
c***********************************************************************
c
      SUBROUTINE tavocn

*     Accumulates various oceanic quantities (mainly forcing
*     fields) during a run, for later computation of time averages.

*     Modules
      USE parameters
#ifndef atmos_only
      USE occonst
      USE ocstate
#endif
      USE intrfac, ONLY : sst, tauxo, tauyo, fnetoc, hmoc, tsbdy, tnbdy

      IMPLICIT NONE

*     Subroutine arguments
*
#ifndef atmos_only
*
*     Local variables
      integer i,j,k
      double precision uvgfac,rhf0hm,uuf(nxpo),tuf(nxpo),utuf(nxpo),
     &                 vvf,tvf,vtvf
*
      uvgfac = ycexp*rdxof0
      rhf0hm = 0.5d0/(fnot*hmoc)

!$OMP PARALLEL DEFAULT (NONE)
!$OMP&         PRIVATE (i,j,k,uuf,tuf,utuf,vvf,tvf,vtvf)
!$OMP&         SHARED  (txocav,tauxo,tyocav,tauyo,wpocav,wekpo,wtocav,
!$OMP&                  wekto,fmocav,fnetoc,sstav,sst,uvgfac,po,qo,
!$OMP&                  rhf0hm,uufo,tufo,utufo,tsbdy,tnbdy,pocav,qocav,
!$OMP&                  vvfo,tvfo,vtvfo)
*
*     Windstress fields & wekpo (p-grid fields)
*     -----------------------------------------
!$OMP DO SCHEDULE (STATIC)
      do j=1,nypo
        do i=1,nxpo
          txocav(i,j) = txocav(i,j) + tauxo(i,j)
          tyocav(i,j) = tyocav(i,j) + tauyo(i,j)
          wpocav(i,j) = wpocav(i,j) + wekpo(i,j)
        enddo
      enddo
!$OMP END DO NOWAIT

*     Ekman velocity, forcing fields & mixed
*     layer temperatures (T-grid fields)
*     --------------------------------------
!$OMP DO SCHEDULE (STATIC)
      do j=1,nyto
        do i=1,nxto
          wtocav(i,j) = wtocav(i,j) + wekto(i,j)
          fmocav(i,j) = fmocav(i,j) + fnetoc(i,j)
          sstav(i,j) = sstav(i,j) + sst(i,j)
        enddo
      enddo
!$OMP END DO NOWAIT

*     Fields for mixed layer temperature advection calculations
*     ---------------------------------------------------------
*     C-grid advection, second order accurate,
*     borrowed from subroutine omladv in omlsubs.f

*     Zonal advection
*     ---------------
*     In finite box case, no normal flux for temperature, so deem u to
*     be zero. Deem boundary temperature to be that of internal point.

!$OMP DO SCHEDULE (STATIC)
      do j=1,nyto
*       Western boundary
#  ifdef cyclic_ocean
*       Zonally cyclic ocean
        uuf (1) = -uvgfac*( po(1,j+1,1)-po(1,j,1) )
     &            +rhf0hm*(tauyo(1,j+1)+tauyo(1,j))
        tuf (1) = 0.5d0*( sst(1,j) + sst(nxto,j) )
        utuf(1) = uuf(1)*tuf(1)
#  else
*       Finite box ocean (no normal heat flux)
        uuf ( 1 ) = 0.0d0
        tuf ( 1 ) = sst(1,j)
        utuf( 1 ) = 0.0d0
#  endif
*       Inner points.
        do i=2,nxpo-1
          uuf (i) = -uvgfac*( po(i,j+1,1)-po(i,j,1) )
     &              +rhf0hm*(tauyo(i,j+1)+tauyo(i,j))
          tuf (i) = 0.5d0*( sst(i,j) + sst(i-1,j) )
          utuf(i) = uuf(i)*tuf(i)
        enddo
*       Eastern boundary
#  ifdef cyclic_ocean
*       Zonally cyclic ocean
        uuf (nxpo) = uuf (1)
        tuf (nxpo) = tuf (1)
        utuf(nxpo) = utuf(1)
#  else
*       Finite box ocean (no normal heat flux)
        uuf (nxpo) = 0.0d0
        tuf (nxpo) = sst(nxto,j)
        utuf(nxpo) = 0.0d0
#  endif
        do i=1,nxpo
          uufo (i,j) = uufo (i,j) + uuf (i)
          tufo (i,j) = tufo (i,j) + tuf (i)
          utufo(i,j) = utufo(i,j) + utuf(i)
        enddo
      enddo
!$OMP END DO NOWAIT

*     Meridional advection
*     --------------------
*     Points excluding zonal boundaries
!$OMP DO SCHEDULE (STATIC)
      do j=2,nyto
        do i=1,nxto
          vvf  =  uvgfac*( po(i+1,j,1)-po(i,j,1) )
     &           -rhf0hm*(tauxo(i+1,j)+tauxo(i,j))
          tvf  = 0.5d0*( sst(i,j) + sst(i,j-1) )
          vtvf = vvf*tvf
          vvfo (i,j) = vvfo (i,j) + vvf
          tvfo (i,j) = tvfo (i,j) + tvf
          vtvfo(i,j) = vtvfo(i,j) + vtvf
        enddo
      enddo
!$OMP END DO NOWAIT

*     Zonal boundaries including corners; choice of BCs
*     -------------------------------------------------
!$OMP DO SCHEDULE (STATIC)
      do i=1,nxto

*       Southern boundary
#  ifdef sb_hflux
*       Advection consistent with an outflow across the
*       southern boundary equal to the Ekman transport,
*       carrying fluid of a specified temperature tsbdy.
*       p contribution to vm vanishes because p is uniform along bdy.
        vvf  = -rhf0hm*( tauxo(i+1,1) + tauxo(i,1) )
        tvf  = 0.5d0*( sst(i,1) + tsbdy )
        vtvf = vvf*tvf
#  else
*       Advection (no normal mass flux => vm = 0)
        vvf  = 0.0d0
        tvf  = sst(i,1)
        vtvf = 0.0d0
#  endif
        vvfo (i, 1 ) = vvfo (i, 1 ) + vvf
        tvfo (i, 1 ) = tvfo (i, 1 ) + tvf
        vtvfo(i, 1 ) = vtvfo(i, 1 ) + vtvf

*       Northern boundary
#  ifdef nb_hflux
*       Advection consistent with an outflow across the
*       northern boundary equal to the Ekman transport,
*       carrying fluid of a specified temperature tnbdy.
*       p contribution to vp vanishes because p is uniform along bdy.
        vvf  = -rhf0hm*( tauxo(i+1,nyto+1) + tauxo(i,nyto+1) )
        tvf  = 0.5d0*( sst(i,nyto) + tnbdy )
        vtvf = vvf*tvf
#  else
*       Advection (no normal mass flux => vp = 0)
        vvf  = 0.0d0
        tvf  = sst(i,nyto)
        vtvf = 0.0d0
#  endif
        vvfo (i,nypo) = vvfo (i,nypo) + vvf
        tvfo (i,nypo) = tvfo (i,nypo) + tvf
        vtvfo(i,nypo) = vtvfo(i,nypo) + vtvf

      enddo
!$OMP END DO NOWAIT

*     Dynamic pressure and vorticity in QG layers
*     -------------------------------------------
      do k=1,nlo
!$OMP   DO SCHEDULE (STATIC)
        do j=1,nypo
          do i=1,nxpo
            pocav(i,j,k) = pocav(i,j,k) + po(i,j,k)
            qocav(i,j,k) = qocav(i,j,k) + qo(i,j,k)
          enddo
        enddo
!$OMP   END DO NOWAIT
      enddo

!$OMP END PARALLEL

*     Update counter
*     --------------
      nsumoc = nsumoc + 1
#endif

      END SUBROUTINE tavocn
c
c***********************************************************************
c
      SUBROUTINE tavout

*     Computes time averages of various quantities (mainly forcing
*     fields) at the end of a run, and outputs a netCDF dump of them.

*     Modules
      USE parameters
      USE atconst
      USE occonst
      USE nc_subs, ONLY : handle_err

      IMPLICIT NONE

#ifdef use_netcdf
      INCLUDE 'netcdf.inc'
#endif /* use_netcdf */
*
*     Subroutine arguments
*
*     Local parameters
      character (len=*), parameter :: subnam = 'tavout'

*     Local variables
      integer i,j,k
#ifndef ocean_only
      double precision rnsat,uptpat(nxpa,nyta),vptpat(nxta,nypa)
#endif
#ifndef atmos_only
      double precision rnsoc,uptpoc(nxpo,nyto),vptpoc(nxto,nypo)
#endif

#ifdef use_netcdf
      integer status,xdims,dims(2),pdims(3),start(2),
     &        count(2),startp(3),countp(3)
#  ifndef ocean_only
      integer xapdim,yapdim,ladim,xatdim,yatdim,xap_id,xat_id,
     &        yap_id,yat_id,la_id,ast_id,wekta_id,fmat_id,
     &        txa_id,tya_id,pa_id,qa_id,utat_id,vtat_id
      double precision xxa(nxpa),yya(nypa),tmpa(nla)
#  endif
#  ifndef atmos_only
      integer xopdim,yopdim,lodim,xotdim,yotdim,xop_id,xot_id,
     &        yop_id,yot_id,lo_id,sst_id,wekto_id,fmoc_id,txo_id,
     &        tyo_id,wekpo_id,po_id,qo_id,utoc_id,vtoc_id
      double precision xxo(nxpo),yyo(nypo),tmpo(nlo)
#  endif
#endif /* use_netcdf */

#ifndef ocean_only
      if ( nsumat.eq.0 ) then
        rnsat = 0.0d0
       else
        rnsat = 1.0d0/dble(nsumat)
      endif
#endif
#ifndef atmos_only
      if ( nsumoc.eq.0 ) then
        rnsoc = 0.0d0
       else
        rnsoc = 1.0d0/dble(nsumoc)
      endif
#endif

#ifndef ocean_only
*     Atmosphere fields
*     =================
!$OMP PARALLEL DEFAULT (NONE)
!$OMP&         PRIVATE (i,j,k)
!$OMP&         SHARED  (rnsat,txatav,tyatav,wtatav,fmatav,astav,
!$OMP&                  patav,qatav,uufa,tufa,utufa,uptpat,
!$OMP&                  vvfa,tvfa,vtvfa,vptpat)

*     Windstress fields (p-grid fields)
*     ---------------------------------
!$OMP DO SCHEDULE (STATIC)
      do j=1,nypa
        do i=1,nxpa
          txatav(i,j) = rnsat*txatav(i,j)
          tyatav(i,j) = rnsat*tyatav(i,j)
        enddo
      enddo
!$OMP END DO NOWAIT

*     Ekman velocity, forcing fields & mixed
*     layer temperatures (T-grid fields)
*     --------------------------------------
!$OMP DO SCHEDULE (STATIC)
      do j=1,nyta
        do i=1,nxta
          wtatav(i,j) = rnsat*wtatav(i,j)
          fmatav(i,j) = rnsat*fmatav(i,j)
          astav(i,j) = rnsat*astav(i,j)
        enddo
      enddo
!$OMP END DO NOWAIT

*     Dynamic pressure and vorticity in QG layers
*     -------------------------------------------
      do k=1,nla
!$OMP   DO SCHEDULE (STATIC)
        do j=1,nypa
          do i=1,nxpa
            patav(i,j,k) = rnsat*patav(i,j,k)
            qatav(i,j,k) = rnsat*qatav(i,j,k)
          enddo
        enddo
!$OMP   END DO NOWAIT
      enddo

*     Fields for mixed layer temperature flux calculations
*     ----------------------------------------------------
*     Compute means, and derive eddy fluxes
!$OMP DO SCHEDULE (STATIC)
      do j=1,nyta
        do i=1,nxpa
          uufa (i,j) = rnsat*uufa (i,j)
          tufa (i,j) = rnsat*tufa (i,j)
          utufa(i,j) = rnsat*utufa(i,j)
          uptpat(i,j) = utufa(i,j) - uufa(i,j)*tufa(i,j)
        enddo
      enddo
!$OMP END DO NOWAIT
!$OMP DO SCHEDULE (STATIC)
      do j=1,nypa
        do i=1,nxta
          vvfa (i,j) = rnsat*vvfa (i,j)
          tvfa (i,j) = rnsat*tvfa (i,j)
          vtvfa(i,j) = rnsat*vtvfa(i,j)
          vptpat(i,j) = vtvfa(i,j) - vvfa(i,j)*tvfa(i,j)
        enddo
      enddo
!$OMP END DO NOWAIT

!$OMP END PARALLEL
#endif

#ifndef atmos_only
*     Ocean fields
*     ============
!$OMP PARALLEL DEFAULT (NONE)
!$OMP&         PRIVATE (i,j,k)
!$OMP&         SHARED  (rnsoc,txocav,tyocav,wpocav,wtocav,fmocav,
!$OMP&                  sstav,pocav,qocav,uufo,tufo,utufo,uptpoc,
!$OMP&                  vvfo,tvfo,vtvfo,vptpoc)

*     Windstress fields & wekpo (p-grid fields)
*     -----------------------------------------
!$OMP DO SCHEDULE (STATIC)
      do j=1,nypo
        do i=1,nxpo
          txocav(i,j) = rnsoc*txocav(i,j)
          tyocav(i,j) = rnsoc*tyocav(i,j)
          wpocav(i,j) = rnsoc*wpocav(i,j)
        enddo
      enddo
!$OMP END DO NOWAIT

*     Ekman velocity, forcing fields & mixed
*     layer temperatures (T-grid fields)
*     --------------------------------------
!$OMP DO SCHEDULE (STATIC)
      do j=1,nyto
        do i=1,nxto
          wtocav(i,j) = rnsoc*wtocav(i,j)
          fmocav(i,j) = rnsoc*fmocav(i,j)
          sstav(i,j) = rnsoc*sstav(i,j)
        enddo
      enddo
!$OMP END DO NOWAIT

*     Dynamic pressure and vorticity in QG layers
*     -------------------------------------------
      do k=1,nlo
!$OMP   DO SCHEDULE (STATIC)
        do j=1,nypo
          do i=1,nxpo
            pocav(i,j,k) = rnsoc*pocav(i,j,k)
            qocav(i,j,k) = rnsoc*qocav(i,j,k)
          enddo
        enddo
!$OMP   END DO NOWAIT
      enddo

*     Fields for mixed layer temperature flux calculations
*     ----------------------------------------------------
*     Compute means, and derive eddy fluxes
!$OMP DO SCHEDULE (STATIC)
      do j=1,nyto
        do i=1,nxpo
          uufo (i,j) = rnsoc*uufo (i,j)
          tufo (i,j) = rnsoc*tufo (i,j)
          utufo(i,j) = rnsoc*utufo(i,j)
          uptpoc(i,j) = utufo(i,j) - uufo(i,j)*tufo(i,j)
        enddo
      enddo
!$OMP END DO NOWAIT
!$OMP DO SCHEDULE (STATIC)
      do j=1,nypo
        do i=1,nxto
          vvfo (i,j) = rnsoc*vvfo (i,j)
          tvfo (i,j) = rnsoc*tvfo (i,j)
          vtvfo(i,j) = rnsoc*vtvfo(i,j)
          vptpoc(i,j) = vtvfo(i,j) - vvfo(i,j)*tvfo(i,j)
        enddo
      enddo
!$OMP END DO NOWAIT

!$OMP END PARALLEL
#endif

#ifdef use_netcdf
!!    Initialise netCDF file
#  ifndef atmos_only
      status = nf_def_dim(tavncid, 'xpo', nxpo, xopdim)
      if ( status.ne.NF_NOERR ) call handle_err (status, subnam)
      status = nf_def_dim(tavncid, 'ypo', nypo, yopdim)
      if ( status.ne.NF_NOERR ) call handle_err (status, subnam)
      status = nf_def_dim(tavncid, 'zo', nlo, lodim)
      if ( status.ne.NF_NOERR ) call handle_err (status, subnam)

      status = nf_def_dim(tavncid, 'xto', nxto, xotdim)
      if ( status.ne.NF_NOERR ) call handle_err (status, subnam)
      status = nf_def_dim(tavncid, 'yto', nyto, yotdim)
      if ( status.ne.NF_NOERR ) call handle_err (status, subnam)
#  endif /* not atmos_only */

#  ifndef ocean_only
      status = nf_def_dim(tavncid, 'xpa', nxpa, xapdim)
      if ( status.ne.NF_NOERR ) call handle_err (status, subnam)
      status = nf_def_dim(tavncid, 'ypa', nypa, yapdim)
      if ( status.ne.NF_NOERR ) call handle_err (status, subnam)
      status = nf_def_dim(tavncid, 'za', nla, ladim)
      if ( status.ne.NF_NOERR ) call handle_err (status, subnam)

      status = nf_def_dim(tavncid, 'xta', nxta, xatdim)
      if ( status.ne.NF_NOERR ) call handle_err (status, subnam)
      status = nf_def_dim(tavncid, 'yta', nyta, yatdim)
      if ( status.ne.NF_NOERR ) call handle_err (status, subnam)
#  endif /* not ocean_only */

#  ifndef atmos_only
!!  Define a 1-dimensional variable called
!!  'xpo' which stores x gridpoints (p-grid)
      xdims = xopdim
      status = nf_def_var(tavncid, 'xpo', NF_DOUBLE, 1, xdims, xop_id)
      if ( status.ne.NF_NOERR ) call handle_err (status, subnam)
      status = nf_put_att_text(tavncid,xop_id,'units',2,'km')
      if ( status.ne.NF_NOERR ) call handle_err (status, subnam)
      status = nf_put_att_text(tavncid,xop_id,'long_name',
     &            21, 'Ocean X axis (p-grid)')

!!  Define a 1-dimensional variable called
!!  'xto' which stores x gridpoints (T-grid)
      xdims = xotdim
      status = nf_def_var(tavncid, 'xto', NF_DOUBLE, 1, xdims, xot_id)
      if ( status.ne.NF_NOERR ) call handle_err (status, subnam)
      status = nf_put_att_text(tavncid,xot_id,'units',2,'km')
      if ( status.ne.NF_NOERR ) call handle_err (status, subnam)
      status = nf_put_att_text(tavncid,xot_id,'long_name',
     &            21, 'Ocean X axis (T-grid)')
      if ( status.ne.NF_NOERR ) call handle_err (status, subnam)

!!  Define a 1-dimensional variable called
!!  'ypo' which stores y gridpoints (p-grid)
      xdims = yopdim
      status = nf_def_var(tavncid, 'ypo', NF_DOUBLE, 1, xdims, yop_id)
      if ( status.ne.NF_NOERR ) call handle_err (status, subnam)
      status = nf_put_att_text(tavncid,yop_id,'units',2,'km')
      if ( status.ne.NF_NOERR ) call handle_err (status, subnam)
      status = nf_put_att_text(tavncid,yop_id,'long_name',
     &            21, 'Ocean Y axis (p-grid)')
      if ( status.ne.NF_NOERR ) call handle_err (status, subnam)

!!  Define a 1-dimensional variable called
!!  'yto' which stores y gridpoints (T-grid)
      xdims = yotdim
      status = nf_def_var(tavncid, 'yto', NF_DOUBLE, 1, xdims, yot_id)
      if ( status.ne.NF_NOERR ) call handle_err (status, subnam)
      status = nf_put_att_text(tavncid,yot_id,'units',2,'km')
      if ( status.ne.NF_NOERR ) call handle_err (status, subnam)
      status = nf_put_att_text(tavncid,yot_id,'long_name',
     &            21, 'Ocean Y axis (T-grid)')
      if ( status.ne.NF_NOERR ) call handle_err (status, subnam)

!!  Define a 1-dimensional variable called 'zo'
!!  to store layer thicknesses
      xdims = lodim
      status = nf_def_var(tavncid, 'zo', NF_DOUBLE, 1, xdims, lo_id)
      if ( status.ne.NF_NOERR ) call handle_err (status, subnam)
      status = nf_put_att_text(tavncid,lo_id,'units',2,'km')
      if ( status.ne.NF_NOERR ) call handle_err (status, subnam)
      status = nf_put_att_text(tavncid,lo_id,'long_name',
     &            16, 'Ocean depth axis')
      if ( status.ne.NF_NOERR ) call handle_err (status, subnam)
#  endif /* not atmos_only */

#  ifndef ocean_only
!!  Define a 1-dimensional variable called
!!  'xpa' which stores x gridpoints (p-grid)
      xdims = xapdim
      status = nf_def_var(tavncid, 'xpa', NF_DOUBLE, 1, xdims, xap_id)
      if ( status.ne.NF_NOERR ) call handle_err (status, subnam)
      status = nf_put_att_text(tavncid,xap_id,'units',2, 'km')
      if ( status.ne.NF_NOERR ) call handle_err (status, subnam)
      status = nf_put_att_text(tavncid,xap_id,'long_name',
     &            26, 'Atmosphere X axis (p-grid)')
      if ( status.ne.NF_NOERR ) call handle_err (status, subnam)

!!  Define a 1-dimensional variable called
!!  'xta' which stores x gridpoints (T-grid)
      xdims = xatdim
      status = nf_def_var(tavncid, 'xta', NF_DOUBLE, 1, xdims, xat_id)
      if ( status.ne.NF_NOERR ) call handle_err (status, subnam)
      status = nf_put_att_text(tavncid,xat_id,'units',2,'km')
      if ( status.ne.NF_NOERR ) call handle_err (status, subnam)
      status = nf_put_att_text(tavncid,xat_id,'long_name',
     &            26, 'Atmosphere X axis (T-grid)')
      if ( status.ne.NF_NOERR ) call handle_err (status, subnam)

!!  Define a 1-dimensional variable called
!!  'ypa' which stores y gridpoints (p-grid)
      xdims = yapdim
      status = nf_def_var(tavncid, 'ypa', NF_DOUBLE, 1, xdims, yap_id)
      if ( status.ne.NF_NOERR ) call handle_err (status, subnam)
      status = nf_put_att_text(tavncid,yap_id,'units',2,'km')
      if ( status.ne.NF_NOERR ) call handle_err (status, subnam)
      status = nf_put_att_text(tavncid,yap_id,'long_name',
     &            26, 'Atmosphere Y axis (p-grid)')
      if ( status.ne.NF_NOERR ) call handle_err (status, subnam)

!!  Define a 1-dimensional variable called
!!  'yta' which stores y gridpoints (T-grid)
      xdims = yatdim
      status = nf_def_var(tavncid, 'yta', NF_DOUBLE, 1, xdims, yat_id)
      if ( status.ne.NF_NOERR ) call handle_err (status, subnam)
      status = nf_put_att_text(tavncid,yat_id,'units',2,'km')
      if ( status.ne.NF_NOERR ) call handle_err (status, subnam)
      status = nf_put_att_text(tavncid,yat_id,'long_name',
     &            26, 'Atmosphere Y axis (T-grid)')
      if ( status.ne.NF_NOERR ) call handle_err (status, subnam)

!!  Define a 1-dimensional variable called 'za'
!!  to store layer thicknesses
!!  p-grid only
      xdims = ladim
      status = nf_def_var(tavncid, 'za', NF_DOUBLE, 1, xdims, la_id)
      if ( status.ne.NF_NOERR ) call handle_err (status, subnam)
      status = nf_put_att_text(tavncid,la_id,'units',2,'km')
      if ( status.ne.NF_NOERR ) call handle_err (status, subnam)
      status = nf_put_att_text(tavncid,la_id,'long_name',
     &            22, 'Atmosphere height axis')
      if ( status.ne.NF_NOERR ) call handle_err (status, subnam)
#  endif /* not ocean_only */

#  ifndef atmos_only
!!  Define a 2-dimensional variable called 'sst'
!!  to store ocean surface temperatures
      dims(1) = xotdim
      dims(2) = yotdim
      status = nf_def_var(tavncid, 'sst', NF_DOUBLE, 2, dims, sst_id)
      if ( status.ne.NF_NOERR ) call handle_err (status, subnam)
      status = nf_put_att_text(tavncid,sst_id,'units',1,'K')
      if ( status.ne.NF_NOERR ) call handle_err (status, subnam)
      status = nf_put_att_text(tavncid,sst_id,'long_name',
     &            25, 'Ocean surface temperature')
      if ( status.ne.NF_NOERR ) call handle_err (status, subnam)

!!  Define a 2-dimensional variable called 'wekto'
!!  to store ocean Ekman velocity (on T-grid)
      dims(1) = xotdim
      dims(2) = yotdim
      status = nf_def_var(tavncid, 'wekto', NF_DOUBLE,
     &                    2, dims, wekto_id)
      if ( status.ne.NF_NOERR ) call handle_err (status, subnam)
      status = nf_put_att_text(tavncid,wekto_id,'units',3,'m/s')
      if ( status.ne.NF_NOERR ) call handle_err (status, subnam)
      status = nf_put_att_text(tavncid,wekto_id,'long_name',
     &            29, 'Ocean Ekman velocity (T-grid)')
      if ( status.ne.NF_NOERR ) call handle_err (status, subnam)

!!  Define a 2-dimensional variable called 'fnetoc'
!!  to store ocean forcing
      dims(1) = xotdim
      dims(2) = yotdim
      status = nf_def_var(tavncid, 'fnetoc', NF_DOUBLE,
     &                    2, dims, fmoc_id)
      if ( status.ne.NF_NOERR ) call handle_err (status, subnam)
      status = nf_put_att_text(tavncid,fmoc_id,'units',5,'W/m^2')
      if ( status.ne.NF_NOERR ) call handle_err (status, subnam)
      status = nf_put_att_text(tavncid,fmoc_id,'long_name',
     &            30, 'Ocean surface forcing (T-grid)')
      if ( status.ne.NF_NOERR ) call handle_err (status, subnam)

!!  Define a 2-dimensional variable called 'tauxo'
!!  to store ocean dynamic windstress
      dims(1) = xopdim
      dims(2) = yopdim
      status = nf_def_var(tavncid, 'tauxo', NF_DOUBLE, 2, dims, txo_id)
      if ( status.ne.NF_NOERR ) call handle_err (status, subnam)
      status = nf_put_att_text(tavncid,txo_id,'units',7,'m^2/s^2')
      if ( status.ne.NF_NOERR ) call handle_err (status, subnam)
      status = nf_put_att_text(tavncid,txo_id,'long_name',
     &            31, 'Ocean surface x stress (p-grid)')
      if ( status.ne.NF_NOERR ) call handle_err (status, subnam)

!!  Define a 2-dimensional variable called 'tauyo'
!!  to store ocean dynamic windstress
      dims(1) = xopdim
      dims(2) = yopdim
      status = nf_def_var(tavncid, 'tauyo', NF_DOUBLE, 2, dims, tyo_id)
      if ( status.ne.NF_NOERR ) call handle_err (status, subnam)
      status = nf_put_att_text(tavncid,tyo_id,'units',7,'m^2/s^2')
      if ( status.ne.NF_NOERR ) call handle_err (status, subnam)
      status = nf_put_att_text(tavncid,tyo_id,'long_name',
     &            31, 'Ocean surface y stress (p-grid)')
      if ( status.ne.NF_NOERR ) call handle_err (status, subnam)

!!  Define a 2-dimensional variable called 'wekpo'
!!  to store ocean Ekman velocity (on p-grid)
      dims(1) = xopdim
      dims(2) = yopdim
      status = nf_def_var(tavncid, 'wekpo', NF_DOUBLE,
     &                    2, dims, wekpo_id)
      if ( status.ne.NF_NOERR ) call handle_err (status, subnam)
      status = nf_put_att_text(tavncid,wekpo_id,'units',3,'m/s')
      if ( status.ne.NF_NOERR ) call handle_err (status, subnam)
      status = nf_put_att_text(tavncid,wekpo_id,'long_name',
     &            29, 'Ocean Ekman velocity (p-grid)')
      if ( status.ne.NF_NOERR ) call handle_err (status, subnam)

!!  Define a 3-dimensional variable called 'po'
!!  to store ocean dynamic pressure
      pdims(1) = xopdim
      pdims(2) = yopdim
      pdims(3) = lodim
      status = nf_def_var(tavncid, 'po', NF_DOUBLE, 3, pdims, po_id)
      if ( status.ne.NF_NOERR ) call handle_err (status, subnam)
      status = nf_put_att_text(tavncid,po_id,'units',7,'m^2/s^2')
      if ( status.ne.NF_NOERR ) call handle_err (status, subnam)
      status = nf_put_att_text(tavncid,po_id,'long_name',
     &            22, 'Ocean dynamic pressure')
      if ( status.ne.NF_NOERR ) call handle_err (status, subnam)

!!  Define a 3-dimensional variable called 'qo'
!!  to store ocean vorticity
      pdims(1) = xopdim
      pdims(2) = yopdim
      pdims(3) = lodim
      status = nf_def_var(tavncid, 'qo', NF_DOUBLE, 3, pdims, qo_id)
      if ( status.ne.NF_NOERR ) call handle_err (status, subnam)
      status = nf_put_att_text(tavncid,qo_id,'units',4,'s^-1')
      if ( status.ne.NF_NOERR ) call handle_err (status, subnam)
      status = nf_put_att_text(tavncid,qo_id,'long_name',
     &            15, 'Ocean vorticity')
      if ( status.ne.NF_NOERR ) call handle_err (status, subnam)

!!  Define a 2-dimensional variable called 'uptpoc'
!!  to store ocean <u'T'>
      dims(1) = xopdim
      dims(2) = yotdim
      status = nf_def_var(tavncid, 'uptpoc', NF_DOUBLE,
     &                    2, dims, utoc_id)
      if ( status.ne.NF_NOERR ) call handle_err (status, subnam)
      status = nf_put_att_text(tavncid,utoc_id,'units',5,'K.m/s')
      if ( status.ne.NF_NOERR ) call handle_err (status, subnam)
      status = nf_put_att_text(tavncid,utoc_id,'long_name',
     &            12, 'Ocean <upTp>')
      if ( status.ne.NF_NOERR ) call handle_err (status, subnam)

!!  Define a 2-dimensional variable called 'vptpoc'
!!  to store ocean <v'T'>
      dims(1) = xotdim
      dims(2) = yopdim
      status = nf_def_var(tavncid, 'vptpoc', NF_DOUBLE,
     &                    2, dims, vtoc_id)
      if ( status.ne.NF_NOERR ) call handle_err (status, subnam)
      status = nf_put_att_text(tavncid,vtoc_id,'units',5,'K.m/s')
      if ( status.ne.NF_NOERR ) call handle_err (status, subnam)
      status = nf_put_att_text(tavncid,vtoc_id,'long_name',
     &            12, 'Ocean <vpTp>')
      if ( status.ne.NF_NOERR ) call handle_err (status, subnam)
#  endif /* not atmos_only */

#  ifndef ocean_only
!!  Define a 2-dimensional variable called 'ast'
!!  to store atmos surface temperatures
      dims(1) = xatdim
      dims(2) = yatdim
      status = nf_def_var(tavncid, 'ast', NF_DOUBLE, 2, dims, ast_id)
      if ( status.ne.NF_NOERR ) call handle_err (status, subnam)
      status = nf_put_att_text(tavncid,ast_id,'units',1,'K')
      if ( status.ne.NF_NOERR ) call handle_err (status, subnam)
      status = nf_put_att_text(tavncid,ast_id,'long_name',
     &            25, 'Atmos surface temperature')
      if ( status.ne.NF_NOERR ) call handle_err (status, subnam)

!!  Define a 2-dimensional variable called 'wekta'
!!  to store atmos Ekman velocity (on T-grid)
      dims(1) = xatdim
      dims(2) = yatdim
      status = nf_def_var(tavncid, 'wekta', NF_DOUBLE,
     &                    2, dims, wekta_id)
      if ( status.ne.NF_NOERR ) call handle_err (status, subnam)
      status = nf_put_att_text(tavncid,wekta_id,'units',3,'m/s')
      if ( status.ne.NF_NOERR ) call handle_err (status, subnam)
      status = nf_put_att_text(tavncid,wekta_id,'long_name',
     &            29, 'Atmos Ekman velocity (T-grid)')
      if ( status.ne.NF_NOERR ) call handle_err (status, subnam)

!!  Define a 2-dimensional variable called 'fnetat'
!!  to store atmos forcing
      dims(1) = xatdim
      dims(2) = yatdim
      status = nf_def_var(tavncid, 'fnetat', NF_DOUBLE,
     &                    2, dims, fmat_id)
      if ( status.ne.NF_NOERR ) call handle_err (status, subnam)
      status = nf_put_att_text(tavncid,fmat_id,'units',5,'W/m^2')
      if ( status.ne.NF_NOERR ) call handle_err (status, subnam)
      status = nf_put_att_text(tavncid,fmat_id,'long_name',
     &            30, 'Atmos surface forcing (T-grid)')
      if ( status.ne.NF_NOERR ) call handle_err (status, subnam)

!!  Define a 2-dimensional variable called 'tauxa'
!!  to store atmos stress
      dims(1) = xapdim
      dims(2) = yapdim
      status = nf_def_var(tavncid, 'tauxa', NF_DOUBLE, 2, dims, txa_id)
      if ( status.ne.NF_NOERR ) call handle_err (status, subnam)
      status = nf_put_att_text(tavncid,txa_id,'units',7,'m^2/s^2')
      if ( status.ne.NF_NOERR ) call handle_err (status, subnam)
      status = nf_put_att_text(tavncid,txa_id,'long_name',
     &            31, 'Atmos surface x stress (p-grid)')
      if ( status.ne.NF_NOERR ) call handle_err (status, subnam)

!!  Define a 2-dimensional variable called 'tauya'
!!  to store atmos stress
      dims(1) = xapdim
      dims(2) = yapdim
      status = nf_def_var(tavncid, 'tauya', NF_DOUBLE, 2, dims, tya_id)
      if ( status.ne.NF_NOERR ) call handle_err (status, subnam)
      status = nf_put_att_text(tavncid,tya_id,'units',7,'m^2/s^2')
      if ( status.ne.NF_NOERR ) call handle_err (status, subnam)
      status = nf_put_att_text(tavncid,tya_id,'long_name',
     &            31, 'Atmos surface y stress (p-grid)')
      if ( status.ne.NF_NOERR ) call handle_err (status, subnam)

!!  Define a 3-dimensional variable called 'pa'
!!  to store atmos dynamic pressure
      pdims(1) = xapdim
      pdims(2) = yapdim
      pdims(3) = ladim
      status = nf_def_var(tavncid, 'pa', NF_DOUBLE, 3, pdims, pa_id)
      if ( status.ne.NF_NOERR ) call handle_err (status, subnam)
      status = nf_put_att_text(tavncid,pa_id,'units',7,'m^2/s^2')
      if ( status.ne.NF_NOERR ) call handle_err (status, subnam)
      status = nf_put_att_text(tavncid,pa_id,'long_name',
     &            22, 'Atmos dynamic pressure')
      if ( status.ne.NF_NOERR ) call handle_err (status, subnam)

!!  Define a 3-dimensional variable called 'qa'
!!  to store atmos vorticity
      pdims(1) = xapdim
      pdims(2) = yapdim
      pdims(3) = ladim
      status = nf_def_var(tavncid, 'qa', NF_DOUBLE, 3, pdims, qa_id)
      if ( status.ne.NF_NOERR ) call handle_err (status, subnam)
      status = nf_put_att_text(tavncid,qa_id,'units',4,'s^-1')
      if ( status.ne.NF_NOERR ) call handle_err (status, subnam)
      status = nf_put_att_text(tavncid,qa_id,'long_name',
     &            15, 'Atmos vorticity')
      if ( status.ne.NF_NOERR ) call handle_err (status, subnam)

!!  Define a 2-dimensional variable called 'uptpat'
!!  to store atmos <u'T'>
      dims(1) = xapdim
      dims(2) = yatdim
      status = nf_def_var(tavncid, 'uptpat', NF_DOUBLE,
     &                    2, dims, utat_id)
      if ( status.ne.NF_NOERR ) call handle_err (status, subnam)
      status = nf_put_att_text(tavncid,utat_id,'units',5,'K.m/s')
      if ( status.ne.NF_NOERR ) call handle_err (status, subnam)
      status = nf_put_att_text(tavncid,utat_id,'long_name',
     &            12, 'Atmos <upTp>')
      if ( status.ne.NF_NOERR ) call handle_err (status, subnam)

!!  Define a 2-dimensional variable called 'vptpat'
!!  to store atmos <v'T'>
      dims(1) = xatdim
      dims(2) = yapdim
      status = nf_def_var(tavncid, 'vptpat', NF_DOUBLE,
     &                    2, dims, vtat_id)
      if ( status.ne.NF_NOERR ) call handle_err (status, subnam)
      status = nf_put_att_text(tavncid,vtat_id,'units',5,'K.m/s')
      if ( status.ne.NF_NOERR ) call handle_err (status, subnam)
      status = nf_put_att_text(tavncid,vtat_id,'long_name',
     &            12, 'Atmos <vpTp>')
      if ( status.ne.NF_NOERR ) call handle_err (status, subnam)
#  endif /* not ocean_only */

!! Leave definition mode: entering data mode.
      status = nf_enddef(tavncid)
      if ( status.ne.NF_NOERR ) call handle_err (status, subnam)

#  ifndef atmos_only
!!    Oceanic axes
!!    ------------
!!    Convert all x- and y-axes into km
!!    Calculate x gridpoints and store in 'x' arrays
!!    p-grid points
      do i=1,nxpo
        xxo(i) = 1.0d-3*( xpo(i) - xpo(1) )
      enddo
      status = nf_put_vara_double(tavncid, xop_id, 1, nxpo, xxo)
      if ( status.ne.NF_NOERR ) call handle_err (status, subnam)
!!    T-grid points
      do i=1,nxto
        xxo(i) = 1.0d-3*( xto(i) - xpo(1) )
      enddo
      status = nf_put_vara_double(tavncid, xot_id, 1, nxto, xxo)
      if ( status.ne.NF_NOERR ) call handle_err (status, subnam)

!!    Calculate y gridpoints and store in 'y' arrays
!!    p-grid points
      do j=1,nypo
        yyo(j) = 1.0d-3*( ypo(j) - ypo(1) )
      enddo
      status = nf_put_vara_double(tavncid, yop_id, 1, nypo, yyo)
      if ( status.ne.NF_NOERR ) call handle_err (status, subnam)
!!    T-grid points
      do j=1,nyto
        yyo(j) = 1.0d-3*( yto(j) - ypo(1) )
      enddo
      status = nf_put_vara_double(tavncid, yot_id, 1, nyto, yyo)
      if ( status.ne.NF_NOERR ) call handle_err (status, subnam)

!!    Convert mid-layer depth into km and store in 'z'
      tmpo(1) = 0.5d-3*hoc(1)
      do k=2,nlo
        tmpo(k) = tmpo(k-1) + 0.5d-3*( hoc(k-1) + hoc(k) )
      enddo
      status = nf_put_var_double(tavncid, lo_id, tmpo)
      if ( status.ne.NF_NOERR ) call handle_err (status, subnam)
#  endif /* not atmos_only */

#  ifndef ocean_only
!!    Atmospheric axes
!!    ----------------
!!    Convert all x- and y-axes into km
!!    Calculate x gridpoints and store in 'x' arrays
!!    p-grid points
      do i=1,nxpa
        xxa(i) = 1.0d-3*xpa(i)
      enddo
      status = nf_put_vara_double(tavncid, xap_id, 1, nxpa, xxa)
      if ( status.ne.NF_NOERR ) call handle_err (status, subnam)
!!    T-grid points
      do i=1,nxta
        xxa(i) = 1.0d-3*xta(i)
      enddo
      status = nf_put_vara_double(tavncid, xat_id, 1, nxta, xxa)
      if ( status.ne.NF_NOERR ) call handle_err (status, subnam)

!!    Calculate y gridpoints and store in 'y' arrays
!!    p-grid points
      do j=1,nypa
        yya(j) = 1.0d-3*ypa(j)
      enddo
      status = nf_put_vara_double(tavncid, yap_id, 1, nypa, yya)
      if ( status.ne.NF_NOERR ) call handle_err (status, subnam)
!!    T-grid points
      do j=1,nyta
        yya(j) = 1.0d-3*yta(j)
      enddo
      status = nf_put_vara_double(tavncid, yat_id, 1, nyta, yya)
      if ( status.ne.NF_NOERR ) call handle_err (status, subnam)

!!    Convert mid-layer height into km and store in 'z'
      tmpa(1) = 0.5d-3*hat(1)
      do k=2,nla
        tmpa(k) = tmpa(k-1) + 0.5d-3*( hat(k-1) + hat(k) )
      enddo
      status = nf_put_var_double(tavncid, la_id, tmpa)
      if ( status.ne.NF_NOERR ) call handle_err (status, subnam)
#  endif /* not ocean_only */

#  ifndef atmos_only
!!    Output oceanic data
!!    -------------------
      start(1) = 1
      start(2) = 1
      count(1) = nxto
      count(2) = nyto
      status=nf_put_vara_double(tavncid, sst_id, start, count, sstav)
      if ( status.ne.NF_NOERR ) call handle_err (status, subnam)
      status=nf_put_vara_double(tavncid, wekto_id, start, count, wtocav)
      if ( status.ne.NF_NOERR ) call handle_err (status, subnam)
      status=nf_put_vara_double(tavncid, fmoc_id, start, count, fmocav)
      if ( status.ne.NF_NOERR ) call handle_err (status, subnam)
      count(1) = nxpo
      count(2) = nypo
      status=nf_put_vara_double(tavncid, txo_id, start, count, txocav)
      if ( status.ne.NF_NOERR ) call handle_err (status, subnam)
      status=nf_put_vara_double(tavncid, tyo_id, start, count, tyocav)
      if ( status.ne.NF_NOERR ) call handle_err (status, subnam)
      status=nf_put_vara_double(tavncid, wekpo_id, start, count, wpocav)
      if ( status.ne.NF_NOERR ) call handle_err (status, subnam)
      startp(1) = 1
      startp(2) = 1
      startp(3) = 1
      countp(1) = nxpo
      countp(2) = nypo
      countp(3) = nlo
      status=nf_put_vara_double(tavncid, po_id, startp, countp, pocav)
      if ( status.ne.NF_NOERR ) call handle_err (status, subnam)
      status=nf_put_vara_double(tavncid, qo_id, startp, countp, qocav)
      if ( status.ne.NF_NOERR ) call handle_err (status, subnam)
      count(1) = nxpo
      count(2) = nyto
      status=nf_put_vara_double(tavncid, utoc_id, start, count, uptpoc)
      if ( status.ne.NF_NOERR ) call handle_err (status, subnam)
      count(1) = nxto
      count(2) = nypo
      status=nf_put_vara_double(tavncid, vtoc_id, start, count, vptpoc)
      if ( status.ne.NF_NOERR ) call handle_err (status, subnam)
#  endif /* not atmos_only */

#  ifndef ocean_only
!!    Output atmospheric data
!!    -----------------------
      start(1) = 1
      start(2) = 1
      count(1) = nxta
      count(2) = nyta
      status=nf_put_vara_double(tavncid, ast_id, start, count, astav)
      if ( status.ne.NF_NOERR ) call handle_err (status, subnam)
      status=nf_put_vara_double(tavncid, wekta_id, start, count, wtatav)
      if ( status.ne.NF_NOERR ) call handle_err (status, subnam)
      status=nf_put_vara_double(tavncid, fmat_id, start, count, fmatav)
      if ( status.ne.NF_NOERR ) call handle_err (status, subnam)
      count(1) = nxpa
      count(2) = nypa
      status=nf_put_vara_double(tavncid, txa_id, start, count, txatav)
      if ( status.ne.NF_NOERR ) call handle_err (status, subnam)
      status=nf_put_vara_double(tavncid, tya_id, start, count, tyatav)
      if ( status.ne.NF_NOERR ) call handle_err (status, subnam)
      startp(1) = 1
      startp(2) = 1
      startp(3) = 1
      countp(1) = nxpa
      countp(2) = nypa
      countp(3) = nla
      status=nf_put_vara_double(tavncid, pa_id, startp, countp, patav)
      if ( status.ne.NF_NOERR ) call handle_err (status, subnam)
      status=nf_put_vara_double(tavncid, qa_id, startp, countp, qatav)
      if ( status.ne.NF_NOERR ) call handle_err (status, subnam)
      count(1) = nxpa
      count(2) = nyta
      status=nf_put_vara_double(tavncid, utat_id, start, count, uptpat)
      if ( status.ne.NF_NOERR ) call handle_err (status, subnam)
      count(1) = nxta
      count(2) = nypa
      status=nf_put_vara_double(tavncid, vtat_id, start, count, vptpat)
      if ( status.ne.NF_NOERR ) call handle_err (status, subnam)
#  endif /* not ocean_only */

#endif /* use_netcdf */

      END SUBROUTINE tavout
c
c***********************************************************************
*
      END MODULE timavge
*
c***********************************************************************
