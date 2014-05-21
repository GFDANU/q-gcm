c***********************************************************************
c     Q-GCM Version 1.5.0 : last modified 01/07/2013
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
      MODULE covaria

*     Contains subprograms for computing and outputting
*     (in netCDF format) averages and covariances of
*     various spatially subsampled quantities during a run.
*     Ensure that the internal arrays are only defined if
*     they are actually needed, as they tend to be large.

*     Modules
      USE parameters, ONLY : nmcvat, nvcvat, nmcvoc, nvcvoc

      IMPLICIT NONE

      PRIVATE

      integer, PUBLIC, SAVE :: covncid

#ifdef get_covar

      PUBLIC  :: covini, covatm, covocn, covout
      PRIVATE :: tsampl, psampl, dssp

#  ifndef ocean_only
*     Storage for matrices of covariance, vectors of averages and
*     sums of weights at subsampled points, for atmospheric data

      integer, PRIVATE, SAVE :: nupa,nuta
      double precision, PRIVATE, SAVE :: covpa(nmcvat), avgpa(nvcvat),
     &                                   covta(nmcvat), avgta(nvcvat),
     &                                   swtpa, swtta
#  endif /* not ocean_only */

#  ifndef atmos_only
*     Storage for matrices of covariance, vectors of averages
*     and sums of weights at subsampled points, for oceanic data

      integer, PRIVATE, SAVE :: nupo,nuto
      double precision, PRIVATE, SAVE :: covpo(nmcvoc), avgpo(nvcvoc),
     &                                   covto(nmcvoc), avgto(nvcvoc),
     &                                   swtpo, swtto
#  endif /* not atmos_only */

      CONTAINS
#endif /* get_covar */

c***********************************************************************
c
#ifdef get_covar
*
      SUBROUTINE covini

*     Initialises the storage used for
*     computing time covariance matrices

*     Modules
      USE parameters, ONLY : nmcvat, nvcvat, nmcvoc, nvcvoc

      IMPLICIT NONE

*     Subroutine arguments
*
*     Local variables
      integer k

#  ifndef ocean_only
      nupa = 0
      nuta = 0
      swtpa = 0.0d0
      swtta = 0.0d0

*     Initialise atmosphere fields
*     ----------------------------
!$OMP PARALLEL DEFAULT (NONE)
!$OMP&         PRIVATE (k)
!$OMP&         SHARED  (covpa,covta,avgpa,avgta)

!$OMP DO SCHEDULE (STATIC)
      do k=1,nmcvat
        covpa(k) = 0.0d0
        covta(k) = 0.0d0
      enddo
!$OMP END DO NOWAIT
!$OMP DO SCHEDULE (STATIC)
      do k=1,nvcvat
        avgpa(k) = 0.0d0
        avgta(k) = 0.0d0
      enddo
!$OMP END DO NOWAIT

!$OMP END PARALLEL
#  endif /* not ocean_only */

#  ifndef atmos_only
      nupo = 0
      nuto = 0
      swtpo = 0.0d0
      swtto = 0.0d0

*     Initialise ocean fields
*     -----------------------
!$OMP PARALLEL DEFAULT (NONE)
!$OMP&         PRIVATE (k)
!$OMP&         SHARED  (covpo,covto,avgpo,avgto)

!$OMP DO SCHEDULE (STATIC)
      do k=1,nmcvoc
        covpo(k) = 0.0d0
        covto(k) = 0.0d0
      enddo
!$OMP END DO NOWAIT
!$OMP DO SCHEDULE (STATIC)
      do k=1,nvcvoc
        avgpo(k) = 0.0d0
        avgto(k) = 0.0d0
      enddo
!$OMP END DO NOWAIT

!$OMP END PARALLEL
#  endif /* not atmos_only */

      END SUBROUTINE covini

#endif /* get_covar */
c
c***********************************************************************
c
#ifdef get_covar
*
      SUBROUTINE covatm

*     Accumulates atmospheric time covariance matrices during a run

*     Modules
      USE parameters, ONLY : nxpa, nypa, nxta, nyta, nscvat, nvcvat
      USE atstate, ONLY : pa
      USE intrfac, ONLY : ast

      IMPLICIT NONE

*     Subroutine arguments
*
#  ifndef ocean_only

*     Local variables
      integer ifault
      double precision u(nvcvat)

*     Spatially subsample array pa into vector u
      call psampl (pa(1,1,1), nxpa, nypa, nscvat, u, nvcvat)
*     Update covariance matrix
      call dssp (u, avgpa, covpa, 1.0d0, swtpa, nvcvat, nupa, ifault)
      if ( ifault.ne.0 ) then
        print *,' dssp problem in covatm on pa; ifault = ',ifault
      endif

*     Spatially subsample array ast into vector u
      call tsampl (ast, nxta, nyta, nscvat, u, nvcvat)
*     Update covariance matrix
      call dssp (u, avgta, covta, 1.0d0, swtta, nvcvat, nuta, ifault)
      if ( ifault.ne.0 ) then
        print *,' dssp problem in covatm on ast; ifault = ',ifault
      endif

#  endif /* not ocean_only */

      END SUBROUTINE covatm

#endif /* get_covar */
c
c***********************************************************************
c
#ifdef get_covar
*
      SUBROUTINE covocn

*     Accumulates oceanic time covariance matrices during a run

*     Modules
      USE parameters, ONLY : nxpo, nypo, nxto, nyto, nscvoc, nvcvoc
      USE ocstate, ONLY : po
      USE intrfac, ONLY : sst

      IMPLICIT NONE

*     Subroutine arguments
*
#  ifndef atmos_only

*     Local variables
      integer ifault
      double precision u(nvcvoc)

*     Spatially subsample array po into vector u
      call psampl (po(1,1,1), nxpo, nypo, nscvoc, u, nvcvoc)
*     Update covariance matrix
      call dssp (u, avgpo, covpo, 1.0d0, swtpo, nvcvoc, nupo, ifault)
      if ( ifault.ne.0 ) then
        print *,' dssp problem in covocn on po; ifault = ',ifault
      endif

*     Spatially subsample array sst into vector u
      call tsampl (sst, nxto, nyto, nscvoc, u, nvcvoc)
*     Update covariance matrix
      call dssp (u, avgto, covto, 1.0d0, swtto, nvcvoc, nuto, ifault)
      if ( ifault.ne.0 ) then
        print *,' dssp problem in covocn on sst; ifault = ',ifault
      endif

#  endif /* not atmos_only */

      END SUBROUTINE covocn

#endif /* get_covar */
c
c***********************************************************************
c
#ifdef get_covar
*
      SUBROUTINE covout

*     Writes oceanic and atmospheric time covariance
*     matrices at the end of a run, using netCDF

*     Modules
      USE parameters, ONLY : nmcvat, nvcvat, nmcvoc, nvcvoc
      USE nc_subs, ONLY : handle_err

      IMPLICIT NONE

#  ifdef use_netcdf
      INCLUDE 'netcdf.inc'
#  endif
*
*     Subroutine arguments
*
*     Local parameters
      character (len=*), parameter :: subnam = 'covout'

*     Local variables
#  ifdef use_netcdf
      integer status,idim
#    ifndef ocean_only
      integer ncadim,nvadim,covpa_id,covta_id,
     &        avgpa_id,avgta_id,swtpa_id,swtta_id
#    endif
#    ifndef atmos_only
      integer ncodim,nvodim,covpo_id,covto_id,
     &        avgpo_id,avgto_id,swtpo_id,swtto_id
#    endif

!!    Initialise netCDF file
!!    ----------------------
      status = nf_def_dim(covncid, 's', 1, idim)
      if (status.ne.NF_NOERR) call handle_err (status, subnam)

#    ifndef ocean_only
      status = nf_def_dim(covncid, 'nmcvat', nmcvat, ncadim)
      if (status.ne.NF_NOERR) call handle_err (status, subnam)
      status = nf_def_dim(covncid, 'nvcvat', nvcvat, nvadim)
      if (status.ne.NF_NOERR) call handle_err (status, subnam)
      status=nf_def_var(covncid, 'covpa', NF_FLOAT, 1, ncadim, covpa_id)
      if (status.ne.NF_NOERR) call handle_err (status, subnam)
      status=nf_def_var(covncid, 'covta', NF_FLOAT, 1, ncadim, covta_id)
      if (status.ne.NF_NOERR) call handle_err (status, subnam)
      status=nf_def_var(covncid, 'avgpa', NF_FLOAT, 1, nvadim, avgpa_id)
      if (status.ne.NF_NOERR) call handle_err (status, subnam)
      status=nf_def_var(covncid, 'avgta', NF_FLOAT, 1, nvadim, avgta_id)
      if (status.ne.NF_NOERR) call handle_err (status, subnam)
      status=nf_def_var(covncid, 'swtpa', NF_FLOAT, 1, idim, swtpa_id)
      if (status.ne.NF_NOERR) call handle_err (status, subnam)
      status=nf_def_var(covncid, 'swtta', NF_FLOAT, 1, idim, swtta_id)
      if (status.ne.NF_NOERR) call handle_err (status, subnam)
#    endif /* not ocean_only */
#    ifndef atmos_only
      status = nf_def_dim(covncid, 'nmcvoc', nmcvoc, ncodim)
      if (status.ne.NF_NOERR) call handle_err (status, subnam)
      status = nf_def_dim(covncid, 'nvcvoc', nvcvoc, nvodim)
      if (status.ne.NF_NOERR) call handle_err (status, subnam)
      status=nf_def_var(covncid, 'covpo', NF_FLOAT, 1, ncodim, covpo_id)
      if (status.ne.NF_NOERR) call handle_err (status, subnam)
      status=nf_def_var(covncid, 'covto', NF_FLOAT, 1, ncodim, covto_id)
      if (status.ne.NF_NOERR) call handle_err (status, subnam)
      status=nf_def_var(covncid, 'avgpo', NF_FLOAT, 1, nvodim, avgpo_id)
      if (status.ne.NF_NOERR) call handle_err (status, subnam)
      status=nf_def_var(covncid, 'avgto', NF_FLOAT, 1, nvodim, avgto_id)
      if (status.ne.NF_NOERR) call handle_err (status, subnam)
      status=nf_def_var(covncid, 'swtpo', NF_FLOAT, 1, idim, swtpo_id)
      if (status.ne.NF_NOERR) call handle_err (status, subnam)
      status=nf_def_var(covncid, 'swtto', NF_FLOAT, 1, idim, swtto_id)
      if (status.ne.NF_NOERR) call handle_err (status, subnam)
#    endif /* not atmos_only */

!!    Leave definition mode: entering data mode.
      status = nf_enddef(covncid)
      if (status.ne.NF_NOERR) call handle_err (status, subnam)

!!    Write data to arrays
!!    --------------------
#    ifndef ocean_only
      status = nf_put_vara_double(covncid, covpa_id, 1, nmcvat, covpa)
      if (status.ne.NF_NOERR) call handle_err (status, subnam)
      status = nf_put_vara_double(covncid, covta_id, 1, nmcvat, covta)
      if (status.ne.NF_NOERR) call handle_err (status, subnam)
      status = nf_put_vara_double(covncid, avgpa_id, 1, nvcvat, avgpa)
      if (status.ne.NF_NOERR) call handle_err (status, subnam)
      status = nf_put_vara_double(covncid, avgta_id, 1, nvcvat, avgta)
      if (status.ne.NF_NOERR) call handle_err (status, subnam)
      status = nf_put_vara_double(covncid, swtpa_id, 1, 1, swtpa)
      if (status.ne.NF_NOERR) call handle_err (status, subnam)
      status = nf_put_vara_double(covncid, swtta_id, 1, 1, swtta)
      if (status.ne.NF_NOERR) call handle_err (status, subnam)
#    endif /* not ocean_only */
#    ifndef atmos_only
      status = nf_put_vara_double(covncid, covpo_id, 1, nmcvoc, covpo)
      if (status.ne.NF_NOERR) call handle_err (status, subnam)
      status = nf_put_vara_double(covncid, covto_id, 1, nmcvoc, covto)
      if (status.ne.NF_NOERR) call handle_err (status, subnam)
      status = nf_put_vara_double(covncid, avgpo_id, 1, nvcvoc, avgpo)
      if (status.ne.NF_NOERR) call handle_err (status, subnam)
      status = nf_put_vara_double(covncid, avgto_id, 1, nvcvoc, avgto)
      if (status.ne.NF_NOERR) call handle_err (status, subnam)
      status = nf_put_vara_double(covncid, swtpo_id, 1, 1, swtpo)
      if (status.ne.NF_NOERR) call handle_err (status, subnam)
      status = nf_put_vara_double(covncid, swtto_id, 1, 1, swtto)
      if (status.ne.NF_NOERR) call handle_err (status, subnam)
#    endif /* not atmos_only */
#  endif /* use_netcdf */

      END SUBROUTINE covout

#endif /* get_covar */
c
c***********************************************************************
c
#ifdef get_covar
*
      SUBROUTINE tsampl (datat, nx, ny, nsi, dsamt, nsvec)
*
*     Given an array datat(nx,ny) of data tabulated at T points, and
*     a subsampling interval nsi, computes the values of dsamt(nsvec),
*     a vector containing the required subsample of the data
*     This will work for both channel and box configurations

      IMPLICIT NONE

*     Subroutine arguments
      integer nx,ny,nsi,nsvec
      double precision datat(nx,ny),dsamt(nsvec)
*
*     Local variables
      integer is,js,ivs,id,id1,id2,jd,jd1,jd2
      double precision sumd,sumi

      if ( nsi.gt.1 ) then
*       Subsample the data by computing the arithmetic mean
!$OMP   PARALLEL DO DEFAULT (NONE)
!$OMP&           PRIVATE (is,js,id1,id2,id,jd1,jd2,jd,sumd,sumi,ivs)
!$OMP&           SHARED  (nx,ny,nsi,datat,dsamt)
!$OMP&           SCHEDULE (STATIC)
        do js=1,(ny/nsi)
          jd1 = 1 + (js-1)*nsi
          jd2 = js*nsi
          do is=1,(nx/nsi)
            id1 = 1 + (is-1)*nsi
            id2 = is*nsi
            sumd = 0.0d0
            do jd=jd1,jd2
              sumi = 0.0d0
              do id=id1,id2
                sumi = sumi + datat(id,jd)
              enddo
              sumd = sumd + sumi
            enddo
            ivs = (js-1)*(nx/nsi) + is
            dsamt(ivs) = sumd
          enddo
        enddo
!$OMP   END PARALLEL DO
       else if ( nsi.eq.1 ) then
*       Just copy the data into the subsample array
!$OMP   PARALLEL DO DEFAULT (NONE)
!$OMP&           PRIVATE (id,jd,ivs)
!$OMP&           SHARED  (nx,ny,datat,dsamt)
!$OMP&           SCHEDULE (STATIC)
        do jd=1,ny
          do id=1,nx
            ivs = (jd-1)*nx + id
            dsamt(ivs) = datat(id,jd)
          enddo
        enddo
!$OMP   END PARALLEL DO
       else
        print *,' WARNING: tsampl called with invalid nsi = ',nsi
      endif
      ivs = (nx/nsi)*(ny/nsi)
      if ( ivs.ne.nsvec ) then
        print *,' WARNING: inconsistent subsample vector in tsampl'
        print *,' nsvec, final ivs = ',nsvec,ivs
      endif

      END SUBROUTINE tsampl

#endif /* get_covar */
c
c***********************************************************************
c
#ifdef get_covar
*
      SUBROUTINE psampl (datap, nx, ny, nsi, dsamp, nsvec)
*
*     Given an array datap(nx,ny) of data tabulated at p points, and
*     a subsampling interval nsi, computes the values of dsamp(nvec),
*     a vector containing the required subsample of data

      IMPLICIT NONE

*     Subroutine arguments
      integer nx,ny,nsi,nsvec
      double precision datap(nx,ny),dsamp(nsvec)
*
*     Local variables
      integer is,js,ivs,id,id1,id2,jd,jd1,jd2
      double precision sumd,sumi,sums,sumn

!$OMP PARALLEL DO DEFAULT (NONE)
!$OMP&         PRIVATE (is,js,id1,id2,id,jd1,jd2,jd,
!$OMP&                  sumd,sumi,sums,sumn,ivs)
!$OMP&         SHARED  (nx,ny,nsi,datap,dsamp)
!$OMP&         SCHEDULE (STATIC)
      do js=1,(ny/nsi)
        jd1 = 1 + (js-1)*nsi
        jd2 = 1 + js*nsi
        do is=1,(nx/nsi)
          id1 = 1 + (is-1)*nsi
          id2 = 1 + is*nsi
*         Sum over inner rows of averaging area
          sumd = 0.0d0
          do jd=jd1+1,jd2-1
            sumi = 0.5d0*datap(id1,jd)
            do id=id1+1,id2-1
              sumi = sumi + datap(id,jd)
            enddo
            sumi = sumi + 0.5d0*datap(id2,jd)
            sumd = sumd + sumi
          enddo
*         Contributions from S & N boundaries of averaging area
          sums = 0.5d0*datap(id1,jd1)
          sumn = 0.5d0*datap(id1,jd2)
          do id=id1+1,id2-1
            sums = sums + datap(id,jd1)
            sumn = sumn + datap(id,jd2)
          enddo
          sums = sums + 0.5d0*datap(id2,jd1)
          sumn = sumn + 0.5d0*datap(id2,jd2)
          ivs = (js-1)*(nx/nsi) + is
          dsamp(ivs) = sumd + 0.5d0*( sums + sumn )
        enddo
      enddo
!$OMP END PARALLEL DO
      ivs = (nx/nsi)*(ny/nsi)
      if ( ivs.ne.nsvec ) then
        print *,' WARNING: inconsistent subsample vector in psampl'
        print *,' nsvec, final ivs = ',nsvec,ivs
      endif

      END SUBROUTINE psampl

#endif /* get_covar */
c
c***********************************************************************
c
#ifdef get_covar
c
      SUBROUTINE dssp (x, xmean, xssp, wt, sumwt, nvar, nunit, ifault)
c
c       Algorithm AS 41 j.r.statist.soc.c. (1971) vol. 20 no.2
c
c       This subroutine updates the mean vector xmean (length nvar)
c       and the matrix of corrected sums of squares and products xssp
c       (length nvar(nvar+1)/2, stored by lower triangle), when a
c       data vector x (length nvar) with weight wt is either included
c       (wt.gt.0) or excluded (wt.lt.0).  sumwt is the current sum of
c       weights on entry and the updated sum on exit and nunit is
c       the current and updated sample size.  ifault=0 indicates normal
c       exit,  ifault=1 indicates zero or negative value of sumwt,
c       ifault=2 indicates zero or negative nunit, ifault=3 indicates
c       nvar.lt.1.  Note that x, xmean, xssp, wt and sumwt are double
c       precision and must be declared as such in the calling program.

      IMPLICIT NONE

*     Subroutine arguments
      double precision x(*), xmean(*), xssp(*)
      double precision wt, sumwt
      integer nvar,nunit,ifault

*     Local parameters
      double precision co
      parameter ( co=0.0d0 )

*     Local variables
      integer i,j,k
      double precision b, c
c
c          Check variates, weights and sample size
c
      ifault = 0
      if ( nvar.lt.1 ) then
        ifault = 3
        return
      endif
      if ( wt.lt.0.0d0 ) then
        nunit = nunit-1
       else if ( wt.eq.0.0d0 ) then
        return
       else if ( wt.gt.0.0d0 ) then
        nunit = nunit+1
      endif
      sumwt = sumwt+wt
      if ( sumwt.le.co ) then
        ifault = 1
        return
      endif
      b = wt/sumwt
      if ( nunit.lt.1 ) then
        ifault = 2
       else if ( nunit.eq.1 ) then
c
c       Initialise means and ssp for sample size  =  1
c
        do i = 1,nvar
          if ( wt.lt.co ) then
            xmean(i) = xmean(i) + b*(x(i)-xmean(i))
           else
            xmean(i) = x(i)
          endif
          x(i) = co
          k = (i*(i-1))/2
          do j = 1,i
            k = k + 1
            xssp(k) = co
          enddo
        enddo
       else if ( nunit.gt.1 ) then
c
c       Update means and ssp for sample size greater than 1
c
        c = wt - b*wt
!$OMP   PARALLEL DEFAULT (NONE)
!$OMP&           PRIVATE (i,j,k)
!$OMP&           SHARED  (nvar,x,xmean,b,xssp,c)

!$OMP   DO SCHEDULE (STATIC)
        do i = 1,nvar
          x(i) = x(i) - xmean(i)
          xmean(i) = xmean(i) + b*x(i)
        enddo
!$OMP   END DO
!$OMP   DO SCHEDULE (RUNTIME)
        do i = 1,nvar
          k = (i*(i-1))/2
          do j = 1,i
            k = k + 1
            xssp(k) = xssp(k) + c*x(i)*x(j)
          enddo
        enddo
!$OMP   END DO NOWAIT

!$OMP   END PARALLEL

      endif

      END SUBROUTINE dssp

#endif /* get_covar */
c
c***********************************************************************
*
      END MODULE covaria
*
c***********************************************************************
