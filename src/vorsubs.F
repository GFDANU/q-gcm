c***********************************************************************
c     Q-GCM Version 1.5.0 : last modified 18/02/2013
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
      MODULE vorsubs

*     Contains subroutines qcomp, merqcy, ocqbdy and atqzbd which
*     compute potential vorticity q from dynamic pressure p,
*     respectively at internal points (qcomp), along meridional
*     boundaries in cyclic domains (merqcy), on solid oceanic
*     boundaries (ocqbdy), and on zonal atmospheric boundaries (atqzbd).

*     Modules

      IMPLICIT NONE

      PRIVATE

      PUBLIC :: qcomp, merqcy
#ifndef atmos_only
      PUBLIC :: ocqbdy
#endif
#ifndef ocean_only
      PUBLIC :: atqzbd
#endif

      CONTAINS

c***********************************************************************
c
      SUBROUTINE qcomp (q, p, aaa, yprel, dxm2,
     &                  nxp, nyp, nl, ddyn, kbot)

*     Compute pot. vort. q from p at internal points (equation 7.15)
*     Works for both atmosphere and ocean

*     Input arguments:
*     p     : 3-D array of dynamic pressure at p points
*     aaa   : A(nl,nl) matrix linking pressures and eta
*     yprel : y-coordinate of p points rel. to central latitude
*     dxm2  : dx^-2
*     nxp   : x-dimension of p array
*     nyp   : y-dimension of p array
*     nl    : no. of QG layers
*     ddyn  : dynamic (rescaled) topography
*     kbot  : index of layer containing topography

*     Output arguments:
*     q    : 3-D array of vorticity at p points

*     Modules
      USE parameters

      IMPLICIT NONE

*     Subroutine arguments
      integer nxp,nyp,nl,kbot
      double precision q(nxp,nyp,nl),p(nxp,nyp,nl),aaa(nl,nl),
     &                 yprel(nyp),dxm2,ddyn(nxp,nyp)
*
*     Local variables
      integer i,j,k
      double precision dx2fac,betay

      dx2fac = dxm2/fnot

!$OMP PARALLEL DEFAULT (NONE)
!$OMP&         PRIVATE (i,j,k,betay)
!$OMP&         SHARED  (nxp,nyp,nl,yprel,q,dx2fac,p,aaa,ddyn,kbot)

*     k = 1 layer (top(ocean) or bottom(atmos.))
!$OMP DO SCHEDULE (STATIC)
      do j=2,nyp-1
        betay = beta*yprel(j)
        do i=2,nxp-1
          q(i,j, 1) = dx2fac*(  p(i,j-1,1) + p(i-1,j,1) + p(i+1,j,1)
     &                        + p(i,j+1,1) - 4.0d0*p(i,j,1) ) + betay
     &               - fnot*( aaa(1,1)*p(i,j,1) + aaa(1,2)*p(i,j,2) )
        enddo
      enddo
!$OMP END DO NOWAIT
*     Intermediate layers
      do k=2,nl-1
!$OMP   DO SCHEDULE (STATIC)
        do j=2,nyp-1
          betay = beta*yprel(j)
          do i=2,nxp-1
            q(i,j,k) = dx2fac*(  p(i,j-1,k) + p(i-1,j,k) + p(i+1,j,k)
     &                         + p(i,j+1,k) - 4.0d0*p(i,j,k) ) + betay
     &                - fnot*( aaa(k,k-1)*p(i,j,k-1) + aaa(k,k)*p(i,j,k)
     &                        +aaa(k,k+1)*p(i,j,k+1) )
          enddo
        enddo
!$OMP   END DO NOWAIT
      enddo
*     k = nl layer (bottom(ocean) or top(atmos.))
!$OMP DO SCHEDULE (STATIC)
      do j=2,nyp-1
        betay = beta*yprel(j)
        do i=2,nxp-1
          q(i,j,nl) = dx2fac*(  p(i,j-1,nl) + p(i-1,j,nl) + p(i+1,j,nl)
     &                        + p(i,j+1,nl) - 4.0d0*p(i,j,nl) ) + betay
     &               - fnot*(  aaa(nl,nl-1)*p(i,j,nl-1)
     &                       + aaa(nl, nl )*p(i,j, nl ) )
        enddo
      enddo
!$OMP END DO

*     k = kbot layer (bottom - for topography)
!$OMP DO SCHEDULE (STATIC)
      do j=2,nyp-1
        do i=2,nxp-1
          q(i,j,kbot) = q(i,j,kbot) + ddyn(i,j)
        enddo
      enddo
!$OMP END DO NOWAIT

!$OMP END PARALLEL

      END SUBROUTINE qcomp
c
c***********************************************************************
c
      SUBROUTINE merqcy (q, p, aaa, yprel, dxm2,
     &                   nxp, nyp, nl, ddyn, kbot)

*     Compute pot. vort. on meridional boundaries for a cyclic domain
*     Works for both atmosphere and ocean

*     Input arguments:
*     p     : 3-D array of dynamic pressure at p points
*     aaa   : A(nl,nl) matrix linking pressures and eta
*     yprel : y-coordinate of p points rel. to central latitude
*     dxm2  : dx^-2
*     nxp   : x-dimension of p array
*     nyp   : y-dimension of p array
*     nl    : no. of QG layers
*     ddyn  : dynamic (rescaled) topography
*     kbot  : index of layer containing topography

*     Output arguments:
*     q    : 3-D array of vorticity at p points

*     Modules
      USE parameters

      IMPLICIT NONE

*     Subroutine arguments
      integer nxp,nyp,nl,kbot
      double precision q(nxp,nyp,nl),p(nxp,nyp,nl),aaa(nl,nl),
     &                 yprel(nyp),dxm2,ddyn(nxp,nyp)
*
*     Local variables
      integer j,k
      double precision dx2fac,betay

      dx2fac = dxm2/fnot

!$OMP PARALLEL DEFAULT (NONE)
!$OMP&         PRIVATE (j,k,betay)
!$OMP&         SHARED  (nxp,nyp,nl,yprel,q,dx2fac,p,aaa,kbot,ddyn)

*     k = 1 layer (top(ocean) or bottom(atmos.))
!$OMP DO SCHEDULE (STATIC)
      do j=2,nyp-1
        betay = beta*yprel(j)
*       Western boundary (periodic)
        q( 1 ,j, 1 ) = dx2fac*
     &                  (  p(1,j-1,1) + p(nxp-1,j,1) + p(2,j,1)
     &                   + p(1,j+1,1) - 4.0d0*p(1,j,1) ) + betay
     &                 - fnot*( aaa(1,1)*p(1,j,1) + aaa(1,2)*p(1,j,2) )
*       Eastern boundary (periodic)
        q(nxp,j, 1 ) = q( 1 ,j, 1 )
      enddo
!$OMP END DO NOWAIT

*     Intermediate layers
      do k=2,nl-1
!$OMP   DO SCHEDULE (STATIC)
        do j=2,nyp-1
          betay = beta*yprel(j)
*         Western boundary (periodic)
          q( 1 ,j,k) = dx2fac*
     &                (  p( 1 ,j-1,k) + p(nxp-1,j,k) + p(2,j,k)
     &                 + p( 1 ,j+1,k) - 4.0d0*p( 1 ,j,k) ) + betay
     &                - fnot*( aaa(k,k-1)*p(1,j,k-1) + aaa(k,k)*p(1,j,k)
     &                        +aaa(k,k+1)*p(1,j,k+1) )
*         Eastern boundary (periodic)
          q(nxp,j,k) = q( 1 ,j,k)
        enddo
!$OMP   END DO NOWAIT
      enddo

*     k = nl layer (bottom(ocean) or top(atmos.))
!$OMP DO SCHEDULE (STATIC)
      do j=2,nyp-1
        betay = beta*yprel(j)
*       Western boundary (periodic)
        q( 1 ,j,nl) = dx2fac*
     &                (  p( 1 ,j-1,nl) + p(nxp-1,j,nl) + p(2,j,nl)
     &                 + p( 1 ,j+1,nl) - 4.0d0*p( 1 ,j,nl) ) + betay
     &               - fnot*(  aaa(nl,nl-1)*p(1,j,nl-1)
     &                       + aaa(nl, nl )*p(1,j, nl ) )
*       Eastern boundary (periodic)
        q(nxp,j,nl) = q( 1 ,j,nl)
      enddo
!$OMP END DO

*     k = kbot layer (bottom - for topography)
!$OMP DO SCHEDULE (STATIC)
      do j=2,nyp-1
        q( 1 ,j,kbot) = q( 1 ,j,kbot) + ddyn(1,j)
*       Eastern boundary (periodic)
        q(nxp,j,kbot) = q( 1 ,j,kbot)
      enddo
!$OMP END DO NOWAIT

!$OMP END PARALLEL

      END SUBROUTINE merqcy
c
c***********************************************************************
c
#ifndef atmos_only

      SUBROUTINE ocqbdy (qo, po)
*
*     Derive oceanic boundary pot. vorticity qo from dynamic pressure po
*     Equation (7.15), simplified because the tangential derivative p_tt
*     is known to be zero on all solid boundaries, and the normal derivative
*     p_nn is given by the mixed condition (2.27) (see also Appendix C)
*     Does all solid boundaries, i.e. zonal boundaries, and
*     also meridional boundaries if ocean is not cyclic.

*     Modules
      USE parameters
      USE occonst

      IMPLICIT NONE

*     Subroutine arguments
      double precision qo(nxpo,nypo,nlo),po(nxpo,nypo,nlo)
*
*     Local variables
      integer i,j,k
      double precision betays,betayn,f0Am,f0Ac,f0Ap,betay,bcfaco

*     Version with nondimensional bccooc
      bcfaco = bccooc*dxom2/( 0.5d0*bccooc + 1.0d0 )/fnot

*     Zonal boundaries (incl. end/corner pts) - mixed BCs
*     ===================================================
      betays = beta*yporel(  1 )
      betayn = beta*yporel(nypo)

!$OMP PARALLEL DEFAULT (NONE)
!$OMP&         PRIVATE (i,j,k,f0Am,f0Ac,f0Ap,betay)
!$OMP&         SHARED  (qo,po,bcfaco,betays,betayn,amatoc,yporel,ddynoc)

*     Top layer
*     ---------
      f0Ac = fnot*amatoc(1,1)
      f0Ap = fnot*amatoc(1,2)
!$OMP DO SCHEDULE (STATIC)
      do i=1,nxpo
        qo(i,  1 ,1) =  bcfaco*( po(i,   2  ,1) - po(i,  1 ,1) )
     &                - ( f0Ac*po(i,  1 ,1) + f0Ap*po(i,  1 ,2) )
     &                + betays
        qo(i,nypo,1) =  bcfaco*( po(i,nypo-1,1) - po(i,nypo,1) )
     &                - ( f0Ac*po(i,nypo,1) + f0Ap*po(i,nypo,2) )
     &                + betayn
      enddo
!$OMP END DO NOWAIT

*     Intermediate layers
*     -------------------
      do k=2,nlo-1
        f0Am = fnot*amatoc(k,k-1)
        f0Ac = fnot*amatoc(k, k )
        f0Ap = fnot*amatoc(k,k+1)
!$OMP   DO SCHEDULE (STATIC)
        do i=1,nxpo
          qo(i,  1 ,k) =  bcfaco*( po(i,   2  ,k) - po(i,  1 ,k) )
     &                  - (  f0Am*po(i,  1 ,k-1) + f0Ac*po(i,  1 ,k)
     &                     + f0Ap*po(i,  1 ,k+1) )
     &                  + betays
          qo(i,nypo,k) =  bcfaco*( po(i,nypo-1,k) - po(i,nypo,k) )
     &                  - (  f0Am*po(i,nypo,k-1) + f0Ac*po(i,nypo,k)
     &                     + f0Ap*po(i,nypo,k+1) )
     &                  + betayn
        enddo
!$OMP   END DO NOWAIT
      enddo

*     Bottom layer (includes topography)
*     ----------------------------------
      f0Am = fnot*amatoc(nlo,nlo-1)
      f0Ac = fnot*amatoc(nlo, nlo )
!$OMP DO SCHEDULE (STATIC)
      do i=1,nxpo
        qo(i,  1 ,nlo) =  bcfaco*( po(i,   2  ,nlo) - po(i,  1 ,nlo) )
     &                - ( f0Am*po(i,  1 ,nlo-1) + f0Ac*po(i,  1 ,nlo) )
     &                  + betays + ddynoc(i,  1 )
        qo(i,nypo,nlo) =  bcfaco*( po(i,nypo-1,nlo) - po(i,nypo,nlo) )
     &                - ( f0Am*po(i,nypo,nlo-1) + f0Ac*po(i,nypo,nlo) )
     &                  + betayn + ddynoc(i,nypo)
      enddo
!$OMP END DO NOWAIT

#  ifndef cyclic_ocean
*     Meridional (internal) boundaries - mixed BCs
*     ============================================
*     Top layer
*     ---------
      f0Ac = fnot*amatoc(1,1)
      f0Ap = fnot*amatoc(1,2)
!$OMP DO SCHEDULE (STATIC)
      do j=2,nypo-1
        betay = beta*yporel(j)
        qo(  1 ,j, 1 ) =  bcfaco*( po(   2  ,j,1) - po(  1 ,j,1) )
     &                  - ( f0Ac*po(  1 ,j,1) + f0Ap*po(  1 ,j,2) )
     &                  + betay
        qo(nxpo,j, 1 ) =  bcfaco*( po(nxpo-1,j,1) - po(nxpo,j,1) )
     &                  - ( f0Ac*po(nxpo,j,1) + f0Ap*po(nxpo,j,2) )
     &                  + betay
      enddo
!$OMP END DO NOWAIT

*     Intermediate layers
*     -------------------
      do k=2,nlo-1
        f0Am = fnot*amatoc(k,k-1)
        f0Ac = fnot*amatoc(k, k )
        f0Ap = fnot*amatoc(k,k+1)
!$OMP   DO SCHEDULE (STATIC)
        do j=2,nypo-1
          betay = beta*yporel(j)
          qo(  1 ,j,k) =  bcfaco*( po(   2  ,j,k) - po(  1 ,j,k) )
     &                  - (  f0Am*po(  1 ,j,k-1) + f0Ac*po(  1 ,j,k)
     &                     + f0Ap*po(  1 ,j,k+1) )
     &                  + betay
          qo(nxpo,j,k) =  bcfaco*( po(nxpo-1,j,k) - po(nxpo,j,k) )
     &                  - (  f0Am*po(nxpo,j,k-1) + f0Ac*po(nxpo,j,k)
     &                     + f0Ap*po(nxpo,j,k+1) )
     &                  + betay
        enddo
!$OMP   END DO NOWAIT
      enddo

*     Bottom layer (includes topography)
*     ----------------------------------
      f0Am = fnot*amatoc(nlo,nlo-1)
      f0Ac = fnot*amatoc(nlo, nlo )
!$OMP DO SCHEDULE (STATIC)
      do j=2,nypo-1
        betay = beta*yporel(j)
        qo(  1 ,j,nlo) =  bcfaco*( po(   2  ,j,nlo) - po(  1 ,j,nlo) )
     &                - ( f0Am*po(  1 ,j,nlo-1) + f0Ac*po(  1 ,j,nlo) )
     &                  + betay + ddynoc(  1 ,j)
        qo(nxpo,j,nlo) =  bcfaco*( po(nxpo-1,j,nlo) - po(nxpo,j,nlo) )
     &                - ( f0Am*po(nxpo,j,nlo-1) + f0Ac*po(nxpo,j,nlo) )
     &                  + betay + ddynoc(nxpo,j)
      enddo
!$OMP END DO NOWAIT
#  endif

!$OMP END PARALLEL

      END SUBROUTINE ocqbdy

#endif /* not atmos_only */
c
c***********************************************************************
c
#ifndef ocean_only

      SUBROUTINE atqzbd (qa, pa)

*     Derive atmospheric zonal boundary pot. vort. qa from pressure pa
*     Equation (7.15), simplified because the tangential derivative p_tt
*     is known to be zero on all solid boundaries, and the normal derivative
*     p_nn is given by the mixed condition (2.27) (see also Appendix C)

*     Modules
      USE parameters
      USE atconst

      IMPLICIT NONE

*     Subroutine arguments
      double precision qa(nxpa,nypa,nla),pa(nxpa,nypa,nla)
*
*     Local variables
      integer i,k
      double precision zbfaca,betays,betayn,f0Am,f0Ac,f0Ap

*     Version with nondimensional bccoat
      zbfaca = bccoat*dxam2/( 0.5d0*bccoat + 1.0d0 )/fnot

*     Zonal boundaries (incl ends) - mixed BCs
*     ========================================
      betays = beta*yparel(  1 )
      betayn = beta*yparel(nypa)

!$OMP PARALLEL DEFAULT (NONE)
!$OMP&         PRIVATE (i,k,f0Am,f0Ac,f0Ap)
!$OMP&         SHARED  (qa,pa,zbfaca,amatat,betays,betayn,ddynat)

*     Bottom layer (includes topography)
*     ----------------------------------
      f0Ac = fnot*amatat(1,1)
      f0Ap = fnot*amatat(1,2)
!$OMP DO SCHEDULE (STATIC)
      do i=1,nxpa
        qa(i,  1 ,1) =  zbfaca*( pa(i,   2  ,1) - pa(i,  1 ,1) )
     &                - ( f0Ac*pa(i,  1 ,1) + f0Ap*pa(i,  1 ,2) )
     &                + betays + ddynat(i,  1 )
        qa(i,nypa,1) =  zbfaca*( pa(i,nypa-1,1) - pa(i,nypa,1) )
     &                - ( f0Ac*pa(i,nypa,1) + f0Ap*pa(i,nypa,2) )
     &                + betayn + ddynat(i,nypa)
      enddo
!$OMP END DO NOWAIT

*     Intermediate layers
*     -------------------
      do k=2,nla-1
        f0Am = fnot*amatat(k,k-1)
        f0Ac = fnot*amatat(k, k )
        f0Ap = fnot*amatat(k,k+1)
!$OMP   DO SCHEDULE (STATIC)
        do i=1,nxpa
          qa(i,  1 ,k) =  zbfaca*( pa(i,   2  ,k) - pa(i,  1 ,k) )
     &                  - (  f0Am*pa(i,  1 ,k-1) + f0Ac*pa(i,  1 ,k)
     &                     + f0Ap*pa(i,  1 ,k+1) )
     &                  + betays
          qa(i,nypa,k) =  zbfaca*( pa(i,nypa-1,k) - pa(i,nypa,k) )
     &                  - (  f0Am*pa(i,nypa,k-1) + f0Ac*pa(i,nypa,k)
     &                     + f0Ap*pa(i,nypa,k+1) )
     &                  + betayn
        enddo
!$OMP   END DO NOWAIT
      enddo

*     Top layer
*     ---------
      f0Am = fnot*amatat(nla,nla-1)
      f0Ac = fnot*amatat(nla, nla )
!$OMP DO SCHEDULE (STATIC)
      do i=1,nxpa
        qa(i,  1 ,nla) =  zbfaca*( pa(i,   2  ,nla) - pa(i,  1 ,nla) )
     &                - ( f0Am*pa(i,  1 ,nla-1) + f0Ac*pa(i,  2 ,nla) )
     &                  + betays
        qa(i,nypa,nla) =  zbfaca*( pa(i,nypa-1,nla) - pa(i,nypa,nla) )
     &                - ( f0Am*pa(i,nypa,nla-1) + f0Ac*pa(i,nypa,nla) )
     &                  + betayn
      enddo
!$OMP END DO NOWAIT

!$OMP END PARALLEL

      END SUBROUTINE atqzbd

#endif /* not ocean_only */
c
c***********************************************************************
c
      END MODULE vorsubs
c
c***********************************************************************
