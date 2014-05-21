c***********************************************************************
c     Q-GCM Version 1.5.0 : last modified 27/07/2012
c***********************************************************************
c
c     Copyright 2012 Jeff Blundell, Andy Hogg and Bill Dewar.
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
      MODULE intsubs

*     Contains subroutines xintt and xintp for performing integrals
*     of fields tabulated at T- and p-gridpoints respectively

*     Modules

      IMPLICIT NONE

      PRIVATE

      PUBLIC :: xintt, xintp

      CONTAINS

c***********************************************************************
c
      SUBROUTINE xintt (xant, valt, nxt, nyt)
*
*     Computes area integral of a field valt tabulated at T points
*     Modified version with reduced accumulator error

      IMPLICIT NONE

*     Subroutine arguments
      integer nxt,nyt
      double precision xant,valt(nxt,nyt)
*
*     Local variables
      integer i,j
      double precision sumt,sumi
*
      sumt = 0.0d0
!$OMP PARALLEL DEFAULT (NONE)
!$OMP&         PRIVATE (i,j,sumi)
!$OMP&         SHARED  (nxt,nyt,sumt,valt)

!$OMP DO SCHEDULE (STATIC)
!$OMP&   REDUCTION(+:sumt)
      do j=1,nyt
        sumi = 0.0d0
        do i=1,nxt
          sumi = sumi + valt(i,j)
        enddo
        sumt = sumt + sumi
      enddo
!$OMP END DO NOWAIT

!$OMP END PARALLEL
      xant = sumt

      END SUBROUTINE xintt
c
c***********************************************************************
c
      SUBROUTINE xintp (xanp, valp, nxp, nyp)
*
*     Computes area integral of a field valp tabulated at p points
*     Includes weighted contribution from boundary points
*     Modified version with reduced accumulator error

      IMPLICIT NONE

*     Subroutine arguments
      integer nxp,nyp
      double precision xanp,valp(nxp,nyp)
*
*     Local variables
      integer i,j
      double precision sump,sumi,xxs,xxn
*
      sump = 0.0d0
      xxs = 0.5d0*valp(1, 1 )
      xxn = 0.5d0*valp(1,nyp)

!$OMP PARALLEL DEFAULT (NONE)
!$OMP&         PRIVATE (i,j,sumi)
!$OMP&         SHARED  (nxp,nyp,valp,sump,xxs,xxn)

*     Inner points + 0.5d0*( W & E boundaries)
!$OMP DO SCHEDULE (STATIC)
!$OMP&   REDUCTION(+:sump)
      do j=2,nyp-1
        sumi = 0.5d0*valp(1,j)
        do i=2,nxp-1
          sumi = sumi + valp(i,j)
        enddo
        sumi = sumi + 0.5d0*valp(nxp,j)
        sump = sump + sumi
      enddo
!$OMP END DO NOWAIT

*     N & S boundary contributions from inner
*     points (do corners separately because
*     they need an additional factor of 0.5)
!$OMP DO SCHEDULE (STATIC)
!$OMP&   REDUCTION(+:xxs,xxn)
      do i=2,nxp-1
        xxs = xxs + valp(i, 1 )
        xxn = xxn + valp(i,nyp)
      enddo
!$OMP END DO NOWAIT

!$OMP END PARALLEL

      xxs = xxs + 0.5d0*valp(nxp, 1 )
      xxn = xxn + 0.5d0*valp(nxp,nyp)

      xanp = sump + 0.5d0*( xxs + xxn )

      END SUBROUTINE xintp
c
c***********************************************************************
c
      END MODULE intsubs
c
c***********************************************************************
