c***********************************************************************
c     Q-GCM Version 1.5.0 : last modified 21/06/2013
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
      MODULE ocstate

*     Contains the oceanic state (i.e. dynamic pressure and
*     vorticity) at the current and previous time levels, the
*     entrainment between the oceanic layers, and the vertical
*     Ekman velocity in the ocean at both p and T gridpoints.

      USE parameters, ONLY : nxpo, nypo, nxto, nyto, nlo

      IMPLICIT NONE

      PUBLIC
      SAVE

#ifndef atmos_only

      double precision po(nxpo,nypo,nlo),pom(nxpo,nypo,nlo),
     &                 qo(nxpo,nypo,nlo),qom(nxpo,nypo,nlo),
     &                 entoc(nxpo,nypo),wekpo(nxpo,nypo),
     &                 wekto(nxto,nyto)

*     po  is the dynamic pressure in each layer at the current  time (m^2 s^-2)
*     pom is the dynamic pressure in each layer at the previous time (m^2 s^-2)
*     qo  is the vorticity in each layer at the current  time (s^-1)
*     qom is the vorticity in each layer at the previous time (s^-1)
*     entoc is the entrainment between the ocean layers (at p points) (m s^-1)
*     wekpo is the ocean Ekman velocity at p gridpoints (m s^-1)
*     wekto is the ocean Ekman velocity at T gridpoints (m s^-1)

#endif /* not atmos_only */

      END MODULE ocstate
c
c***********************************************************************
