c***********************************************************************
c     Q-GCM Version 1.5.0 : last modified 25/07/2012
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
      MODULE atstate

*     Contains the atmospheric state (i.e. dynamic pressure and
*     vorticity) at the current and previous time levels, the
*     entrainment between the atmospheric layers, and the vertical
*     Ekman velocity in the atmosphere at both p and T gridpoints.

      USE parameters, ONLY : nxpa, nypa, nxta, nyta, nla

      IMPLICIT NONE

      PUBLIC
      SAVE

      double precision pa(nxpa,nypa,nla),pam(nxpa,nypa,nla),
     &                 qa(nxpa,nypa,nla),qam(nxpa,nypa,nla),
     &                 entat(nxpa,nypa),entat1(nxpa,nypa),
     &                 wekpa(nxpa,nypa),wekta(nxta,nyta),
     &                 precip_i(nxta,nyta,nla),
     &                 hum(nxta,nyta,nla),humm(nxta,nyta,nla)
*
*     pa  is the dynamic pressure in each layer at the current  time (m^2 s^-1)
*     pam is the dynamic pressure in each layer at the previous time (m^2 s^-1)
*     qa  is the vorticity in each layer at the current  time (s^-1)
*     qam is the vorticity in each layer at the previous time (s^-1)
*     entat is the entrainment between the atmos. layers (at p points) (m s^-1)
*     We are assuming all entrainment is across interface 1 (not any longer! SK: 07/22/2020
*     wekpa is the atmos. Ekman velocity at p gridpoints (m s^-1)
*     wekta is the atmos. Ekman velocity at T gridpoints (m s^-1)
*
*     entat0 is the lower atmo layer latent heat release in the moist model
*     (added by SK 10/22/2020)
*
*     hum, humm are interior layers' specific humidities (kg/kg)
*
*
      END MODULE atstate
c
c***********************************************************************

