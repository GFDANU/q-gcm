c***********************************************************************
c     Q-GCM Version 1.5.0 : last modified 20/07/2012
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
      MODULE intrfac

*     Storage for mixed layer parameters and
*     variables, for both atmosphere and ocean.
*     Includes dynamic stresses (tau/rho), the horizontal
*     Ekman velocity components in the atmospheric mixed layer,
*     and averaged stresses along the sides of atmos. p-cells.

      USE parameters, ONLY : nxpa, nypa, nxta, nyta,
     &                       nxpo, nypo, nxto, nyto

      IMPLICIT NONE

      PUBLIC
      SAVE

      double precision sst(nxto,nyto),sstm(nxto,nyto),sstbar(nyto),
     &                 ast(nxta,nyta),astm(nxta,nyta),astbar(nyta),
     &                 tauxo(nxpo,nypo),tauyo(nxpo,nypo),
     &                 tauxa(nxpa,nypa),tauya(nxpa,nypa),
     &                 uekat(nxpa,nyta),vekat(nxta,nypa),
     &                 txatpb(nxpa,0:nypa),tyatpb(0:nxpa,nypa),
     &                 fnetoc(nxto,nyto),fnetat(nxta,nyta),
     &                 hmixa(nxta,nyta),hmixam(nxta,nyta),hmadmp,
     &                 hmoc,hmat,hmamin,cdat,raoro,st2d,st4d,
     &                 at2d,at4d,ahmd,tsbdy,tnbdy
*
*     sst  is the current  sea-surface temperature (anomaly) (K)
*     sstm is the previous sea-surface temperature (anomaly) (K)
*     sstbar is the radiative equilibrium SST anomaly (K)
*     (a function of latitude only)
*     ast  is the current  atmos. m. l. temperature (anomaly) (K)
*     astm is the previous atmos. m. l. temperature (anomaly) (K)
*     astbar is the radiative equilibrium a.m.l. temp. anomaly (K)
*     (a function of latitude only)
*
*     tauxo, tauyo are dynamic stress components in oceanic m. l. (m^2 s^-2)
*     tauxa, tauya are dynamic stress components in atmos.  m. l. (m^2 s^-2)
*     (stresses are tabulated at p points)
*     uekat, vekat are the Ekman velocity components in the
*     atmos. mixed layer, tabulated at centres of T-cell faces,
*     for use in advection terms in subroutine amladf (m s^-1)
*     txatpb, tyatpb are dynamic stress components averaged
*     along the sides of atmos p-cells, for use in efficiently
*     computing wekpa within xforc, and also enabling accurate
*     computation of the atmos. wind work (u.tau) term (m^2 s^-2)
*
*     fnetoc is the net oceanic mixed layer forcing at T pts (W m^-2)
*     fnetat is the net atmos.  mixed layer forcing at T pts (W m^-2)
*     hmixa  is the current  thickness of the atmos. m. l. (m)
*     hmixam is the previous thickness of the atmos. m. l. (m)
*     (mixed layer thickness is tabulated at T points)
*     hmadmp is the atmospheric m. l. damping coefft. (W m^-3)
*     hmoc is the (fixed) oceanic mixed layer depth (m)
*     hmat is the (fixed) atmospheric mixed layer depth to which
*     hmixa is initialised, and towards which it relaxes. (m)
*     hmamin is the minimum permitted atmos. m. l. thickness (m)
*     cdat = Quadratic drag coefft. (nondim.)
*     raoro = density ratio (atmospheric/oceanic) (nondim.)
*
*     st2d is the horizontal sst Del-sqd diffusivity (m^2 s^-1)
*     st4d is the horizontal sst Del-4th diffusivity (m^4 s^-1)
*     at2d is the horizontal ast Del-sqd diffusivity (m^2 s^-1)
*     at4d is the horizontal ast Del-4th diffusivity (m^4 s^-1)
*     ahmd is the horizontal atmos. m. l. diffusivity (m^2 s^-1)
*
*     tsbdy is the (relative) southern boundary temp-
*     erature specified for the oceanic mixed later (K),
*     only active if the sb_hflux option is activated.
*     tnbdy is the (relative) northern boundary temp-
*     erature specified for the oceanic mixed later (K),
*     only active if the nb_hflux option is activated.
*
      END MODULE intrfac
c
c***********************************************************************
