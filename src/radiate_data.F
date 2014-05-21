c***********************************************************************
c     Q-GCM Version 1.5.0 : last modified 19/07/2012
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
      MODULE radiate

*     Contains radiation parameters (supplied and derived) needed during run

      USE parameters, ONLY : nla

      IMPLICIT NONE

      PUBLIC
      SAVE

      double precision fsbar,fspamp,fspco,zm,zopt(nla),gamma,xlamda,
     &                 rrcpat,rrcpoc,Aup(nla,nla-1),Adown(nla,nla-1),
     &                 Bmup,Bup(nla),B1down,Cmup,Cup(nla),C1down,
     &                 D0up,Dmup,Dup(nla),Dmdown,rbetat(nla-1),
     &                 aface(nla-1),bface,cface,dface
*
*     fsbar = background forcing amplitude (<0)       (W m^-2)
*     fspamp = perturb.  forcing amplitude (>0)       (W m^-2)
*     fspco = signed perturbation coefft (sign as f0) (W m^-2)
*     zm = optical depth in atmos. mixed layer        (m)
*     gamma = adiabatic lapse rate                    (K m^-1)
*     zopt(k) = optical depth in atmos. layer k       (m)
*     xlamda = (sensible + latent) transfer coefft.   (W m^-2 K^-1)
*     rrcpat = 1/(rho*cp) for atmosphere              (J^-1 m^3 K)
*     rrcpoc = 1/(rho*cp) for ocean                   (J^-1 m^3 K)
*     Coefficients for internal QG interface etas:
*     Aup(k,l) = layer k upward radiation coefft.     (W m^-3)
*     (multiplying eta of l-th interface)
*     Adown(k,l) = layer k downward radiation coefft. (W m^-3)
*     (multiplying eta of l-th interface)
*     Coefficients for atmos. m.l. eta:
*     Bmup = a.m.l.  upward radiation coefft.         (W m^-3)
*     Bup(k) = layer k upward radiation coefft.       (W m^-3)
*     B1down = layer 1 downward radiation coefft.     (W m^-3)
*     Coefficients for atmos. topography:
*     Cmup = a.m.l.  upward radiation coefft.         (W m^-3)
*     Cup(k) = layer k upward radiation coefft.       (W m^-3)
*     C1down = layer 1 downward radiation coefft.     (W m^-3)
*     Coefficients for m.l. temperature perturbation:
*     D0up = o.m.l. upward radiation coefft.          (W m^-2 K^-1)
*     Dmup = a.m.l. upward radiation coefft.          (W m^-2 K^-1)
*     Dup(k) = layer k upward radiation coefft.       (W m^-2 K^-1)
*     Dmdown = a.m.l. downward radiation coefft.      (W m^-2 K^-1)
*     rbetat(k) = multiplier of fs' for radiative balance
*     initialisation of interface displacement eta(k)
*     aface, bface, cface, dface = coefficients for calculating atmos.
*     entrainment (now assumed to be entirely across interface 1)
*     from eta, etam, topography and aTm' respectively
*
c***********************************************************************
*
      END MODULE radiate
*
c***********************************************************************
