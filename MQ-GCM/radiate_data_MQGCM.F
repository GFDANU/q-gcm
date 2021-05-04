c***********************************************************************
c     Q-GCM Version 1.5.0 : last modified 19/07/2012
c
c     Modified 08/15/2020 to include entrainment formulation for both
c     interfaces (Sergey Kravtsov)
c
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
     &                 aface(nla-1),bface,cface,dface,
     &                 aface1(nla-1),dface1,eface,eface1,
     &                 Emup,Eup(nla),f2,dfacel,dface1l,afacel(nla-1),
     &                 aface1l(nla-1),epsfac
      double precision Flambar,Fmebar,Febar(nla-1)
c-------------------------------------------------------------------------
c
c     added aface1 and dface1 for the entrainment radiation/diffusion coefficients
c     at the second interface; eface and eface1 are entrainment coefficients 
c     associated with o.m.l. radiation in atmospheric layers; Emup, Eup are radiation 
c     coefficients for SST contribution to radiation; f2 is the ratio of layer-2/3 to
c     layer-1/2 non-diffusive entrainment
c
c     Sergey Kravtsov 08/10/2020-08/16/2020
c
c     Added 4 versions of aface, dface and eface coefficients, representing 
c     various combinations of interior entrainment signs for the use with the
c     moist version of the model
c     
c     Sergey Kravtsov 10/21/2020
c
c------------------------------------------------------------------------
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
*     Emup = o.m.l. upward radiation coefft. in a.m.l.(W m^-2 K^-1)
*     Eup(k) = --------"-----------  in atmo. layer k (W m^-2 K^-1)
*     Dmdown = a.m.l. downward radiation coefft.      (W m^-2 K^-1)
*     rbetat(k) = multiplier of fs' for radiative balance
*     initialisation of interface displacement eta(k)
*     aface, bface, cface, dface, eface  = coefficients for calculating atmos.
*     entrainment across 1/2 from eta, etam, topography aTm', oTm', respectively
*     aface1, dface1, eface1 are coefficients for interface 2/3 
*     f2 is the fraction of e1' entrainment forcing 2/3 interface (e2'=f2*e1')
*     (where e1' and e2' is radiation/heat-exchange driven entrainment,
*      that is all entrainment less vertical diffusion)     
*
*     dfacel, dface1l, afacel, aface1l are analogous to dface and aface,
*     but act over land (08/16/2020)
*
*     the attenuation factor
*     epsfac=(1-zopt(3))*(1-zopt(2))*(1-zopt(3))
*
*
*     Flambar,Fmebar,Febar(nla-1) - climatological convective fluxes
*
c***********************************************************************
*
      END MODULE radiate
*
c***********************************************************************
