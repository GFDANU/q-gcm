c***********************************************************************
c     Q-GCM Version 1.5.0 : last modified 21/07/2013
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
      MODULE athomog

*     Contains the homogeneous solutions for the atmosphere,
*     and quantities used to ensure that constraints
*     on the overall solution are correctly satisfied

      USE parameters, ONLY : nypa,nla

      IMPLICIT NONE

      PUBLIC
      SAVE

#ifndef ocean_only

*     Variables and arrays containing the homogeneous
*     solutions themselves, and other time-invariant
*     quantities such as area and near-boundary
*     integrals used in the constraint equations

      double precision pch1at(nypa,nla-1),hc1sat(nla-1),hc1nat(nla-1),
     &                 pch2at(nypa,nla-1),hc2sat(nla-1),hc2nat(nla-1),
     &                 pbhat(nypa),aipcha(nla-1),hbsiat,aipbha

*     Variables and arrays containing time-varying
*     quantities such as area and (near-)boundary
*     integrals used in the constraint equations

      double precision xan(nla-1),dpiat(nla-1),dpiatp(nla-1)
      double precision atmcs(nla),atmcn(nla),atmcsp(nla),atmcnp(nla),
     &                 enisat(nla-1),eninat(nla-1),
     &                 ajisat(nla),ajinat(nla),ap5sat(nla),ap5nat(nla),
     &                 txisat,txinat

*     pch1at, pch2at are the two homogeneous solutions for
*     each baroclinic mode, computed in subroutine homsol.
*     Each has a value of 1 on one zonal boundary,
*     and is exponentially small on the other.
*     They are functions of y only (tabulated at p points).
*     aipcha is the area integral of these solutions for each
*     mode (both solutions having the same area integral)
*
*     pbhat is the homogeneous barotropic mode solution, computed in
*     subroutine homsol. There is only one barotropic homogeneous
*     mode, which is a function of y only. aipbha is its area integral
*
*     hc1sat, hc2sat, hc1nat, hc2nat are boundary integrals of
*     the homogeneous baroclinic solutions above, which arise
*     in the momentum constraint equations in a channel model
*
*     hbsiat is the inverse of a simpler equivalent of the
*     above quantities which occurs for the barotropic mode
*
*     xan is the area integral of the entrainment across each interface
*
*     dpiat, dpiatp contain the current and previous values of the area
*     integral of the pressure difference across each internal interface.
*     This quantity (related to the interface displacement eta)
*     occurs in the mass continuity constraints in the atmosphere
*
*     atmcs, atmcn are vectors containing the quantities (one for each
*     layer) required for the momemtum constraints at the southern and
*     northern boundaries respectively, at the current time. atmcsp,
*     atmcnp are these quantities at the previous time level
*
*     enisat, eninat are the integrals of the entrainment along the
*     southern and northern boundaries respectively for each interface
*
*     ajisat, ajinat are the integrals of the Jacobian advection terms
*     in the atmosphere near the southern and northern boundaries (the
*     non-cancelling part of the Jacobian integral) for each layer
*
*     ap5sat, ap5nat are the integrals of ah4at*d5p/dy5 along the strips in
*     the atmosphere near the southern and northern boundaries for each layer
*
*     txisat, txinat are the integrals of the atmospheric windstress
*     component tauxa along the southern and northern boundaries of the
*     internal q domain (i.e. 1/2 gridlength in from physical boundaries)

#endif /* not ocean_only */

      END MODULE athomog
c
c***********************************************************************
