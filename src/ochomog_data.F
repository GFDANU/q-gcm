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
      MODULE ochomog

*     Contains the homogeneous solutions for the ocean,
*     and quantities used to ensure that constraints
*     on the overall solution are correctly satisfied

      USE parameters, ONLY : nxpo,nypo,nlo

      IMPLICIT NONE

      PUBLIC
      SAVE

#ifndef atmos_only

*     Variables and arrays containing the homogeneous
*     solutions themselves, and other time-invariant
*     quantities such as area and near-boundary
*     integrals used in the constraint equations,
*     and matrices used in the finite box case

#  ifdef cyclic_ocean
*     Zonally cyclic ocean
      double precision pch1oc(nypo,nlo-1),hc1soc(nlo-1),hc1noc(nlo-1),
     &                 pch2oc(nypo,nlo-1),hc2soc(nlo-1),hc2noc(nlo-1),
     &                 pbhoc(nypo),aipcho(nlo-1),hbsioc,aipbho
#  else
*     Finite box ocean
      double precision ochom(nxpo,nypo,nlo-1),aipohs(nlo-1),
     &                 cdiffo(nlo,nlo-1),cdhoc(nlo-1,nlo-1),
     &                 cdhlu(nlo-1,nlo-1)
      integer ipivch(nlo-1)
#  endif /* cyclic_ocean */

*     Variables and arrays containing time-varying
*     quantities such as area and (near-)boundary
*     integrals used in the constraint equations

      double precision xon(nlo-1),dpioc(nlo-1),dpiocp(nlo-1)
#  ifdef cyclic_ocean
*     Zonally cyclic ocean
      double precision ocncs(nlo),ocncn(nlo),ocncsp(nlo),ocncnp(nlo),
     &                 enisoc(nlo-1),eninoc(nlo-1),
     &                 ajisoc(nlo),ajinoc(nlo),
     &                 ap3soc(nlo),ap3noc(nlo),ap5soc(nlo),ap5noc(nlo),
     &                 txisoc,txinoc,bdrins,bdrinn
#  else
*     Finite box ocean
#  endif /* cyclic_ocean */

*     Zonally cyclic ocean
*     --------------------
*     pch1oc, pch2oc are the two homogeneous solutions for
*     each baroclinic mode, computed in subroutine homsol.
*     Each has a value of 1 on one zonal boundary,
*     and is exponentially small on the other.
*     They are functions of y only (tabulated at p points).
*     aipcho is the area integral of these solutions for each
*     mode (both solutions having the same area integral)
*
*     pbhoc is the homogeneous barotropic mode solution, computed in
*     subroutine homsol. There is only one barotropic homogeneous
*     mode, which is a function of y only. aipbho is its area integral
*
*     hc1soc, hc2soc, hc1noc, hc2noc are boundary integrals of
*     the homogeneous baroclinic solutions above, which arise
*     in the momentum constraint equations in a channel model
*
*     hbsioc is the inverse of a simpler equivalent of the
*     above quantities which occurs for the barotropic mode
*
*     ocncs, ocncn are vectors containing the quantities (one for each
*     layer) required for the momemtum constraints at the southern and
*     northern boundaries respectively, at the current time. ocncsp,
*     ocncnp are these quantities at the previous time level
*
*     enisoc, eninoc are the integrals of the entrainment along the
*     southern and northern boundaries respectively for each interface
*
*     txisoc, txinoc are the integrals of the oceanic windstress
*     component tauxo along the southern and northern boundaries of the
*     internal q domain (i.e. 1/2 gridlength in from physical boundaries)
*
*     ajisoc, ajinoc are the integrals of the Jacobian advection terms
*     in the ocean near the southern and northern boundaries (the
*     non-cancelling part of the Jacobian integral) for each layer
*
*     ap3soc, ap3noc are the integrals of ah2oc*d3p/dy3 along the strips in
*     the ocean near the southern and northern boundaries for each layer
*
*     ap5soc, ap5noc are the integrals of ah4oc*d5p/dy5 along the strips in
*     the ocean near the southern and northern boundaries for each layer
*
*     txisoc, txinoc are the integrals of the oceanic windstress
*     component tauxo along the southern and northern boundaries of the
*     internal q domain (i.e. 1/2 gridlength in from physical boundaries)
*
*     bdrins, bdrinn are integrals of boundary drag near zonal boundaries

*     Finite box ocean
*     ----------------
*     ochom contains the homogeneous baroclinic solutions
*
*     aipohs contains the area integrals of these solutions
*
*     cdiffo and cdhoc contain (constant) matrices used in the
*     mass constraint equation applied in subroutine ocinvq
*     cdiffo(m,k) contains the coefficient of the mode m
*     contribution to the interface k constraint equation
*
*     cdhlu contains the LU factorisation of cdhoc
*     ipivch contains the pivot indices for cdhlu
*
*     All the above are time-invariant,
*     and are computed in subroutine homsol
*
*     Constrain quantities for both configurations
*     --------------------------------------------
*     xon is the area integral of the entrainment across the
*     ocean interfaces; only xon(1) is nonzero (in the ocean
*     this is the only interface with any entrainment)
*
*     dpioc, dpiocp contain the current and previous values of the area
*     integral of the pressure difference across each internal interface.
*     This quantity (related to the interface displacement eta)
*     occurs in the mass continuity constraints in the ocean

#endif /* not atmos_only */

      END MODULE ochomog
c
c***********************************************************************
