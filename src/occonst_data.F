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
      MODULE occonst

*     Contains the parameters of oceanic grid, layers,
*     timestepping and oceanic boundary condition.
*     Also oceanic vertical eigenmode details

      USE parameters, ONLY : nxpo, nypo, nxto, nyto, nlo

      IMPLICIT NONE

      PUBLIC
      SAVE

      double precision dxo,dyo,hdxom1,dxom2,rdxof0,xlo,ylo,
     &                 xpo(nxpo),ypo(nypo),yporel(nypo),
     &                 xto(nxto),yto(nyto),ytorel(nyto),
     &                 rhooc,cpoc,gpoc(nlo-1),hoc(nlo),hto,tabsoc(nlo),
     &                 tmbaro,toc(nlo),ah2oc(nlo),ah4oc(nlo),delek,
     &                 dto,tdto,bccooc,ycexp,amatoc(nlo,nlo),
     &                 cphsoc(nlo),rdefoc(nlo),rdm2oc(nlo),
     &                 ctl2moc(nlo,nlo),ctm2loc(nlo,nlo),
     &                 dtopoc(nxpo,nypo),davgoc,ddynoc(nxpo,nypo)
*
*     dxo, dyo are the oceanic gridlengths (m)
*     hdxom1 = 0.5d0/dxo      (m^-1)
*     dxom2 = 1.0d0/(dxo*dxo) (m^-2)
*     rdxof0 = 1.0d0/(dxo*fnot) (m^-1 s)
*     xlo = x-range of ocean domain (m)
*     ylo = y-range of ocean domain (m)
*     xpo = vector of p-gridpoint positions in ocean
*           (including offset w.r.t. atmosphere) (m)
*     ypo = vector of p-gridpoint positions in ocean
*           (including offset w.r.t. atmosphere) (m)
*     xto = vector of T-gridpoint positions in ocean
*           (including offset w.r.t. atmosphere) (m)
*     yto = vector of T-gridpoint positions in ocean
*           (including offset w.r.t. atmosphere) (m)
*     yporel, ytorel are gridpoint positions relative to midlatitude
*
*     rhooc = oceanic density (kg m^-3)
*     cpoc = oceanic specific heat capacity (J kg^-1 K^-1)
*     gpoc = reduced gravity across oceanic interfaces (m s^-2)
*     hoc = unperturbed oceanic layer thicknesses (m)
*     hto = total ocean depth = hoc(1) + ... + hoc(nlo) (m)
*     tabsoc is the absolute (potential) temp. of each ocean layer (K)
*     tmbaro is the mean ocean mixed layer absolute temperature (K)
*     toc is the temperature anomaly of each oceanic layer,
*     relative to the mean state radiative equilibrium (K)
*     ah2oc = Del-sqd damping coefft for ocean (m^2 s^-1)
*     ah4oc = Del-4th damping coefft for ocean (m^4 s^-1)
*     delek = bottom Ekman layer thickness (m)
*     dto = ocean timestep (s)
*     tdto = 2.0d0*dto (s)
*
*     bccooc = mixed boundary coefft. for ocean (nondim.)
*     bccooc = 0  =>  free slip/no stress
*     bccooc -> Inf  =>  no slip
*
*     ycexp = coupling coefft. y
*
*     amatoc is the oceanic A matrix linking pressures and eta
*     cphsoc = phase speeds of oceanic modes (m s^-1)
*     rdefoc = deformation radii of oceanic modes (m)
*     rdm2oc = 1.0/(deformation radii of oceanic modes)**2 (m^-2)
*     cl2moc = coefficients for conversion from layers to modes
*     ctl2moc = transpose of matrix for conversion from layers to modes
*     cm2loc = coefficients for conversion from modes to layers
*     ctm2loc = transpose of matrix for conversion from modes to layers
*     (transposes are held for reasons of computational efficiency)
*     (All the above are got by solving an eigenvalue/vector
*      equation using subroutine eigmod, called from the main program)
*
*     dtopoc is the (nxo,nyo) array of topography at the
*     bottom of the ocean, tabulated at p points (m)
*     davgoc is its average value (m)
*     ddynoc is the dynamic (rescaled) topography (s^-1)
*
      END MODULE occonst
c
c***********************************************************************
