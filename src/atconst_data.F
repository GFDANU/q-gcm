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
      MODULE atconst

*     Contains the parameters of atmospheric grid, layers,
*     timestepping and atmos. (zonal) boundary condition.
*     Also atmospheric vertical eigenmode details

      USE parameters, ONLY : nxpa, nypa, nxta, nyta, nla

      IMPLICIT NONE

      PUBLIC
      SAVE

      double precision dxa,dya,hdxam1,dxam2,rdxaf0,xla,yla,
     &                 xpa(nxpa),ypa(nypa),yparel(nypa),
     &                 xta(nxta),yta(nyta),ytarel(nyta),
     &                 rhoat,cpat,gpat(nla-1),hat(nla),hta,tabsat(nla),
     &                 tmbara,tat(nla),ah4at(nla),dta,tdta,bccoat,
     &                 xc1ast(nxta,nyta),xcexp,amatat(nla,nla),
     &                 cphsat(nla),rdefat(nla),rdm2at(nla),
     &                 ctl2mat(nla,nla),ctm2lat(nla,nla),
     &                 dtopat(nxpa,nypa),davgat,ddynat(nxpa,nypa)
*
*     dxa, dya are the atmospheric gridlengths (m)
*     hdxam1 = 0.5d0/dxa      (m^-1)
*     dxam2 = 1.0d0/(dxa*dxa) (m^-2)
*     rdxaf0 = 1.0d0/(dxa*fnot) (m^-1 s)
*     xla = x-range of atmos. domain (m)
*     yla = y-range of atmos. domain (m)
*     xpa = vector of p-gridpoint positions in atmosphere (m)
*     ypa = vector of p-gridpoint positions in atmosphere (m)
*     xta = vector of T-gridpoint positions in atmosphere (m)
*     yta = vector of T-gridpoint positions in atmosphere (m)
*     yparel, ytarel are gridpoint positions relative to midlatitude
*
*     rhoat = atmospheric density (kg m^-3)
*     cpat = atmospheric specific heat capacity (J kg^-1 K^-1)
*     gpat = reduced gravity across atmospheric interfaces (m s^-2)
*     hat = unperturbed atmospheric layer thicknesses (m)
*     hta = total atmos. depth = hat(1) + ... + hat(nla) (m)
*     tabsat is the absolute (potential) temp. of each atmos. layer (K)
*     tmbara is the mean atmos. mixed layer absolute temperature (K)
*     tat is the temperature anomaly of each atmospheric layer,
*     relative to the mean state radiative equilibrium (K)
*     ah4at = Del-4th damping coefft for atmosphere (m^4 s^-1)
*     dta = atmos. timestep (s)
*     tdta = 2.0d0*dta (s)
*
*     bccoat = mixed boundary coefft. for atmosphere (nondim.)
*     bccoat = 0  =>  free slip/no stress
*     bccoat -> Inf  =>  no slip
*     xc1ast = an array to store info for xcexp runs
*     xcexp = coupling coefft. x
*
*     amatat is the atmospheric A matrix linking pressures and eta
*     cphsat = phase speeds of atmospheric modes (m s^-1)
*     rdefat = deformation radii of atmospheric modes (m)
*     rdm2at = 1.0/(deformation radii of atmos. modes)**2 (m^-2)
*     cl2mat = coefficients for conversion from layers to modes
*     ctl2mat = transpose of matrix for conversion from layers to modes
*     cm2lat = coefficients for conversion from modes to layers
*     ctm2lat = transpose of matrix for conversion from modes to layers
*     (transposes are held for reasons of computational efficiency)
*     (All the above are got by solving an eigenvalue/vector
*      equation using subroutine eigmod, called from the main program)
*
*     dtopat is the (nxpa,nypa) array of topography at the
*     bottom of the atmosphere, tabulated at p points (m)
*     davgat is its average value (m)
*     ddynat is the dynamic (rescaled) topography (s^-1)
*
      END MODULE atconst
c
c***********************************************************************
