***********************************************************************
c     Q-GCM Version 1.5.0 : last modified 31/12/2013
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
      MODULE monitor

*     Storage for diagnostic quantities for monitoring/debugging.
*     These quantities are mainly computed by subroutine monnc_comp
*     in monitor_diag.F, but also at various other points in the
*     code, for technical reasons. These values are output as time
*     series, by calling subroutine monnc_out from the main program.
*     Output interval is nocmon atmospheric timesteps = dgnday days
*     This replaces the former common block monitor.cmn

      USE parameters, ONLY : nla, nlo

      IMPLICIT NONE

      PUBLIC
      SAVE

#ifndef ocean_only
      integer atstpos(nla)
      double precision entmat(nla-1),enamat(nla-1),etamat(nla-1),
     &                 ermasa(nla-1),emfrat(nla-1),kealat(nla),
     &                 pavgat(nla),qavgat(nla),ah4dat(nla),
     &                 ddtkeat(nla),ddtpeat(nla-1),et2mat(nla-1),
     &                 ugminat(nla),ugmaxat(nla),
     &                 vgminat(nla),vgmaxat(nla),cnqgat(nla)
#endif /* not ocean_only */

#ifndef atmos_only
      integer ocjpos(nlo)
      double precision entmoc,enamoc,etamoc(nlo-1),
     &                 osfmin(nlo),osfmax(nlo),occirc(nlo),occtot,
     &                 ermaso(nlo-1),emfroc(nlo-1),kealoc(nlo),
     &                 pavgoc(nlo),qavgoc(nlo),ah2doc(nlo),ah4doc(nlo),
     &                 ddtkeoc(nlo),ddtpeoc(nlo-1),et2moc(nlo-1),
     &                 ugminoc(nlo),ugmaxoc(nlo),
     &                 vgminoc(nlo),vgmaxoc(nlo),cnqgoc(nlo)
#endif /* not atmos_only */

      double precision wetmoc,wetmat,wepmoc,wepmat,
     &                 watmoc,watmat,wapmoc,wapmat,
     &                 tmlmoc,sstmin,sstmax,tmlmat,
     &                 astmin,astmax,tmaooc,ttmads,
     &                 vfmads,ttmdfs,ttmadn,vfmadn,ttmdfn,slhfav,
     &                 oradav,arocav,arlaav,cfraoc,cfraat,centoc,
     &                 hmlmat,centat,hfmloc,hcmlat,
     &                 pkenoc,pkenat(nla-1),utauoc,utauat,
     &                 btdgoc,olrtop,atstval(nla),ocjval(nlo),
     &                 umminoc,ummaxoc,vmminoc,vmmaxoc,cnmloc,
     &                 umminat,ummaxat,vmminat,vmmaxat,cnmlat
*
*     wetmoc, wetmat are mean Ekman velocities at T points (should be zero)
*     wepmoc, wepmat are mean Ekman velocities at p points (should be zero)
*     watmoc, watmat are mean modulus of Ekman velocities at T points
*     wapmoc, wapmat are mean modulus of Ekman velocities at p points
*     All the above are in (m s^-1).
*     Computed in monnc_comp
*
*     entmat, enamat are the mean values of entrainment and
*     |entrainment| for the atmosphere; entmoc, enamoc are the
*     corresponding quantities for the ocean. All in (m s^-1).
*     Computed in monnc_comp
*
*     etamoc, etamat are mean values of eta at p points. (m)
*     Computed in monnc_comp
*
*     tmlmoc, tmlmat are mean mixed layer (rel.) temperatures
*     tmaooc is the mean atmos. m.l. (rel) temp over the ocean
*     hmlmat is mean mixed layer thickness (m)
*     Computed in monnc_comp
*
*     sstmin, sstmax are min. and max.
*     ocean mixed layer (rel.) temperatures
*     Computed in monnc_comp
*
*     astmin, astmax are min. and max.
*     atmos. mixed layer (rel.) temperatures
*     Computed in monnc_comp
*
*     ttmads, ttmdfs are the mean advective and diffusive
*     contributions to dT/dt in the oceanic mixed layer (K s^-1)
*     vfmads is the mean advective velocity into the oceanic mixed
*     layer. The above are produced by the modified southern boundary
*     condition. ttmadn, ttmdfn and vfmadn are the corresponding
*     quantities for the modified northern boundary condition.
*     Quantities are zero if the modified b.c. is not activated.
*     Computed in omladf
*
*     slhfav, oradav are heat fluxes averaged over ocean (W m^-2)
*     arocav, arlaav are atmospheric radiative fluxes (upwards +ve)
*     averaged over ocean and land respectively (W m^-2)
*     All computed in xforc, so no duplication.
*
*     cfraoc, cfraat are the fractions of cells
*     convecting in the ocean and the atmosphere.
*     Computed in oml and aml
*
*     centoc, centat are the (area integrals of)
*     convective entrainment in the ocean and atmosphere.
*     Computed in oml and aml
*
*     ermasa, emfrat are the absolute and fractional
*     mass errors at each atmospheric interface
*     Computed in atinvq
*
*     ermaso, emfroc are the absolute and fractional
*     mass errors at each oceanic interface
*     Computed in ocinvq
*
*     kealoc(k), kealat(k) are average KE of each layer (J m^-2)
*     = 0.5*rho0*H(k) * Integ [ ug(k)**2 + vg(k)**2 ] dA / Area
*     where ug, vg are geostrophic velocity components, and no
*     distinction is made for the mixed layer embedded in layer 1
*     Computed in monnc_comp
*
*     pavgoc, pavgat are average dynamic pressure of each layer (m^2 s^-2)
*     Computed in monnc_comp
*
*     qavgoc, qavgat are average vorticity of each layer (s^-1)
*     Computed in monnc_comp
*
*     hfmloc is the heat flux at the bottom of the ocean m.l. (W m^-2)
*     = Integ (oTm'*Wekman) dA / Area
*     Computed in monnc_comp
*
*     hcmlat is the total heat content of the atmos. m.l. (J m^-2)
*     = Integ (aTm'*ahm) dA / Area
*     Computed in monnc_comp
*
*     ah2doc(k) is the area-averaged 2nd order dissipation in each layer (W m^-2)
*     = ( rho0*Ah(k)*H(k) )
*      * Integ ( u(k)*Del-sqd(u(k)) + v(k)*Del-sqd(v(k)) ) dA / Area
*     Computed in monnc_comp
*
*     ah4doc(k), ah4dat(k) are area-averaged 4th-order dissipation in each layer (W m^-2)
*     = ( rho0*Ah(k)*H(k) )
*      * Integ ( u(k)*Del-4th(u(k)) + v(k)*Del-4th(v(k)) ) dA / Area
*     Computed in monnc_comp
*
*     ddtkeoc(k), ddtkeat(k) are the instantaneous rate of change
*     of the (area averaged) KE in each layer (W m^-2)
*     = rho0*H(k) * Integ ( ug(k)*d/dt[ug(k)] + vg(k)*d/dt[vg(k)] ) dA / Area
*     where ug, vg are geostrophic velocity components, and no
*     distinction is made for the mixed layer embedded in layer 1
*     Computed in monnc_comp
*
*     ddtpeoc(k), ddtpeat(k) are the instantaneous rate of change
*     of the (area averaged) PE at each interface (W m^-2)
*     = rho0*gpr(k)*[ Integ ( eta(k)*d/dt[eta(k)] ) dA /Area ]
*     Computed in monnc_comp
*
*     pkenoc, pkenat are part of the KE exchanges between layers (W m^-2)
*     = rho0*gpr(k)*Integ ( eta1*entrainment ) dA /Area
*     Computed in monnc_comp
*
*     utauoc, utauat are the KE exchange between ocean & atmosphere (W m^-2)
*     = rho*Integ ( v1*tauy - u1*taux ) dA / Area
*     Computed in monnc_comp
*
*     btdgoc is the mean energy dissipation due
*     to drag in the bottom Ekman layer (W m^-2)
*     = 0.5*rhooc*delek*|f0|*Integ ( u^2(nlo) + v^2(nlo) ) dA /Area
*     Computed in monnc_comp
*
*     et2moc, et2mat are the mean values of eta^2 at each interface
*     Computed in monnc_comp (usually) and constr (initially)
*
*     osfmin, osfmax are the min and max values of the ocean
*     volume transport stream function (in Sverdrups) in each
*     layer, where stream function psi = hoc(k)*(p - pbdy)/f0
*     Sign convention is the geophysical one, as in Gill 1982,
*     which is the opposite of usual fluid dynamics convention.
*     Here u = -d(psi)/dy, v = d(psi)/dx
*     Computed in monnc_comp
*
*     occirc is the zonal volume transport (in Sverdrups)
*     in each layer, and occtot is the total over all layers.
*     Only non-zero in the cyclic ocean case
*     Computed in monnc_comp
*
*     olrtop is the mean outgoing longwave radiation
*     perturbation at the top of the atmosphere (W m^-2)
*     Computed in monnc_comp
*
*     umminoc and ummaxoc are the minimum and maximum u-components
*     of the mixed layer velocity on the ocean, and ugminoc(k) and
*     ugmaxoc(k) are the corresponding values for the geostrophic
*     flow in the oceanic Q-G layers. The v-components are similarly
*     computed, as are all the equivalent atmospheric values.
*     Computed in subroutines courat, couroc in monitor_diag.F
*
*     cnqgoc and cnmloc are the maximum Courant numbers in the
*     oceanic Q-G layers and the oceanic mixed layer respectively.
*     cnqgat and cnmlat are the atmospheric equivalents.
*     Computed in subroutines courat, couroc in monitor_diag.F
*
      END MODULE monitor
c
c***********************************************************************
