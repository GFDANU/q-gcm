c***********************************************************************
c     Q-GCM Version 1.5.0 : last modified 02/01/2014
c***********************************************************************
c
c     Copyright 2014 Jeff Blundell, Andy Hogg and Bill Dewar.
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
      MODULE mondiag

*     Contains subprograms for computing and outputting (in netCDF
*     format) the instantaneous values of various diagnostic quantities
*     required for monitoring/debugging the Q-GCM climate model.
*     A few quantities (e.g. convection) need to be computed
*     as part of the actual timestepping; otherwise they are
*     all computed in subroutine monnc_comp for efficiency.
*     Output interval is nocmon atmospheric timesteps = dgnday days
*     (controlled by when monnc_comp is called from the main program).
*     This is intended to be a short interval; many time samples can
*     be afforded because all the output fields are low dimensional.

      IMPLICIT NONE

      PRIVATE

      PUBLIC  :: monnc_comp, monnc_init, monnc_out
#  ifndef cyclic_ocean
      PRIVATE :: del4bx
#  endif
      PRIVATE :: del4ch, genint
#  ifndef ocean_only
      PRIVATE :: courat
#  endif
#  ifndef atmos_only
      PRIVATE :: couroc
#  endif

      integer, PUBLIC, SAVE :: monncid

#ifdef use_netcdf
*     netCDF identifiers for monitoring variables
      integer, PRIVATE, SAVE :: time_id,
     &                  wetmoc_id, wetmat_id, wepmoc_id, wepmat_id,
     &                  watmoc_id, watmat_id, wapmoc_id, wapmat_id,
     &                  entmoc_id, entmat_id, enamoc_id, enamat_id,
     &                  etamoc_id, etamat_id, sstmin_id, sstmax_id,
     &                  tmlmoc_id, astmin_id, astmax_id, tmlmat_id,
     &                  tmaooc_id, ttmads_id, ttmdfs_id, vfmads_id,
     &                  ttmadn_id, ttmdfn_id, vfmadn_id,
     &                  slhfav_id, oradav_id, arocav_id, arlaav_id,
     &                  cfraoc_id, cfraat_id, centoc_id, centat_id,
     &                  ermasa_id, emfrat_id, hmlmat_id, kealoc_id,
     &                  kealat_id, pavgoc_id, pavgat_id, qavgoc_id,
     &                  qavgat_id, hfmloc_id, hcmlat_id, ermaso_id,
     &                  ah2doc_id, ah4doc_id, ah4dat_id, emfroc_id,
     &                  ddtkeoc_id, ddtkeat_id, ddtpeoc_id, ddtpeat_id,
     &                  pkenoc_id, pkenat_id, utauoc_id, utauat_id,
     &                  btdgoc_id, et2mat_id, et2moc_id, atstval_id,
     &                  atstpos_id, ocjval_id, ocjpos_id, osfmin_id,
     &                  osfmax_id, occirc_id, occtot_id, olrtop_id
*     New identifiers for velocity component extrema and Courant numbers
      integer, PRIVATE, SAVE ::
     &                  umminoc_id, ummaxoc_id, vmminoc_id, vmmaxoc_id,
     &                  ugminoc_id, ugmaxoc_id, vgminoc_id, vgmaxoc_id,
     &                  umminat_id, ummaxat_id, vmminat_id, vmmaxat_id,
     &                  ugminat_id, ugmaxat_id, vgminat_id, vgmaxat_id,
     &                  cnmloc_id, cnqgoc_id, cnmlat_id, cnqgat_id
*     (Replaces the former common block monnc.cmn)
#endif /* use_netcdf */

      CONTAINS

c***********************************************************************
*
      SUBROUTINE monnc_comp

*     Computes as many as possible of the diagnostic quantities
*     required for monitoring/debugging the QG climate model.
*     A few quantities (e.g. convection) need to be computed
*     as part of the actual timestepping; otherwise try to do
*     all calculations in this routine for efficiency.
*     The routine should be called after xforc; at this point
*     in the timestepping all quantities are consistent.
*     The results are all written to the data module monitor,
*     defined in monitor_data.F .
*
*     N.B. use included integration routine genint for all area
*     integrations, rather than the familiar routines xintt and
*     xintp, to enable efficient parallelisation using "orphaned
*     directives", within a parallel context defined in this
*     subroutine, and to handle integrations on T-, p- and mixed
*     grids with only one routine, and suitable boundary weights.

*     Modules
      USE parameters, ONLY : nxpa, nypa, nxta, nyta, nla, nxaooc,
     &                       nyaooc, nxpo, nypo, nxto, nyto,
     &                       nx1, ny1, atnorm, ocnorm, fnot
      USE atconst, ONLY : gpat, rhoat, rdxaf0, dxam2,
     &                    ah4at, hat, cpat, dta, davgat
      USE occonst, ONLY : gpoc, rhooc, rdxof0, dxom2, ah2oc,
     &                    ah4oc, hoc, delek, cpoc, dto
#ifndef ocean_only
      USE atstate, ONLY : pa, pam, qa, wekta, wekpa, entat
#endif
#ifndef atmos_only
      USE ocstate, ONLY : po, pom, qo, wekto, wekpo, entoc
#endif
      USE intrfac, ONLY : sst, ast, tauxo, tauyo,
     &                    tauxa, tauya, hmixa, hmat
      USE radiate, ONLY : Aup, Bup, Cup, Dup
      USE monitor

      IMPLICIT NONE

*     Subroutine arguments
*
*     Local parameters
*
*     Local variables
      integer i,j,k
      double precision etaint,et2now,et2dot,pkeint,tauxav,tauyav,
     &                 ugeos,ugdot,utaux,uke,ukedot,u2diss,u4diss,ujet,
     &                 vgeos,vgdot,vtauy,vke,vkedot,v2diss,v4diss,
     &                 pint,qint,rgpoc,rgpat

#ifndef ocean_only
*     Define atmosphere variables
      double precision etaat(nxpa,nypa),etaatdot(nxpa,nypa),
     &                 ugat(nxpa,nyta),vgat(nxta,nypa),
     &                 attwk1(nxta,nypa),atpwk1(nxpa,nypa),
     &                 attwk2(nxta,nypa),atpwk2(nxpa,nypa),
     &                 attwk3(nxta,nypa),atpwk3(nxpa,nypa),
     &                 ujeta(nyta),paref(nla)
#endif /* not ocean_only */
#ifndef atmos_only
*     Define ocean variables
      double precision etaoc(nxpo,nypo),etaocdot(nxpo,nypo),
     &                 ugoc(nxpo,nyto),vgoc(nxto,nypo),
     &                 octwk1(nxto,nypo),ocpwk1(nxpo,nypo),
     &                 octwk2(nxto,nypo),ocpwk2(nxpo,nypo),
     &                 octwk3(nxto,nypo),ocpwk3(nxpo,nypo),
     &                 octwk4(nxto,nypo),ocpwk4(nxpo,nypo),
     &                 ujeto(nyto),pomin,pomax,poref(nlo),psiext
#endif /* not atmos_only */

*     Define reference pressures on equatorward side of domain
*     --------------------------------------------------------
#ifndef ocean_only
      if ( fnot.gt.0.0d0 ) then
        do k=1,nla
          paref(k) = pa(1,1,k)
        enddo
       else if ( fnot.lt.0.0d0 ) then
        do k=1,nla
          paref(k) = pa(1,nypa,k)
        enddo
      endif
#endif /* not ocean_only */
#ifndef atmos_only
      if ( fnot.gt.0.0d0 ) then
        do k=1,nlo
          poref(k) = po(1,1,k)
        enddo
       else if ( fnot.lt.0.0d0 ) then
        do k=1,nlo
          poref(k) = po(1,nypo,k)
        enddo
      endif
#endif /* not atmos_only */

#ifndef ocean_only
*     Atmosphere diagnostics
*     ======================

!$OMP PARALLEL DEFAULT (NONE)
!$OMP&         PRIVATE (i,j,k,rgpoc,rgpat)
!$OMP&         SHARED  (etaint,et2now,et2dot,pkeint,pint,qint,utaux,
!$OMP&                  uke,ukedot,u4diss,vtauy,vke,vkedot,v4diss)
!$OMP&         SHARED  (wekta,wetmat,attwk1,watmat,wekpa,wepmat,atpwk1,
!$OMP&                  wapmat,entat,entmat,enamat,etaat,etaatdot,
!$OMP&                  atpwk2,etamat,pkenat,ddtkeat,ddtpeat,gpat,rhoat,
!$OMP&                  et2mat,rdxaf0,pa,tauxa,tauya,attwk2,utauat,pam,
!$OMP&                  dxam2,ugat,vgat,atpwk3,ujeta,atstpos,atstval,
!$OMP&                  attwk3,qa,pavgat,qavgat,ah4dat,ah4at,hat,ast,
!$OMP&                  kealat,astmin,astmax,tmlmat,hmixa,hmlmat,
!$OMP&                  hcmlat,cpat,dta)

*     Ekman velocity diagnostics
*     --------------------------
*     Mean value of Wekman at T points
      call genint (wekta, nxta, nyta, 1.0d0, 1.0d0, wetmat)
*     Mean value of abs( Wekman ) at T points
!$OMP DO SCHEDULE (STATIC)
      do j=1,nyta
        do i=1,nxta
          attwk1(i,j) = abs( wekta(i,j) )
        enddo
      enddo
!$OMP END DO
      call genint (attwk1, nxta, nyta, 1.0d0, 1.0d0, watmat)
*     Mean value of Wekman at p points
      call genint (wekpa, nxpa, nypa, 0.5d0, 0.5d0, wepmat)
*     Mean value of abs( Wekman ) at p points
!$OMP DO SCHEDULE (STATIC)
      do j=1,nypa
        do i=1,nxpa
          atpwk1(i,j) = abs( wekpa(i,j) )
        enddo
      enddo
!$OMP END DO
      call genint (atpwk1, nxpa, nypa, 0.5d0, 0.5d0, wapmat)
!$OMP SINGLE
      wetmat = wetmat*atnorm
      watmat = watmat*atnorm
      wepmat = wepmat*atnorm
      wapmat = wapmat*atnorm
!$OMP END SINGLE NOWAIT

*     Entrainment diagnostics
*     -----------------------
*     We are assuming all entrainment is across interface 1.
      call genint (entat(1,1), nxpa, nypa, 0.5d0, 0.5d0, entmat(1))
!$OMP DO SCHEDULE (STATIC)
      do j=1,nypa
        do i=1,nxpa
          atpwk1(i,j) = abs( entat(i,j) )
        enddo
      enddo
!$OMP END DO
      call genint (atpwk1, nxpa, nypa, 0.5d0, 0.5d0, enamat(1))
!$OMP SINGLE
      entmat(1) = entmat(1)*atnorm
      enamat(1) = enamat(1)*atnorm
      do k=2,nla-1
        entmat(k) = 0.0d0
        enamat(k) = 0.0d0
      enddo
!$OMP END SINGLE NOWAIT

*     Interface displacement eta diagnostics
*     --------------------------------------
*     Loop over interfaces
      do k=1,nla-1
        rgpat = 1.0d0/gpat(k)
!$OMP   DO SCHEDULE (STATIC)
*       Infer eta and its instantaneous time derivative at p points
        do j=1,nypa
          do i=1,nxpa
            etaat(i,j) = rgpat*( pa(i,j,k) - pa(i,j,k+1) )
            etaatdot(i,j) = (rgpat/dta)*(  pa (i,j,k) - pa (i,j,k+1)
     &                                   - pam(i,j,k) + pam(i,j,k+1) )
            atpwk1(i,j) = etaat(i,j)*etaat(i,j)
            atpwk2(i,j) = etaat(i,j)*etaatdot(i,j)
            atpwk3(i,j) = etaat(i,j)*entat(i,j)
          enddo
        enddo
!$OMP   END DO
*       Compute integral of eta
        call genint (etaat, nxpa, nypa, 0.5d0, 0.5d0, etaint)
*       Compute integral of eta^2
        call genint (atpwk1, nxpa, nypa, 0.5d0, 0.5d0, et2now)
*       Compute integral of eta*d(eta)/dt
        call genint (atpwk2, nxpa, nypa, 0.5d0, 0.5d0, et2dot)
!$OMP   SINGLE
        etamat(k) = etaint*atnorm
        et2now = et2now*atnorm
        pkenat(k) = 0.0d0
        ddtpeat(k) = rhoat*gpat(k)*et2dot
        et2mat(k) = et2now
!$OMP   END SINGLE NOWAIT
*       Compute integral of eta*entrainment "e"
*       We are assuming all entrainment is across interface 1.
        if ( k.eq.1 ) then
          call genint (atpwk3, nxpa, nypa, 0.5d0, 0.5d0, pkeint)
!$OMP     SINGLE
          pkenat(1) = rhoat*gpat(1)*pkeint*atnorm
!$OMP     END SINGLE NOWAIT
        endif
      enddo

*     Atmosphere energetics
*     =====================

*     KE exchange between ocean & atmosphere
*     --------------------------------------
*     Infer u1 geostrophically
*     Compute required integrand
!$OMP DO SCHEDULE (STATIC)
!$OMP&   PRIVATE (ugeos,tauxav)
      do j=1,nyta
        do i=1,nxpa
          ugeos = -rdxaf0*( pa(i,j+1,1) - pa(i,j,1) )
          tauxav = 0.5d0*( tauxa(i,j+1) + tauxa(i,j) )
          atpwk1(i,j) = ugeos*tauxav
        enddo
      enddo
!$OMP END DO
      call genint (atpwk1, nxpa, nyta, 0.5d0, 1.0d0, utaux)
*     Infer v1 geostrophically
*     Compute required integrand
!$OMP DO SCHEDULE (STATIC)
!$OMP&   PRIVATE (vgeos,tauyav)
      do j=1,nypa
        do i=1,nxta
          vgeos = rdxaf0*( pa(i+1,j,1) - pa(i,j,1) )
          tauyav = 0.5d0*( tauya(i+1,j) + tauya(i,j) )
          attwk1(i,j) = vgeos*tauyav
        enddo
      enddo
!$OMP END DO
      call genint (attwk1, nxta, nypa, 1.0d0, 0.5d0, vtauy)
!$OMP SINGLE
      utauat = rhoat*( vtauy + utaux )*atnorm
!$OMP END SINGLE NOWAIT

*     Layer kinetic energy and its time derivative
*     --------------------------------------------
      do k=1,nla
*       Infer lagged u and v geostrophically; derive Del-4th
*       ----------------------------------------------------
!$OMP   DO SCHEDULE (STATIC)
        do j=1,nyta
          do i=1,nxpa
            ugat(i,j) = -rdxaf0*( pam(i,j+1,k) - pam(i,j,k) )
          enddo
        enddo
!$OMP   END DO
*       Differentiate contents of ugat to get Del-4th(lagged u)
*       Return this field in array atpwk1. atpwk3 is workspace
*       (which contains Del-sqd(lagged u))
        call del4ch (ugat, nxpa, nyta, atpwk3, dxam2, atpwk1)
*       Infer lagged v geostrophically
!$OMP   DO SCHEDULE (STATIC)
        do j=1,nypa
          do i=1,nxta
            vgat(i,j) = rdxaf0*( pam(i+1,j,k) - pam(i,j,k) )
          enddo
        enddo
!$OMP   END DO
*       Differentiate contents of vgat to get Del-4th(lagged v).
*       Return this field in array attwk2. attwk3 is workspace
*       (which contains Del-sqd(lagged v))
        call del4ch (vgat, nxta, nypa, attwk3, dxam2, attwk2)
*       Infer current u and v geostrophically
*       -------------------------------------
*       Compute required integrands and zonal jet position
!$OMP   DO SCHEDULE (STATIC)
!$OMP&     PRIVATE (ujet,ugeos,ugdot)
        do j=1,nyta
          ujet = 0.0d0
          do i=1,nxpa
            ugeos = -rdxaf0*( pa(i,j+1,k) - pa(i,j,k) )
            ugdot = -(rdxaf0/dta)*(  pa (i,j+1,k) - pa (i,j,k)
     &                             - pam(i,j+1,k) + pam(i,j,k) )
            ujet = ujet + ugeos
            atpwk1(i,j) = ugeos*atpwk1(i,j)
            atpwk2(i,j) = ugeos*ugeos
            atpwk3(i,j) = ugeos*ugdot
          enddo
*         The above loop double counts the end values, which
*         are equal, either by cyclicity or because ugeos = 0
          ujet = ujet - ugeos
          ujeta(j) = abs( ujet )/dble(nxta)
        enddo
!$OMP   END DO
*       Find position and value of avged max in each layer
!$OMP   SINGLE
        atstpos(k) = 0
        atstval(k) = 0.0d0
        do j=1,nyta
          if ( ujeta(j).gt.atstval(k) ) then
            atstpos(k) = j
            atstval(k) = ujeta(j)
          endif
        enddo
!$OMP   END SINGLE NOWAIT
        call genint (atpwk1, nxpa, nyta, 0.5d0, 1.0d0, u4diss)
        call genint (atpwk2, nxpa, nyta, 0.5d0, 1.0d0, uke)
        call genint (atpwk3, nxpa, nyta, 0.5d0, 1.0d0, ukedot)
*       Compute required integrands
!$OMP   DO SCHEDULE (STATIC)
!$OMP&     PRIVATE (vgeos,vgdot)
        do j=1,nypa
          do i=1,nxta
            vgeos = rdxaf0*( pa(i+1,j,k) - pa(i,j,k) )
            vgdot = (rdxaf0/dta)*(  pa (i+1,j,k) - pa (i,j,k)
     &                            - pam(i+1,j,k) + pam(i,j,k) )
            attwk1(i,j) = vgeos*attwk2(i,j)
            attwk2(i,j) = vgeos*vgeos
          enddo
        enddo
!$OMP   END DO
        call genint (attwk1, nxta, nypa, 1.0d0, 0.5d0, v4diss)
        call genint (attwk2, nxta, nypa, 1.0d0, 0.5d0, vke)
        call genint (attwk3, nxta, nypa, 1.0d0, 0.5d0, vkedot)
        call genint (pa(1,1,k), nxpa, nypa, 0.5d0, 0.5d0, pint)
        call genint (qa(1,1,k), nxpa, nypa, 0.5d0, 0.5d0, qint)
!$OMP   SINGLE
        pavgat(k) = pint*atnorm
        qavgat(k) = qint*atnorm
        ah4dat(k) = rhoat*ah4at(k)*hat(k)*( u4diss + v4diss )*atnorm
        kealat(k) = 0.5d0*rhoat*hat(k)*( uke + vke )*atnorm
        ddtkeat(k) = rhoat*hat(k)*( ukedot + vkedot )*atnorm
!$OMP   END SINGLE NOWAIT
      enddo

*     Mixed layer temperature & thickness diagnostics
*     ===============================================

*     Compute mean atmos. mixed layer temperature and thickness
*     ---------------------------------------------------------
      call genint (ast, nxta, nyta, 1.0d0, 1.0d0, tmlmat)
      call genint (hmixa, nxta, nyta, 1.0d0, 1.0d0, hmlmat)
*     Compute total heat content of atmos m.l.
*     Find extrema of rel. atmos. m.l. temperature
*     --------------------------------------------
!$OMP SINGLE
      astmin =  1.0d30
      astmax = -1.0d30
!$OMP END SINGLE
!$OMP DO SCHEDULE (STATIC)
!$OMP&   REDUCTION(min:astmin) REDUCTION(max:astmax)
      do j=1,nyta
        do i=1,nxta
          astmin = min( astmin, ast(i,j) )
          astmax = max( astmax, ast(i,j) )
          attwk1(i,j) = ast(i,j)*hmixa(i,j)
        enddo
      enddo
!$OMP END DO
      call genint (attwk1, nxta, nyta, 1.0d0, 1.0d0, hcmlat)
!$OMP SINGLE
      tmlmat = tmlmat*atnorm
      hmlmat = hmlmat*atnorm
      hcmlat = rhoat*cpat*hcmlat*atnorm
!$OMP END SINGLE NOWAIT

!$OMP END PARALLEL

*     Compute mean atmos. mixed layer temperature over ocean
*     ------------------------------------------------------
      tmaooc = 0.0d0
      do j=ny1,ny1+nyaooc-1
        do i=nx1,nx1+nxaooc-1
          tmaooc = tmaooc + ast(i,j)
        enddo
      enddo
      tmaooc = tmaooc/dble( nxaooc*nyaooc )

*     Compute mean outgoing long wave radiation
*     -----------------------------------------
      olrtop = Bup(nla)*(hmlmat-hmat)
     &        + Cup(nla)*davgat + Dup(nla)*tmlmat
      do i=1,nla-1
        olrtop = olrtop + Aup(nla,i)*etamat(i)
      enddo

*     Velocity component extrema and max. Courant
*     number in atmospheric layers (mixed and Q-G)
*     ============================================
      call courat

#endif /* not ocean_only */

#ifndef atmos_only
*     Ocean diagnostics
*     =================

!$OMP PARALLEL DEFAULT (NONE)
!$OMP&         PRIVATE (i,j,k,rgpoc,rgpat)
!$OMP&         SHARED  (etaint,et2now,et2dot,pkeint,pint,qint,
!$OMP&                  utaux,uke,ukedot,u2diss,u4diss,
!$OMP&                  vtauy,vke,vkedot,v2diss,v4diss)
!$OMP&         SHARED  (wekto,wetmoc,octwk1,watmoc,wekpo,wepmoc,ocpwk1,
!$OMP&                  wapmoc,entoc,entmoc,enamoc,gpoc,etaoc,etaocdot,
!$OMP&                  po,ocpwk2,pkenoc,rhooc,etamoc,ddtkeoc,ddtpeoc,
!$OMP&                  et2moc,rdxof0,tauxo,tauyo,utauoc,pom,dxom2,ugoc,
!$OMP&                  vgoc,octwk2,octwk3,pomin,pomax,ocpwk3,ujeto,
!$OMP&                  ocjpos,ocjval,qo,pavgoc,qavgoc,ah2doc,ah2oc,
!$OMP&                  ah4doc,ah4oc,hoc,kealoc,psiext,osfmin,osfmax,
!$OMP&                  poref,occirc,btdgoc,delek,sst,sstmin,sstmax,
!$OMP&                  tmlmoc,hfmloc,octwk4,ocpwk4,cpoc,dto)

*     Ekman velocity diagnostics
*     --------------------------
*     Mean value of Wekman at T points
      call genint (wekto, nxto, nyto, 1.0d0, 1.0d0, wetmoc)
*     Mean value of abs( Wekman ) at T points
!$OMP DO SCHEDULE (STATIC)
      do j=1,nyto
        do i=1,nxto
          octwk1(i,j) = abs( wekto(i,j) )
        enddo
      enddo
!$OMP END DO
      call genint (octwk1, nxto, nyto, 1.0d0, 1.0d0, watmoc)
*     Mean value of Wekman at p points
      call genint (wekpo, nxpo, nypo, 0.5d0, 0.5d0, wepmoc)
*     Mean value of abs( Wekman ) at p points
!$OMP DO SCHEDULE (STATIC)
      do j=1,nypo
        do i=1,nxpo
          ocpwk1(i,j) = abs( wekpo(i,j) )
        enddo
      enddo
!$OMP END DO
      call genint (ocpwk1, nxpo, nypo, 0.5d0, 0.5d0, wapmoc)
!$OMP SINGLE
      wetmoc = wetmoc*ocnorm
      watmoc = watmoc*ocnorm
      wepmoc = wepmoc*ocnorm
      wapmoc = wapmoc*ocnorm
!$OMP END SINGLE NOWAIT

*     Entrainment diagnostics
*     -----------------------
      call genint (entoc, nxpo, nypo, 0.5d0, 0.5d0, entmoc)
!$OMP DO SCHEDULE (STATIC)
      do j=1,nypo
        do i=1,nxpo
          ocpwk1(i,j) = abs( entoc(i,j) )
        enddo
      enddo
!$OMP END DO
      call genint (ocpwk1, nxpo, nypo, 0.5d0, 0.5d0, enamoc)
!$OMP SINGLE
      entmoc = entmoc*ocnorm
      enamoc = enamoc*ocnorm
!$OMP END SINGLE NOWAIT

*     Interface displacement eta diagnostics
*     --------------------------------------
*     Loop over interfaces in ocean
      do k=1,nlo-1
        rgpoc = 1.0d0/gpoc(k)
!$OMP   DO SCHEDULE (STATIC)
*       Infer eta and its instantaneous time derivative at p points
        do j=1,nypo
          do i=1,nxpo
            etaoc(i,j) = rgpoc*( po(i,j,k+1) - po(i,j,k) )
            etaocdot(i,j) = (rgpoc/dto)*(  po (i,j,k) - po (i,j,k+1)
     &                                   - pom(i,j,k) + pom(i,j,k+1) )
            ocpwk1(i,j) = etaoc(i,j)*etaoc(i,j)
            ocpwk2(i,j) = etaoc(i,j)*etaocdot(i,j)
            ocpwk3(i,j) = etaoc(i,j)*entoc(i,j)
          enddo
        enddo
!$OMP   END DO
*       Compute integral of eta
        call genint (etaoc, nxpo, nypo, 0.5d0, 0.5d0, etaint)
*       Compute integral of eta^2
        call genint (ocpwk1, nxpo, nypo, 0.5d0, 0.5d0, et2now)
*       Compute integral of eta*d(eta)/dt
        call genint (ocpwk2, nxpo, nypo, 0.5d0, 0.5d0, et2dot)
!$OMP   SINGLE
        etamoc(k) = etaint*ocnorm
        et2now = et2now*ocnorm
        ddtpeoc(k) = rhooc*gpoc(k)*et2dot
        et2moc(k) = et2now
!$OMP   END SINGLE NOWAIT
*       Compute integral of eta*entrainment "e"
*       Only have entrainment across interface 1 in ocean
        if ( k.eq.1 ) then
          call genint (ocpwk3, nxpo, nypo, 0.5d0, 0.5d0, pkeint)
!$OMP     SINGLE
          pkenoc = rhooc*gpoc(1)*pkeint*ocnorm
!$OMP     END SINGLE NOWAIT
        endif
      enddo

*     Ocean energetics
*     ================

*     KE exchange between ocean & atmosphere
*     --------------------------------------
*     Infer u1 geostrophically
*     Compute required integrand
!$OMP DO SCHEDULE (STATIC)
!$OMP&   PRIVATE (ugeos,tauxav)
      do j=1,nyto
        do i=1,nxpo
          ugeos = -rdxof0*( po(i,j+1,1) - po(i,j,1) )
          tauxav = 0.5d0*( tauxo(i,j+1) + tauxo(i,j) )
          ocpwk1(i,j) = ugeos*tauxav
        enddo
      enddo
!$OMP END DO
      call genint (ocpwk1, nxpo, nyto, 0.5d0, 1.0d0, utaux)
*     Infer v1 geostrophically
*     Compute required integrand
!$OMP DO SCHEDULE (STATIC)
!$OMP&   PRIVATE (vgeos,tauyav)
      do j=1,nypo
        do i=1,nxto
          vgeos = rdxof0*( po(i+1,j,1) - po(i,j,1) )
          tauyav = 0.5d0*( tauyo(i+1,j) + tauyo(i,j) )
          octwk1(i,j) = vgeos*tauyav
        enddo
      enddo
!$OMP END DO
      call genint (octwk1, nxto, nypo, 1.0d0, 0.5d0, vtauy)
!$OMP SINGLE
      utauoc = rhooc*( vtauy + utaux )*ocnorm
!$OMP END SINGLE NOWAIT

*     Layer kinetic energy and its time derivative, also stream function
*     ------------------------------------------------------------------
      do k=1,nlo
*       Infer lagged u and v geostrophically; derive Del-4th
*       ----------------------------------------------------
!$OMP   DO SCHEDULE (STATIC)
        do j=1,nyto
          do i=1,nxpo
            ugoc(i,j) = -rdxof0*( pom(i,j+1,k) - pom(i,j,k) )
          enddo
        enddo
!$OMP   END DO
*       Differentiate contents of ugoc to get Del-4th(lagged u)
*       Return this field in array ocpwk2. ocpwk1 is workspace
*       which contains Del-sqd(lagged u)
#  ifdef cyclic_ocean
        call del4ch (ugoc, nxpo, nyto, ocpwk1, dxom2, ocpwk2)
#  else
        call del4bx (ugoc, nxpo, nyto, ocpwk1, dxom2, ocpwk2)
#  endif
!$OMP   DO SCHEDULE (STATIC)
        do j=1,nypo
          do i=1,nxto
            vgoc(i,j) = rdxof0*( pom(i+1,j,k) - pom(i,j,k) )
          enddo
        enddo
!$OMP   END DO
*       Differentiate contents of vgoc to get Del-4th(lagged v)
*       Return this field in array octwk2. octwk1 is workspace
*       which contains Del-sqd(lagged v)
#  ifdef cyclic_ocean
        call del4ch (vgoc, nxto, nypo, octwk1, dxom2, octwk2)
#  else
        call del4bx (vgoc, nxto, nypo, octwk1, dxom2, octwk2)
#  endif
*       Find extrema of po for stream function
*       --------------------------------------
!$OMP   SINGLE
        pomin =  1.0d30
        pomax = -1.0d30
!$OMP   END SINGLE
!$OMP   DO SCHEDULE (STATIC)
!$OMP&     REDUCTION(min:pomin) REDUCTION(max:pomax)
        do j=1,nypo
          do i=1,nxpo
            pomin = min( pomin, po(i,j,k) )
            pomax = max( pomax, po(i,j,k) )
          enddo
        enddo
!$OMP   END DO
*       Infer current u and v geostrophically
*       -------------------------------------
*       Compute required integrands and zonal jet position
!$OMP   DO SCHEDULE (STATIC)
!$OMP&     PRIVATE (ujet,ugeos,ugdot)
        do j=1,nyto
          ujet = 0.0d0
          do i=1,nxpo
            ugeos = -rdxof0*( po(i,j+1,k) - po(i,j,k) )
            ugdot = -(rdxof0/dto)*(  po (i,j+1,k) - pom(i,j,k)
     &                             - pom(i,j+1,k) + pom(i,j,k) )
            ujet = ujet + ugeos
*           Dissipative terms
            ocpwk1(i,j) = ugeos*ocpwk1(i,j)
            ocpwk2(i,j) = ugeos*ocpwk2(i,j)
*           Kinetic energy and its time derivative
            ocpwk3(i,j) = ugeos*ugeos
            ocpwk4(i,j) = ugeos*ugdot
          enddo
*         The above loop double counts the end values, which
*         are equal, either by cyclicity or because ugeos = 0
          ujet = ujet - ugeos
          ujeto(j) = abs( ujet )/dble(nxto)
        enddo
!$OMP   END DO
*       Find position and value of avged max in each layer
!$OMP   SINGLE
        ocjpos(k) = 0
        ocjval(k) = 0.0d0
        do j=1,nyto
          if ( ujeto(j).gt.ocjval(k) ) then
            ocjpos(k) = j
            ocjval(k) = ujeto(j)
          endif
        enddo
!$OMP   END SINGLE NOWAIT
*       Dissipative terms
        call genint (ocpwk1, nxpo, nyto, 0.5d0, 1.0d0, u2diss)
        call genint (ocpwk2, nxpo, nyto, 0.5d0, 1.0d0, u4diss)
*       Kinetic energy and its time derivative
        call genint (ocpwk3, nxpo, nyto, 0.5d0, 1.0d0, uke)
        call genint (ocpwk4, nxpo, nyto, 0.5d0, 1.0d0, ukedot)
*       Compute required integrands
!$OMP   DO SCHEDULE (STATIC)
!$OMP&     PRIVATE (vgeos,vgdot)
        do j=1,nypo
          do i=1,nxto
            vgeos = rdxof0*( po(i+1,j,k) - po(i,j,k) )
            vgdot = (rdxof0/dto)*(  po (i+1,j,k) - po (i,j,k)
     &                            - pom(i+1,j,k) + pom(i,j,k) )
*           Dissipative terms
            octwk1(i,j) = vgeos*octwk1(i,j)
            octwk2(i,j) = vgeos*octwk2(i,j)
*           Kinetic energy and its time derivative
            octwk3(i,j) = vgeos*vgeos
            octwk4(i,j) = vgeos*vgdot
          enddo
        enddo
!$OMP   END DO
        call genint (po(1,1,k), nxpo, nypo, 0.5d0, 0.5d0, pint)
        call genint (qo(1,1,k), nxpo, nypo, 0.5d0, 0.5d0, qint)
*       Dissipative terms
        call genint (octwk1, nxto, nypo, 1.0d0, 0.5d0, v2diss)
        call genint (octwk2, nxto, nypo, 1.0d0, 0.5d0, v4diss)
*       Kinetic energy and its time derivative
        call genint (octwk3, nxto, nypo, 1.0d0, 0.5d0, vke)
        call genint (octwk4, nxto, nypo, 1.0d0, 0.5d0, vkedot)
!$OMP   SINGLE
        pavgoc(k) = pint*ocnorm
        qavgoc(k) = qint*ocnorm
        ah2doc(k) = -rhooc*ah2oc(k)*hoc(k)*( u2diss + v2diss )*ocnorm
        ah4doc(k) =  rhooc*ah4oc(k)*hoc(k)*( u4diss + v4diss )*ocnorm
        kealoc(k) = 0.5d0*rhooc*hoc(k)*( uke + vke )*ocnorm
        ddtkeoc(k) = rhooc*hoc(k)*( ukedot + vkedot )*ocnorm
*       Convert dynamic pressure extrema to volume
*       transport stream function values (in Sverdrups)
        psiext = min( pomin/fnot, pomax/fnot )
        osfmin(k) = 1.0d-6*hoc(k)*( psiext - poref(k)/fnot )
        psiext = max( pomin/fnot, pomax/fnot )
        osfmax(k) = 1.0d-6*hoc(k)*( psiext - poref(k)/fnot )
*       Compute volume transport in each layer (in Sverdrups)
        occirc(k) = 1.0d-6*hoc(k)*( po(1,1,k) - po(1,nypo,k) )/fnot
!$OMP   END SINGLE NOWAIT
      enddo

*     Drag in Ekman layer at ocean bottom
*     -----------------------------------
*     Infer lagged u geostrophically
*     Compute required integrand
!$OMP DO SCHEDULE (STATIC)
!$OMP&   PRIVATE (ugeos)
      do j=1,nyto
        do i=1,nxpo
          ugeos = -rdxof0*( pom(i,j+1,nlo) - pom(i,j,nlo) )
          ocpwk1(i,j) = ugeos*ugeos
        enddo
      enddo
!$OMP END DO
      call genint (ocpwk1, nxpo, nyto, 0.5d0, 1.0d0, u2diss)
*     Infer lagged v geostrophically
*     Compute required integrand
!$OMP DO SCHEDULE (STATIC)
!$OMP&   PRIVATE (vgeos)
      do j=1,nypo
        do i=1,nxto
          vgeos = rdxof0*( pom(i+1,j,nlo) - pom(i,j,nlo) )
          octwk1(i,j) = vgeos*vgeos
        enddo
      enddo
!$OMP END DO
      call genint (octwk1, nxto, nypo, 1.0d0, 0.5d0, v2diss)
!$OMP SINGLE
      btdgoc = 0.5d0*rhooc*delek*abs(fnot)*( u2diss + v2diss )*ocnorm
!$OMP END SINGLE NOWAIT

*     Mixed layer temperature & thickness diagnostics
*     ===============================================

*     Find extrema of rel. ocean m.l. temperature
*     Compute flux of heat at bottom of ocean m.l.
*     --------------------------------------------
!$OMP SINGLE
      sstmin =  1.0d30
      sstmax = -1.0d30
!$OMP END SINGLE
!$OMP DO SCHEDULE (STATIC)
!$OMP&   REDUCTION(min:sstmin) REDUCTION(max:sstmax)
      do j=1,nyto
        do i=1,nxto
          sstmin = min( sstmin, sst(i,j) )
          sstmax = max( sstmax, sst(i,j) )
          octwk1(i,j) = sst(i,j)*wekto(i,j)
        enddo
      enddo
!$OMP END DO
      call genint (octwk1, nxto, nyto, 1.0d0, 1.0d0, hfmloc)
*     Compute mean ocean mixed layer temperature
*     ------------------------------------------
      call genint (sst, nxto, nyto, 1.0d0, 1.0d0, tmlmoc)
!$OMP SINGLE
      tmlmoc = tmlmoc*ocnorm
      hfmloc = rhooc*cpoc*hfmloc*ocnorm
!$OMP END SINGLE NOWAIT

!$OMP END PARALLEL

*     Compute total ocean circulation
*     -------------------------------
      occtot = 0.0d0
      do k=1,nlo
        occtot = occtot + occirc(k)
      enddo

*     Velocity component extrema and max. Courant
*     number in oceanic layers (mixed and Q-G)
*     ===========================================
      call couroc

#endif /* not atmos_only */

*     Check diagnostic values
*     =======================
**    print *,' '

*     Atmosphere diagnostics
*     ----------------------
#ifndef ocean_only
**    write(*,215) '  Atmos. entrainment   entmat = ',
**   &             (entmat(k),k=1,nla-1)
**    write(*,215) '  |Atmos. entrainment| enamat = ',
**   &             (enamat(k),k=1,nla-1)
**    write(*,215) '  Mean atmos. interfac etamat = ',
**   &             (etamat(k),k=1,nla-1)
**    write(*,215) '  Mean atmos. eta^2    et2mat = ',
**   &             (et2mat(k),k=1,nla-1)
**    write(*,215) '  Mean energy transfer utauat = ',utauat
**    write(*,201) '  Atmos. s.t. positn. atstpos = ',
**   &             (atstpos(k),k=1,nla)
**    write(*,215) '  Atmos. storm tracks atstval = ',
**   &             (atstval(k),k=1,nla)
**    write(*,215) '  Mean atmos. Del-4th  ah4dat = ',
**   &             (ah4dat(k),k=1,nla)
**    write(*,215) '  Mean atmos. layer KE kealat = ',
**   &             (kealat(k),k=1,nla)
**    write(*,215) '  Atmos. m.l. T, H, heat cap. = ',
**   &             tmlmat,hmlmat,hcmlat
#endif /* not ocean_only */

*     Ocean diagnostics
*     -----------------
#ifndef atmos_only
**    write(*,215) '  wetmoc,watmoc,wepmoc,wapmoc = ',
**   &             wetmoc,watmoc,wepmoc,wapmoc
**    write(*,215) '  Ocean entr.  entmoc, enamoc = ',entmoc,enamoc
**    write(*,215) '  Mean ocean  interfac etamoc = ',
**   &             (etamoc(k),k=1,nlo-1)
**    write(*,215) '  Mean ocean  eta^2    et2moc = ',
**   &             (et2moc(k),k=1,nlo-1)
**    write(*,215) '  Mean energy transfer utauoc = ',utauoc
**    write(*,201) '  Ocean  jet positions ocjpos = ',
**   &             (ocjpos(k),k=1,nlo)
**    write(*,215) '  Ocean  jet strengths ocjval = ',
**   &             (ocjval(k),k=1,nlo)
**    write(*,215) '  Mean ocean  Del-sqd  ah2doc = ',
**                 (ah2doc(k),k=1,nlo)
**    write(*,215) '  Mean ocean  Del-4th  ah4doc = ',
**   &             (ah4doc(k),k=1,nlo)
**    write(*,215) '  Mean ocean  layer KE kealoc = ',
**   &             (kealoc(k),k=1,nlo)
**    write(*,215) '  Ocean stream func. min (Sv) = ',
**   &             (osfmin(k),k=1,nlo)
**    write(*,215) '  Ocean stream func. max (Sv) = ',
**   &             (osfmax(k),k=1,nlo)
**    write(*,215) '  Ocean layer transport  (Sv) = ',
**   &             (occirc(k),k=1,nlo)
**    write(*,215) '  Bottom drag dissip.  btdgoc = ',btdgoc
**    write(*,215) '  Ocean  m.l. temp, heat flux = ',
**   &             tmlmoc,hfmloc
#endif /* not atmos_only */

**201 format(a,5i13)
**215 format(a,1p,5d13.5)

      END SUBROUTINE monnc_comp
*
c***********************************************************************
*
#ifndef cyclic_ocean
*
      SUBROUTINE del4bx (arr, nx, ny, del2, dxm2, del4)
*
*     Computes Del-4th(arr) for a field arr(nx,ny) tabulated
*     in a finite box. No boundary condition is applied;
*     the boundary values are computed using the appropriate
*     forward or backward difference formulae.
*     del2(nx,ny) is a workspace array used for Del-sqd(arr);
*     Del-4th(arr) is returned in the array del4(nx,ny).
*     arr is returned unchanged.
*
*     N.B. Parallelised using "orphaned directives";
*     only intended for use in the parallel regions
*     established earlier in calling subroutine monnc_comp

      IMPLICIT NONE

*     Subroutine arguments
      integer nx,ny
      double precision arr(nx,ny),del2(nx,ny),del4(nx,ny),dxm2
*
*     Local variables
      integer i,j
*
*     Compute Del-sqd(arr) except at N & S boundaries.
!$OMP DO SCHEDULE (STATIC)
!$OMP&   PRIVATE (i,j)
      do j=2,ny-1
*       Western boundary:
        del2( 1,j) = dxm2*(  arr( 3,j)-2.0d0*arr(  2 ,j)+arr(  1 ,j)
     &                     + arr( 1,j-1)-2.0d0*arr( 1,j)+arr( 1,j+1) )
*       Inner points (standard case)
        do i=2,nx-1
          del2(i,j) = dxm2*(  arr(i,j-1) + arr(i-1,j) + arr(i+1,j)
     &                      + arr(i,j+1) - 4.0d0*arr(i,j) )
        enddo
*       Eastern boundary:
        del2(nx,j) = dxm2*(  arr(nx,j)-2.0d0*arr(nx-1,j)+arr(nx-2,j)
     &                     + arr(nx,j-1)-2.0d0*arr(nx,j)+arr(nx,j+1) )
      enddo
!$OMP END DO NOWAIT
*
*     Compute N & S Del-sqd boundary values
!$OMP SINGLE
      do i=2,nx-1
        del2(i, 1) = dxm2*(  arr(i-1, 1)-2.0d0*arr(i, 1)+arr(i+1,1)
     &                     + arr(i, 3)-2.0d0*arr(i,  2 )+arr(i,  1 ) )
        del2(i,ny) = dxm2*(  arr(i-1,ny)-2.0d0*arr(i,ny)+arr(i+1,ny)
     &                     + arr(i,ny)-2.0d0*arr(i,ny-1)+arr(i,ny-2) )
      enddo
*
*     SW corner
      del2( 1, 1) = dxm2*(  arr( 3, 1) - 2.0d0*arr( 2 , 1 )
     &                    + arr( 1 , 1 )
     &                    + arr( 1, 3) - 2.0d0*arr( 1 , 2 )
     &                    + arr( 1 , 1 ) )
*     SE corner
      del2(nx, 1) = dxm2*(  arr(nx, 1) - 2.0d0*arr(nx-1, 1)
     &                    + arr(nx-2, 1)
     &                    + arr(nx, 3) - 2.0d0*arr( nx , 2)
     &                    + arr( nx , 1) )
*     NW corner
      del2( 1,ny) = dxm2*(  arr( 3,ny) - 2.0d0*arr( 2, ny )
     &                    + arr( 1, ny )
     &                    + arr( 1,ny) - 2.0d0*arr( 1,ny-1)
     &                    + arr( 1,ny-2) )
*     NE corner
      del2(nx,ny) = dxm2*(  arr(nx,ny) - 2.0d0*arr(nx-1,ny)
     &                    + arr(nx-2,ny)
     &                    + arr(nx,ny) - 2.0d0*arr(nx,ny-1)
     &                    + arr(nx,ny-2) )
!$OMP END SINGLE
*
*     Compute Del-4th(arr) except at N & S boundaries.
!$OMP DO SCHEDULE (STATIC)
!$OMP&   PRIVATE (i,j)
      do j=2,ny-1
*       Western boundary:
        del4( 1,j) = dxm2*( del2( 3,j)-2.0d0*del2(  2 ,j)+del2(  1 ,j)
     &                     +del2( 1,j-1)-2.0d0*del2( 1,j)+del2( 1,j+1) )
*       Inner points (standard case)
        do i=2,nx-1
          del4(i,j) = dxm2*(  del2(i,j-1) + del2(i-1,j) + del2(i+1,j)
     &                      + del2(i,j+1) - 4.0d0*del2(i,j) )
        enddo
*       Eastern boundary:
        del4(nx,j) = dxm2*( del2(nx,j)-2.0d0*del2(nx-1,j)+del2(nx-2,j)
     &                     +del2(nx,j-1)-2.0d0*del2(nx,j)+del2(nx,j+1) )
      enddo
!$OMP END DO NOWAIT
*
*     Compute N & S Del-4th boundary values
!$OMP SINGLE
      do i=2,nx-1
        del4(i, 1) = dxm2*( del2(i-1, 1)-2.0d0*del2(i, 1)+del2(i+1,1)
     &                     +del2(i, 3)-2.0d0*del2(i,  2 )+del2(i,  1 ) )
        del4(i,ny) = dxm2*( del2(i-1,ny)-2.0d0*del2(i,ny)+del2(i+1,ny)
     &                     +del2(i,ny)-2.0d0*del2(i,ny-1)+del2(i,ny-2) )
      enddo
*
*     SW corner
      del4( 1, 1) = dxm2*(  del2( 3, 1) - 2.0d0*del2( 2 , 1 )
     &                    + del2( 1 , 1 )
     &                    + del2( 1, 3) - 2.0d0*del2( 1 , 2 )
     &                    + del2( 1 , 1 ) )
*     SE corner
      del4(nx, 1) = dxm2*(  del2(nx, 1) - 2.0d0*del2(nx-1, 1)
     &                    + del2(nx-2, 1)
     &                    + del2(nx, 3) - 2.0d0*del2( nx , 2)
     &                    + del2( nx , 1) )
*     NW corner
      del4( 1,ny) = dxm2*(  del2( 3,ny) - 2.0d0*del2( 2, ny )
     &                    + del2( 1, ny )
     &                    + del2( 1,ny) - 2.0d0*del2( 1,ny-1)
     &                    + del2( 1,ny-2) )
*     NE corner
      del4(nx,ny) = dxm2*(  del2(nx,ny) - 2.0d0*del2(nx-1,ny)
     &                    + del2(nx-2,ny)
     &                    + del2(nx,ny) - 2.0d0*del2(nx,ny-1)
     &                    + del2(nx,ny-2) )
!$OMP END SINGLE
*
      END SUBROUTINE del4bx
*
#endif /* not cyclic_ocean */
*
c***********************************************************************
*
      SUBROUTINE del4ch (arr, nx, ny, del2, dxm2, del4)
*
*     Computes Del-4th(arr) for a field arr(nx,ny) tabulated
*     in a periodic channel. No boundary condition is applied;
*     the boundary values are computed using the appropriate
*     forward or backward difference formulae.
*     del2(nx,ny) is a workspace array used for Del-sqd(arr);
*     we use an externally dimensioned array so that the
*     routine can be used for either atmosphere or ocean.
*     Del-4th(arr) is returned in the array del4(nx,ny).
*     arr is returned unchanged.
*
*     N.B. Parallelised using "orphaned directives";
*     only intended for use in the parallel regions
*     established earlier in calling subroutine monnc_comp

      IMPLICIT NONE

*     Subroutine arguments
      integer nx,ny
      double precision arr(nx,ny),del2(nx,ny),del4(nx,ny),dxm2
*
*     Local variables
      integer i,j
*
*     Compute Del-sqd(arr) except at N & S boundaries.
*     Apply periodic E-W condition
!$OMP DO SCHEDULE (STATIC)
!$OMP&   PRIVATE (i,j)
      do j=2,ny-1
*       Western boundary (periodic)
        del2( 1,j) = dxm2*(  arr(1,j-1) + arr(nx,j) + arr(2,j)
     &                     + arr(1,j+1) - 4.0d0*arr(1,j) )
*       Inner points (standard case)
        do i=2,nx-1
          del2(i,j) = dxm2*(  arr(i,j-1) + arr(i-1,j) + arr(i+1,j)
     &                      + arr(i,j+1) - 4.0d0*arr(i,j) )
        enddo
*       Eastern boundary (periodic)
        del2(nx,j) = dxm2*(  arr(nx,j-1) + arr(nx-1,j) + arr(1,j)
     &                     + arr(nx,j+1) - 4.0d0*arr(nx,j) )
      enddo
!$OMP END DO NOWAIT
*
*     Compute N & S Del-sqd boundary values
*     Apply periodic E-W condition
!$OMP SINGLE
      del2( 1, 1) = dxm2*(  arr( nx, 1) - 2.0d0*arr( 1,  1 )
     &                    + arr( 2 , 1)
     &                    + arr( 1 , 3) - 2.0d0*arr( 1,  2 )
     &                    + arr(1,  1 ) )
      del2( 1,ny) = dxm2*(  arr( nx,ny) - 2.0d0*arr( 1, ny )
     &                    + arr( 2 ,ny)
     &                    + arr( 1 ,ny) - 2.0d0*arr( 1,ny-1)
     &                    + arr(1,ny-2) )
      do i=2,nx-1
        del2(i, 1) = dxm2*(  arr(i-1, 1) - 2.0d0*arr( i,  1 )
     &                     + arr(i+1, 1)
     &                     + arr( i , 3) - 2.0d0*arr( i,  2 )
     &                     + arr(i,  1 ) )
        del2(i,ny) = dxm2*(  arr(i-1,ny) - 2.0d0*arr( i, ny )
     &                     + arr(i+1,ny)
     &                     + arr( i ,ny) - 2.0d0*arr( i,ny-1)
     &                     + arr(i,ny-2) )
      enddo
      del2(nx, 1) = dxm2*(  arr(nx-1, 1) - 2.0d0*arr(nx,  1 )
     &                    + arr( 1 , 1)
     &                    + arr( nx , 3) - 2.0d0*arr(nx,  2 )
     &                    + arr(nx,  1 ) )
      del2(nx,ny) = dxm2*(  arr(nx-1,ny) - 2.0d0*arr(nx, ny )
     &                    + arr( 1 ,ny)
     &                    + arr( nx ,ny) - 2.0d0*arr(nx,ny-1)
     &                    + arr(nx,ny-2) )
!$OMP END SINGLE
*
*     Compute Del-4th(arr) except at N & S boundaries.
*     Apply periodic E-W condition
!$OMP DO SCHEDULE (STATIC)
!$OMP&   PRIVATE (i,j)
      do j=2,ny-1
*       Western boundary (periodic)
        del4( 1,j) = dxm2*(  del2(1,j-1) + del2(nx,j) + del2(2,j)
     &                     + del2(1,j+1) - 4.0d0*del2(1,j) )
*       Inner points (standard case)
        do i=2,nx-1
          del4(i,j) = dxm2*(  del2(i,j-1) + del2(i-1,j) + del2(i+1,j)
     &                      + del2(i,j+1) - 4.0d0*del2(i,j) )
        enddo
*       Eastern boundary (periodic)
        del4(nx,j) = dxm2*(  del2(nx,j-1) + del2(nx-1,j) + del2(1,j)
     &                     + del2(nx,j+1) - 4.0d0*del2(nx,j) )
      enddo
!$OMP END DO NOWAIT
*
*     Compute N & S Del-4th boundary values
*     Apply periodic E-W condition
!$OMP SINGLE
      del4( 1, 1) = dxm2*(  del2( nx, 1) - 2.0d0*del2( 1,  1 )
     &                    + del2( 2 , 1)
     &                    + del2( 1 , 3) - 2.0d0*del2( 1,  2 )
     &                    + del2(1,  1 ) )
      del4( 1,ny) = dxm2*(  del2( nx,ny) - 2.0d0*del2( 1, ny )
     &                    + del2( 2 ,ny)
     &                    + del2( 1 ,ny) - 2.0d0*del2( 1,ny-1)
     &                    + del2(1,ny-2) )
      do i=2,nx-1
        del4(i, 1) = dxm2*(  del2(i-1, 1) - 2.0d0*del2( i,  1 )
     &                     + del2(i+1, 1)
     &                     + del2( i , 3) - 2.0d0*del2( i,  2 )
     &                     + del2(i,  1 ) )
        del4(i,ny) = dxm2*(  del2(i-1,ny) - 2.0d0*del2( i, ny )
     &                     + del2(i+1,ny)
     &                     + del2( i ,ny) - 2.0d0*del2( i,ny-1)
     &                     + del2(i,ny-2) )
      enddo
      del4(nx, 1) = dxm2*(  del2(nx-1, 1) - 2.0d0*del2(nx,  1 )
     &                    + del2( 1 , 1)
     &                    + del2( nx , 3) - 2.0d0*del2(nx,  2 )
     &                    + del2(nx,  1 ) )
      del4(nx,ny) = dxm2*(  del2(nx-1,ny) - 2.0d0*del2(nx, ny )
     &                    + del2( 1 ,ny)
     &                    + del2( nx ,ny) - 2.0d0*del2(nx,ny-1)
     &                    + del2(nx,ny-2) )
!$OMP END SINGLE
*
      END SUBROUTINE del4ch
*
c***********************************************************************
*
      SUBROUTINE genint (val, nx, ny, facwe, facsn, answer)
*
*     Computes area integral of internal values of field
*     val(nx,ny). facwe controls contribution from western
*     & eastern boundaries; facsn controls contribution from
*     southern & northern boundaries. Integral returned as answer
*     Modified version with reduced accumulator error
*
*     N.B. Parallelised using "orphaned directives";
*     only intended for use in the parallel regions
*     established earlier in calling subroutine monnc_comp

      IMPLICIT NONE

*     Subroutine arguments
      integer nx,ny
      double precision val(nx,ny),facwe,facsn,answer
*
*     Local variables
      integer i,j
      double precision sumi,xxs,xxn
*
!$OMP SINGLE
      answer = 0.0d0
!$OMP END SINGLE
*     Inner points + facwe*(W & E boundaries)
!$OMP DO SCHEDULE (STATIC)
!$OMP&   PRIVATE (i,j,sumi)
!$OMP&   REDUCTION (+:answer)
      do j=2,ny-1
        sumi = facwe*val(1,j)
        do i=2,nx-1
          sumi = sumi + val(i,j)
        enddo
        sumi = sumi + facwe*val(nx,j)
        answer = answer + sumi
      enddo
!$OMP END DO

!$OMP SINGLE
*     N & S boundary contributions
      xxs = facwe*val(1, 1)
      xxn = facwe*val(1,ny)
*     Inner points
      do i=2,nx-1
        xxs = xxs + val(i, 1)
        xxn = xxn + val(i,ny)
      enddo
      xxs = xxs + facwe*val(nx, 1)
      xxn = xxn + facwe*val(nx,ny)

      answer = answer + facsn*( xxs + xxn )
!$OMP END SINGLE

      END SUBROUTINE genint
*
c***********************************************************************

#ifndef ocean_only

      SUBROUTINE courat
*
*     Computes the maximum Courant number in each atmospheric
*     layer (mixed layer and all Q-G layers), to check the
*     likely quality of the temporal resolution of the flow.
*     The formulation of the velocity components follows that
*     in the AST advection in subroutine amladf in amlsubs.F,
*     since the advection scheme is likely to give problems if
*     the Courant number is too large. Evaluated at T points.
*     Also computes the min and max u- and v-velocity components
*     in all the layers, evaluated naturally at the face centres.

*     Modules
      USE parameters, ONLY : nxta, nyta, nxpa, nypa, nla
      USE atconst, ONLY : rdxaf0, hdxam1, dta
      USE intrfac, ONLY : uekat, vekat
      USE atstate, ONLY : pa
      USE monitor, ONLY : umminat, ummaxat, vmminat, vmmaxat, cnmlat,
     &                    ugminat, ugmaxat, vgminat, vgmaxat, cnqgat

      IMPLICIT NONE

*     Subroutine arguments
*
*     Local parameters
*
*     Local variables
      double precision umin,umax,vmin,vmax,vsqmax,um,up,vm,vp,velsqd
      integer i,j,k

*     Atmospheric mixed layer
*     =======================
      umin =  1.0d30
      umax = -1.0d30
      vmin =  1.0d30
      vmax = -1.0d30
      vsqmax = -1.0d30

!$OMP PARALLEL DEFAULT (NONE)
!$OMP&         PRIVATE (i,j,um,up,vm,vp,velsqd)
!$OMP&         SHARED  (rdxaf0,pa,uekat,vekat)
!$OMP&         REDUCTION(min:umin,vmin)
!$OMP&         REDUCTION(max:umax,vmax,vsqmax)

!$OMP DO SCHEDULE (STATIC)
      do j=2,nyta-1

*       Initialise the recurrence relation on u by
*       assigning the Western boundary (i=1) u value to up
        up = -rdxaf0*( pa(1,j+1,1) - pa(1, j ,1) ) + uekat(1, j )
        umin = min( umin, up )
        umax = max( umax, up )

*       Western bdy, inner points + Eastern bdy (standard case)
*       -------------------------------------------------------
        do i=1,nxta
          um = up
          up = -rdxaf0*( pa(i+1,j+1,1) - pa(i+1,j,1) ) + uekat(i+1,j)
          vm =  rdxaf0*( pa(i+1, j ,1) - pa(i, j ,1) ) + vekat(i, j )
          vp =  rdxaf0*( pa(i+1,j+1,1) - pa(i,j+1,1) ) + vekat(i,j+1)
          umin = min( umin, up )
          umax = max( umax, up )
          vmin = min( vmin, vm, vp )
          vmax = max( vmax, vm, vp )
          velsqd = (um + up)**2 + (vm + vp)**2
          vsqmax = max( vsqmax, velsqd )
        enddo

      enddo
!$OMP END DO NOWAIT

!$OMP END PARALLEL

*     Zonal boundaries, including corners
*     -----------------------------------
*     (N.B. do in recurrence form for simplicity, by removing
*      //ism on i-loop, unlike (currently) in amlsubs.F)

*     Initialise the recurrence relation on u by
*     assigning the Western boundary (i=1) u value to up
      up = -rdxaf0*( pa(1,2,1) - pa(1,1,1) ) + uekat(1,1)
      umin = min( umin, up )
      umax = max( umax, up )

*     Southern boundary: SW corner, inner points + SE corner (j=1)
*     ------------------------------------------------------------
      do i=1,nxta
        um = up
        up = -rdxaf0*( pa(i+1,2,1) - pa(i+1,1,1) ) + uekat(i+1,1)
        vm =  vekat(i,1)
        vp =  rdxaf0*( pa(i+1,2,1) - pa( i ,2,1) ) + vekat( i ,2)
        umin = min( umin, up )
        umax = max( umax, up )
        vmin = min( vmin, vm, vp )
        vmax = max( vmax, vm, vp )
        velsqd = (um + up)**2 + (vm + vp)**2
        vsqmax = max( vsqmax, velsqd )
      enddo

*     Initialise the recurrence relation on u by
*     assigning the Western boundary (i=1) u value to up
      up = -rdxaf0*( pa(1,nyta+1,1) - pa(1,nyta,1) ) + uekat(1,nyta)
      umin = min( umin, up )
      umax = max( umax, up )

*     Northern boundary: NW corner, inner points + NE corner (j=nyta)
*     ---------------------------------------------------------------
      do i=1,nxta
        um = up
        up = -rdxaf0*( pa(i+1,nyta+1,1) - pa(i+1,nyta,1) )
     &       + uekat(i+1,nyta)
        vm =  rdxaf0*( pa(i+1, nyta ,1) - pa( i ,nyta,1) )
     &       + vekat( i ,nyta)
        vp =  vekat(i,nyta+1)
        umin = min( umin, up )
        umax = max( umax, up )
        vmin = min( vmin, vm, vp )
        vmax = max( vmax, vm, vp )
        velsqd = (um + up)**2 + (vm + vp)**2
        vsqmax = max( vsqmax, velsqd )
      enddo

      umminat = umin
      ummaxat = umax
      vmminat = vmin
      vmmaxat = vmax
      cnmlat = hdxam1*dta*sqrt( vsqmax )

*     Atmosphere Q-G layers
*     =====================
      do k=1,nla
        umin =  1.0d30
        umax = -1.0d30
        vmin =  1.0d30
        vmax = -1.0d30
        vsqmax = -1.0d30

!$OMP   PARALLEL DEFAULT (NONE)
!$OMP&           PRIVATE (i,j,um,up,vm,vp,velsqd)
!$OMP&           SHARED  (k,rdxaf0,pa)
!$OMP&           REDUCTION(min:umin,vmin)
!$OMP&           REDUCTION(max:umax,vmax,vsqmax)

!$OMP   DO SCHEDULE (STATIC)
        do j=2,nyta-1

*         Initialise the recurrence relation on u by
*         assigning the Western boundary (i=1) u value to up
          up = -rdxaf0*( pa(1,j+1,k) - pa(1, j ,k) )
          umin = min( umin, up )
          umax = max( umax, up )

*         Western bdy, inner points + Eastern bdy (standard case)
*         -------------------------------------------------------
          do i=1,nxta
            um = up
            up = -rdxaf0*( pa(i+1,j+1,k) - pa(i+1,j,k) )
            vm =  rdxaf0*( pa(i+1, j ,k) - pa(i, j ,k) )
            vp =  rdxaf0*( pa(i+1,j+1,k) - pa(i,j+1,k) )
            umin = min( umin, up )
            umax = max( umax, up )
            vmin = min( vmin, vm, vp )
            vmax = max( vmax, vm, vp )
            velsqd = (um + up)**2 + (vm + vp)**2
            vsqmax = max( vsqmax, velsqd )
          enddo

        enddo
!$OMP   END DO NOWAIT

!$OMP   END PARALLEL

*       Zonal boundaries, including corners
*       -----------------------------------
*       (N.B. do in recurrence form for simplicity, by removing
*        //ism on i-loop, unlike (currently) in amlsubs.F)

*       Initialise the recurrence relation on u by
*       assigning the Western boundary (i=1) u value to up
        up = -rdxaf0*( pa(1,2,k) - pa(1,1,k) )
        umin = min( umin, up )
        umax = max( umax, up )

*       Southern boundary: SW corner, inner points + SE corner (j=1)
*       ------------------------------------------------------------
        do i=1,nxta
          um = up
          up = -rdxaf0*( pa(i+1,2,k) - pa(i+1,1,k) )
          vm = 0.0d0
          vp =  rdxaf0*( pa(i+1,2,k) - pa( i ,2,k) )
          umin = min( umin, up )
          umax = max( umax, up )
          vmin = min( vmin, vm, vp )
          vmax = max( vmax, vm, vp )
          velsqd = (um + up)**2 + (vm + vp)**2
          vsqmax = max( vsqmax, velsqd )
        enddo

*       Initialise the recurrence relation on u by
*       assigning the Western boundary (i=1) u value to up
        up = -rdxaf0*( pa(1,nypa,k) - pa(1,nyta,k) )
        umin = min( umin, up )
        umax = max( umax, up )

*       Northern boundary: NW corner, inner points + NE corner (j=nyta)
*       ---------------------------------------------------------------
        do i=1,nxta
          um = up
          up = -rdxaf0*( pa(i+1,nyta+1,k) - pa(i+1,nyta,k) )
          vm =  rdxaf0*( pa(i+1, nyta ,k) - pa( i ,nyta,k) )
          vp = 0.0d0
          umin = min( umin, up )
          umax = max( umax, up )
          vmin = min( vmin, vm, vp )
          vmax = max( vmax, vm, vp )
          velsqd = (um + up)**2 + (vm + vp)**2
          vsqmax = max( vsqmax, velsqd )
        enddo

        ugminat(k) = umin
        ugmaxat(k) = umax
        vgminat(k) = vmin
        vgmaxat(k) = vmax
        cnqgat(k) = hdxam1*dta*sqrt( vsqmax )

      enddo

      END SUBROUTINE courat

#endif /* not ocean_only */

c***********************************************************************

#ifndef atmos_only

      SUBROUTINE couroc
*
*     Computes the maximum Courant number in each ocean
*     layer (mixed layer and all Q-G layers), to check the
*     likely quality of the temporal resolution of the flow.
*     The formulation of the velocity components follows that
*     in the SST advection in subroutine omladf in omlsubs.F,
*     since the advection scheme is likely to give problems if
*     the Courant number is too large. Evaluated at T points.
*     Also computes the min and max u- and v-velocity components
*     in all the layers, evaluated naturally at the face centres.

*     Modules
      USE parameters, ONLY : nxto, nyto, nxpo, nypo, nlo, fnot
      USE occonst, ONLY : ycexp, rdxof0, hdxom1, dto
      USE intrfac, ONLY : tauxo, tauyo, hmoc
      USE ocstate, ONLY : po
      USE monitor, ONLY : umminoc, ummaxoc, vmminoc, vmmaxoc, cnmloc,
     &                    ugminoc, ugmaxoc, vgminoc, vgmaxoc, cnqgoc

      IMPLICIT NONE

*     Subroutine arguments
*
*     Local parameters
*
*     Local variables
      double precision uvgfac,rhf0hm,umin,umax,vmin,vmax,vsqmax,
     &                 um,up,vm,vp,velsqd
      integer i,j,k

*     Oceanic mixed layer
*     ===================
      uvgfac = ycexp*rdxof0
      rhf0hm = 0.5d0/(fnot*hmoc)
      umin =  1.0d30
      umax = -1.0d30
      vmin =  1.0d30
      vmax = -1.0d30
      vsqmax = -1.0d30

!$OMP PARALLEL DEFAULT (NONE)
!$OMP&         PRIVATE (i,j,um,up,vm,vp,velsqd)
!$OMP&         SHARED  (uvgfac,po,rhf0hm,tauyo,tauxo)
!$OMP&         REDUCTION(min:umin,vmin)
!$OMP&         REDUCTION(max:umax,vmax,vsqmax)

!$OMP DO SCHEDULE (STATIC)
      do j=2,nyto-1

*       Initialise the recurrence relation on u by
*       assigning the Western boundary (i=1) u value to up
#  ifdef cyclic_ocean
*       Zonally cyclic ocean
        up = -uvgfac*(  po(1,j+1,1) -  po(1,j,1) )
     &       +rhf0hm*( tauyo(1,j+1) + tauyo(1,j) )
#  else
*       Finite box ocean
*       Advection (no normal mass flux => u = 0)
        up = 0.0d0
#  endif
        umin = min( umin, up )
        umax = max( umax, up )

*       Western bdy + inner points (standard case)
*       ------------------------------------------
        do i=1,nxto-1
          um = up
          up = -uvgfac*(  po(i+1,j+1,1) -  po(i+1,j,1) )
     &         +rhf0hm*( tauyo(i+1,j+1) + tauyo(i+1,j) )
          vm =  uvgfac*(  po(i+1, j ,1) -  po(i, j ,1) )
     &         -rhf0hm*( tauxo(i+1, j ) + tauxo(i, j ) )
          vp =  uvgfac*(  po(i+1,j+1,1) -  po(i,j+1,1) )
     &         -rhf0hm*( tauxo(i+1,j+1) + tauxo(i,j+1) )
          umin = min( umin, up )
          umax = max( umax, up )
          vmin = min( vmin, vm, vp )
          vmax = max( vmax, vm, vp )
          velsqd = (um + up)**2 + (vm + vp)**2
          vsqmax = max( vsqmax, velsqd )
        enddo

*       Eastern boundary (i=nxto)
*       -------------------------
        um = up
#  ifdef cyclic_ocean
*       Zonally cyclic ocean
        up = -uvgfac*(  po(nxto+1,j+1,1) -  po(nxto+1,j,1) )
     &       +rhf0hm*( tauyo(nxto+1,j+1) + tauyo(nxto+1,j) )
#  else
*       Finite box ocean
*       Advection (no normal mass flux => up = 0)
        up = 0.0d0
#  endif
        vm =  uvgfac*(  po(nxto+1, j ,1) -  po(nxto, j ,1) )
     &       -rhf0hm*( tauxo(nxto+1, j ) + tauxo(nxto, j ) )
        vp =  uvgfac*(  po(nxto+1,j+1,1) -  po(nxto,j+1,1) )
     &       -rhf0hm*( tauxo(nxto+1,j+1) + tauxo(nxto,j+1) )
        umin = min( umin, up )
        umax = max( umax, up )
        vmin = min( vmin, vm, vp )
        vmax = max( vmax, vm, vp )
        velsqd = (um + up)**2 + (vm + vp)**2
        vsqmax = max( vsqmax, velsqd )

      enddo
!$OMP END DO NOWAIT

!$OMP END PARALLEL

*     Zonal boundaries, including corners (choice of BCs)
*     ---------------------------------------------------
*     (N.B. do in recurrence form for simplicity, by removing
*      //ism on i-loop, unlike (currently) in omlsubs.F)

*     Southwestern corner(i=1, j=1)
*     -----------------------------
#  ifdef cyclic_ocean
*     Zonally cyclic ocean
      um = -uvgfac*(  po(1,2,1) -  po(1,1,1) )
     &     +rhf0hm*( tauyo(1,2) + tauyo(1,1) )
#  else
*     Finite box ocean
*     Advection (no normal mass flux => um = 0)
      um = 0.0d0
#  endif
      up = -uvgfac*(  po(2,2,1) -  po(2,1,1) )
     &     +rhf0hm*( tauyo(2,2) + tauyo(2,1) )
      vp =  uvgfac*(  po(2,2,1) -  po(1,2,1) )
     &     -rhf0hm*( tauxo(2,2) + tauxo(1,2) )
#  ifdef sb_hflux
*     Advection consistent with an outflow across the
*     southern boundary equal to the Ekman transport.
*     p contribution to vm vanishes because p is uniform along bdy.
      vm = -rhf0hm*( tauxo(2,1) + tauxo(1,1) )
#  else
*     Advection (no normal mass flux => vm = 0)
      vm = 0.0d0
#  endif
      umin = min( umin, up )
      umax = max( umax, up )
      vmin = min( vmin, vm, vp )
      vmax = max( vmax, vm, vp )
      velsqd = (um + up)**2 + (vm + vp)**2
      vsqmax = max( vsqmax, velsqd )

*     Southern boundary inner points (j=1)
*     ------------------------------------
      do i=2,nxto-1
        um = up
        up = -uvgfac*(  po(i+1,2,1) -  po(i+1,1,1) )
     &       +rhf0hm*( tauyo(i+1,2) + tauyo(i+1,1) )
        vp =  uvgfac*(  po(i+1,2,1) -  po( i ,2,1) )
     &       -rhf0hm*( tauxo(i+1,2) + tauxo( i ,2) )
#  ifdef sb_hflux
*       Advection consistent with an outflow across the
*       southern boundary equal to the Ekman transport.
*       p contribution to vm vanishes because p is uniform along bdy.
        vm = -rhf0hm*( tauxo(i+1,1) + tauxo(i,1) )
#  else
*       Advection (no normal mass flux => vm = 0)
        vm = 0.0d0
#  endif
        umin = min( umin, up )
        umax = max( umax, up )
        vmin = min( vmin, vm, vp )
        vmax = max( vmax, vm, vp )
        velsqd = (um + up)**2 + (vm + vp)**2
        vsqmax = max( vsqmax, velsqd )
      enddo

*     Southeastern corner (i=nxto, j=1)
*     ---------------------------------
      um = up
#  ifdef cyclic_ocean
*     Zonally cyclic ocean
      up = -uvgfac*(  po(nxto+1,2,1) -  po(nxto+1,1,1) )
     &     +rhf0hm*( tauyo(nxto+1,2) + tauyo(nxto+1,1) )
#  else
*     Finite box ocean
*     Advection (no normal mass flux => up = 0)
      up = 0.0d0
#  endif
      vp =  uvgfac*(  po(nxto+1,2,1) -  po(nxto,2,1) )
     &     -rhf0hm*( tauxo(nxto+1,2) + tauxo(nxto,2) )
#  ifdef sb_hflux
*     Advection consistent with an outflow across the
*     southern boundary equal to the Ekman transport.
*     p contribution to vm vanishes because p is uniform along bdy.
      vm = -rhf0hm*( tauxo(nxto+1,1) + tauxo(nxto,1) )
#  else
*     Advection (no normal mass flux => vm = 0)
      vm = 0.0d0
#  endif
      umin = min( umin, up )
      umax = max( umax, up )
      vmin = min( vmin, vm, vp )
      vmax = max( vmax, vm, vp )
      velsqd = (um + up)**2 + (vm + vp)**2
      vsqmax = max( vsqmax, velsqd )

*     Northwestern corner (i=1, j=nyto)
*     ---------------------------------
#  ifdef cyclic_ocean
*     Zonally cyclic ocean
      um = -uvgfac*(  po(1,nyto+1,1) -  po(1,nyto,1) )
     &     +rhf0hm*( tauyo(1,nyto+1) + tauyo(1,nyto) )
#  else
*     Finite box ocean
*     Advection (no normal mass flux => um = 0)
      um = 0.0d0
#  endif
      up = -uvgfac*(  po(2,nyto+1,1) -  po(2,nyto,1) )
     &     +rhf0hm*( tauyo(2,nyto+1) + tauyo(2,nyto) )
      vm =  uvgfac*(  po(2,nyto,1) -  po(1,nyto,1) )
     &     -rhf0hm*( tauxo(2,nyto) + tauxo(1,nyto) )
#  ifdef nb_hflux
*     Advection consistent with an outflow across the
*     northern boundary equal to the Ekman transport,
*     carrying fluid of a specified temperature tnbdy.
*     p contribution to vp vanishes because p is uniform along bdy.
      vp = -rhf0hm*( tauxo(2,nyto+1) + tauxo(1,nyto+1) )
#  else
*     Advection (no normal mass flux => vp = 0)
      vp = 0.0d0
#  endif
      umin = min( umin, up )
      umax = max( umax, up )
      vmin = min( vmin, vm, vp )
      vmax = max( vmax, vm, vp )
      velsqd = (um + up)**2 + (vm + vp)**2
      vsqmax = max( vsqmax, velsqd )

*     Northern boundary inner points (j=nyto)
*     ---------------------------------------
      do i=2,nxto-1
        um = up
        up = -uvgfac*(  po(i+1,nyto+1,1) -  po(i+1,nyto,1) )
     &       +rhf0hm*( tauyo(i+1,nyto+1) + tauyo(i+1,nyto) )
        vm =  uvgfac*(  po(i+1,nyto,1) -  po( i ,nyto,1) )
     &       -rhf0hm*( tauxo(i+1,nyto) + tauxo( i ,nyto) )
#  ifdef nb_hflux
*       Advection consistent with an outflow across the
*       northern boundary equal to the Ekman transport.
*       p contribution to vp vanishes because p is uniform along bdy.
        vp = -rhf0hm*( tauxo(i+1,nyto+1) + tauxo(i,nyto+1) )
#  else
*       Advection (no normal mass flux => vp = 0)
        vp = 0.0d0
#  endif
        umin = min( umin, up )
        umax = max( umax, up )
        vmin = min( vmin, vm, vp )
        vmax = max( vmax, vm, vp )
        velsqd = (um + up)**2 + (vm + vp)**2
        vsqmax = max( vsqmax, velsqd )
      enddo

*     Northeastern corner (i=nxto, j=nyto)
*     ------------------------------------
      um = up
#  ifdef cyclic_ocean
*     Zonally cyclic ocean
      up = -uvgfac*(  po(nxto+1,nyto+1,1) -  po(nxto+1,nyto,1) )
     &     +rhf0hm*( tauyo(nxto+1,nyto+1) + tauyo(nxto+1,nyto) )
#  else
*     Finite box ocean
*     Advection (no normal mass flux => up = 0)
      up = 0.0d0
#  endif
      vm =  uvgfac*(  po(nxto+1,nyto,1) -  po(nxto,nyto,1) )
     &     -rhf0hm*( tauxo(nxto+1,nyto) + tauxo(nxto,nyto) )
#  ifdef nb_hflux
*     Advection consistent with an outflow across the
*     northern boundary equal to the Ekman transport.
*     p contribution to vp vanishes because p is uniform along bdy.
      vp = -rhf0hm*( tauxo(nxto+1,nyto+1) + tauxo(nxto,nyto+1) )
#  else
*     Advection (no normal mass flux => vp = 0)
      vp = 0.0d0
#  endif
      umin = min( umin, up )
      umax = max( umax, up )
      vmin = min( vmin, vm, vp )
      vmax = max( vmax, vm, vp )
      velsqd = (um + up)**2 + (vm + vp)**2
      vsqmax = max( vsqmax, velsqd )

      umminoc = umin
      ummaxoc = umax
      vmminoc = vmin
      vmmaxoc = vmax
      cnmloc = hdxom1*dto*sqrt( vsqmax )

*     Ocean Q-G layers
*     ================
      do k=1,nlo
        umin =  1.0d30
        umax = -1.0d30
        vmin =  1.0d30
        vmax = -1.0d30
        vsqmax = -1.0d30

!$OMP   PARALLEL DEFAULT (NONE)
!$OMP&           PRIVATE (i,j,um,up,vm,vp,velsqd)
!$OMP&           SHARED  (k,rdxof0,po)
!$OMP&           REDUCTION(min:umin,vmin)
!$OMP&           REDUCTION(max:umax,vmax,vsqmax)

!$OMP   DO SCHEDULE (STATIC)
        do j=2,nyto-1

*         Initialise the recurrence relation on u by
*         assigning the Western boundary (i=1) u value to up
#  ifdef cyclic_ocean
*         Zonally cyclic ocean
          up = -rdxof0*(  po(1,j+1,k) -  po(1,j,k) )
#  else
*         Finite box ocean
*         Advection (no normal mass flux => u = 0)
          up = 0.0d0
#  endif
          umin = min( umin, up )
          umax = max( umax, up )

*         Western bdy + inner points (standard case)
*         ------------------------------------------
          do i=1,nxto-1
            um = up
            up = -rdxof0*( po(i+1,j+1,k) - po(i+1,j,k) )
            vm =  rdxof0*( po(i+1, j ,k) - po(i, j ,k) )
            vp =  rdxof0*( po(i+1,j+1,k) - po(i,j+1,k) )
            umin = min( umin, up )
            umax = max( umax, up )
            vmin = min( vmin, vm, vp )
            vmax = max( vmax, vm, vp )
            velsqd = (um + up)**2 + (vm + vp)**2
            vsqmax = max( vsqmax, velsqd )
          enddo

*         Eastern boundary (i=nxto)
*         -------------------------
          um = up
#  ifdef cyclic_ocean
*         Zonally cyclic ocean
          up = -rdxof0*( po(nxto+1,j+1,k) - po(nxto+1,j,k) )
#  else
*         Finite box ocean
*         Advection (no normal mass flux => up = 0)
          up = 0.0d0
#  endif
          vm =  rdxof0*( po(nxto+1, j ,k) - po(nxto, j ,k) )
          vp =  rdxof0*( po(nxto+1,j+1,k) - po(nxto,j+1,k) )
          umin = min( umin, up )
          umax = max( umax, up )
          vmin = min( vmin, vm, vp )
          vmax = max( vmax, vm, vp )
          velsqd = (um + up)**2 + (vm + vp)**2
          vsqmax = max( vsqmax, velsqd )

        enddo
!$OMP   END DO NOWAIT

!$OMP   END PARALLEL

*       Zonal boundaries, including corners
*       -----------------------------------
*       (N.B. do in recurrence form for simplicity, by removing
*        //ism on i-loop, unlike (currently) in omlsubs.F)

*       Initialise the recurrence relation on u by
*       assigning the Western boundary (i=1) u value to up
#  ifdef cyclic_ocean
*       Zonally cyclic ocean
        up = -rdxof0*( po(1,2,k) - po(1,1,k) )
#  else
*       Finite box ocean
*       Advection (no normal mass flux => u = 0)
        up = 0.0d0
#  endif
        umin = min( umin, up )
        umax = max( umax, up )

*       Southern boundary: SW corner + inner points (j=1)
*       -------------------------------------------------
        do i=1,nxto-1
          um = up
          up = -rdxof0*( po(i+1,2,k) - po(i+1,1,k) )
          vp =  rdxof0*( po(i+1,2,k) - po( i ,2,k) )
*         Advection (no normal mass flux => vm = 0)
          vm = 0.0d0
          umin = min( umin, up )
          umax = max( umax, up )
          vmin = min( vmin, vm, vp )
          vmax = max( vmax, vm, vp )
          velsqd = (um + up)**2 + (vm + vp)**2
          vsqmax = max( vsqmax, velsqd )
        enddo

*       Southeastern corner (i=nxto, j=1)
*       ---------------------------------
        um = up
#  ifdef cyclic_ocean
*       Zonally cyclic ocean
        up = -rdxof0*( po(nxto+1,2,k) - po(nxto+1,1,k) )
#  else
*       Finite box ocean
*       Advection (no normal mass flux => up = 0)
        up = 0.0d0
#  endif
        vp =  rdxof0*( po(nxto+1,2,k) - po(nxto,2,k) )
*       Advection (no normal mass flux => vm = 0)
        vm = 0.0d0
        umin = min( umin, up )
        umax = max( umax, up )
        vmin = min( vmin, vm, vp )
        vmax = max( vmax, vm, vp )
        velsqd = (um + up)**2 + (vm + vp)**2
        vsqmax = max( vsqmax, velsqd )

*       Initialise the recurrence relation on u by
*       assigning the Western boundary (i=1) u value to up
#  ifdef cyclic_ocean
*       Zonally cyclic ocean
        up = -rdxof0*( po(1,nyto+1,k) - po(1,nyto,k) )
#  else
*       Finite box ocean
*       Advection (no normal mass flux => u = 0)
        up = 0.0d0
#  endif
        umin = min( umin, up )
        umax = max( umax, up )

*       Northern boundary: NW corner + inner points (j=nyto)
*       ----------------------------------------------------
        do i=1,nxto-1
          um = up
          up = -rdxof0*( po(i+1,nyto+1,k) - po(i+1,nyto,k) )
          vm =  rdxof0*( po(i+1, nyto ,k) - po( i ,nyto,k) )
*         Advection (no normal mass flux => vp = 0)
          vp = 0.0d0
          umin = min( umin, up )
          umax = max( umax, up )
          vmin = min( vmin, vm, vp )
          vmax = max( vmax, vm, vp )
          velsqd = (um + up)**2 + (vm + vp)**2
          vsqmax = max( vsqmax, velsqd )
        enddo

*       Northeastern corner (i=nxto, j=nyto)
*       ------------------------------------
        um = up
#  ifdef cyclic_ocean
*       Zonally cyclic ocean
        up = -rdxof0*( po(nxto+1,nyto+1,k) - po(nxto+1,nyto,k) )
#  else
*       Finite box ocean
*       Advection (no normal mass flux => up = 0)
        up = 0.0d0
#  endif
        vm =  rdxof0*( po(nxto+1,nyto,k) - po(nxto,nyto,k) )
*       Advection (no normal mass flux => vp = 0)
        vp = 0.0d0
        umin = min( umin, up )
        umax = max( umax, up )
        vmin = min( vmin, vm, vp )
        vmax = max( vmax, vm, vp )
        velsqd = (um + up)**2 + (vm + vp)**2
        vsqmax = max( vsqmax, velsqd )

        ugminoc(k) = umin
        ugmaxoc(k) = umax
        vgminoc(k) = vmin
        vgmaxoc(k) = vmax
        cnqgoc(k) = hdxom1*dto*sqrt( vsqmax )

      enddo

      END SUBROUTINE couroc

#endif /* not atmos_only */

c***********************************************************************

      SUBROUTINE monnc_init

*     Modules
      USE parameters, ONLY : nla, nlo
#ifdef use_netcdf
      USE atconst, ONLY : hat
      USE occonst, ONLY : hoc
#endif /* use_netcdf */
      USE timinfo, ONLY : numoutsteps
      USE nc_subs, ONLY : handle_err

      IMPLICIT NONE

#ifdef use_netcdf
      INCLUDE 'netcdf.inc'
*
*     Subroutine arguments

*     Local parameters
      character (len=*), parameter :: subnam = 'monnc_init'

!! Netcdf variables used locally
      integer ncstat, timedim, ladim, lodim, lamdim, lomdim
      integer edims(2), lo_id, la_id, lom_id, lam_id, k
      double precision tomp(nlo), tamp(nla)

!! Define time dimension:
      ncstat = nf_def_dim(monncid, 'time', numoutsteps, timedim)
      if ( ncstat.ne.NF_NOERR ) call handle_err (ncstat, subnam)

!! Define the zo dimension
#  ifndef atmos_only
      ncstat = nf_def_dim(monncid, 'zo', nlo, lodim)
      if ( ncstat.ne.NF_NOERR ) call handle_err (ncstat, subnam)
!! Define the zom dimension (with length nlo-1)
      ncstat = nf_def_dim(monncid, 'zom', nlo-1, lomdim)
      if ( ncstat.ne.NF_NOERR ) call handle_err (ncstat, subnam)
#  endif /* not atmos_only */

!! Define the za dimension
#  ifndef ocean_only
      ncstat = nf_def_dim(monncid, 'za', nla, ladim)
      if ( ncstat.ne.NF_NOERR ) call handle_err (ncstat, subnam)
!! Define the zam dimension (with length nla-1)
      ncstat = nf_def_dim(monncid, 'zam', nla-1, lamdim)
      if ( ncstat.ne.NF_NOERR ) call handle_err (ncstat, subnam)
#  endif /* not ocean_only */

!! Define a one-dimensional variable called 'time' which
!! stores time data in years
      ncstat=nf_def_var(monncid, 'time', NF_FLOAT, 1, timedim, time_id)
      if ( ncstat.ne.NF_NOERR ) call handle_err (ncstat, subnam)
      ncstat = nf_put_att_text(monncid, time_id, 'units', 5, 'years')
      if ( ncstat.ne.NF_NOERR ) call handle_err (ncstat, subnam)
      ncstat = nf_put_att_text(monncid, time_id,'long_name',4,'Time')
      if ( ncstat.ne.NF_NOERR ) call handle_err (ncstat, subnam)

!! Define a 1-dimensional variable called 'zo'
#  ifndef atmos_only
      ncstat = nf_def_var(monncid, 'zo', NF_FLOAT, 1, lodim, lo_id)
      if ( ncstat.ne.NF_NOERR ) call handle_err (ncstat, subnam)
      ncstat = nf_put_att_text(monncid, lo_id, 'units', 2, 'km')
      if ( ncstat.ne.NF_NOERR ) call handle_err (ncstat, subnam)
      ncstat = nf_put_att_text(monncid, lo_id, 'long_name',
     &            26, 'Ocean mid-layer depth axis')
      if ( ncstat.ne.NF_NOERR ) call handle_err (ncstat, subnam)
#  endif /* not atmos_only */

!! Define a 1-dimensional variable called 'zom'
#  ifndef atmos_only
      ncstat = nf_def_var(monncid, 'zom', NF_FLOAT, 1, lomdim, lom_id)
      if ( ncstat.ne.NF_NOERR ) call handle_err (ncstat, subnam)
      ncstat = nf_put_att_text(monncid, lom_id, 'units', 2, 'km')
      if ( ncstat.ne.NF_NOERR ) call handle_err (ncstat, subnam)
      ncstat = nf_put_att_text(monncid, lom_id, 'long_name',
     &            26, 'Ocean interface depth axis')
      if ( ncstat.ne.NF_NOERR ) call handle_err (ncstat, subnam)
#  endif /* not atmos_only */

!! Define a 1-dimensional variable called 'za'
#  ifndef ocean_only
      ncstat = nf_def_var(monncid, 'za', NF_FLOAT, 1, ladim, la_id)
      if ( ncstat.ne.NF_NOERR ) call handle_err (ncstat, subnam)
      ncstat = nf_put_att_text(monncid, la_id, 'units', 2, 'km')
      if ( ncstat.ne.NF_NOERR ) call handle_err (ncstat, subnam)
      ncstat = nf_put_att_text(monncid, la_id, 'long_name',
     &            22, 'Atmosphere height axis')
      if ( ncstat.ne.NF_NOERR ) call handle_err (ncstat, subnam)
#  endif /* not ocean_only */

!! Define a 1-dimensional variable called 'zam'
#  ifndef ocean_only
      ncstat = nf_def_var(monncid, 'zam', NF_FLOAT, 1, lamdim, lam_id)
      if ( ncstat.ne.NF_NOERR ) call handle_err (ncstat, subnam)
      ncstat = nf_put_att_text(monncid, lam_id, 'units', 2, 'km')
      if ( ncstat.ne.NF_NOERR ) call handle_err (ncstat, subnam)
      ncstat = nf_put_att_text(monncid, lam_id, 'long_name',
     &            22, 'Atmosphere height axis')
      if ( ncstat.ne.NF_NOERR ) call handle_err (ncstat, subnam)
#  endif /* not ocean_only */

#  ifndef ocean_only
*     Atmosphere diagnostics from monnc_comp
*     --------------------------------------
!! Define a one-dimensional variable called 'wetmat'
      ncstat=nf_def_var(monncid,'wetmat',NF_FLOAT,1,timedim,wetmat_id)
      if ( ncstat.ne.NF_NOERR ) call handle_err (ncstat, subnam)
      ncstat = nf_put_att_text(monncid, wetmat_id, 'units', 3, 'm/s')
      if ( ncstat.ne.NF_NOERR ) call handle_err (ncstat, subnam)
      ncstat = nf_put_att_text(monncid, wetmat_id, 'long_name',
     &            43, 'Average atmospheric Ekman velocity (T-grid)')
      if ( ncstat.ne.NF_NOERR ) call handle_err (ncstat, subnam)

!! Define a one-dimensional variable called 'watmat'
      ncstat=nf_def_var(monncid,'watmat',NF_FLOAT,1,timedim,watmat_id)
      if ( ncstat.ne.NF_NOERR ) call handle_err (ncstat, subnam)
      ncstat = nf_put_att_text(monncid, watmat_id, 'units', 3, 'm/s')
      if ( ncstat.ne.NF_NOERR ) call handle_err (ncstat, subnam)
      ncstat = nf_put_att_text(monncid, watmat_id, 'long_name',
     &            44, 'Absolute atmospheric Ekman velocity (T-grid)')
      if ( ncstat.ne.NF_NOERR ) call handle_err (ncstat, subnam)

!! Define a one-dimensional variable called 'wepmat'
      ncstat=nf_def_var(monncid,'wepmat',NF_FLOAT,1,timedim,wepmat_id)
      if ( ncstat.ne.NF_NOERR ) call handle_err (ncstat, subnam)
      ncstat = nf_put_att_text(monncid, wepmat_id, 'units', 3, 'm/s')
      if ( ncstat.ne.NF_NOERR ) call handle_err (ncstat, subnam)
      ncstat = nf_put_att_text(monncid, wepmat_id, 'long_name',
     &            43, 'Average atmospheric Ekman velocity (p-grid)')
      if ( ncstat.ne.NF_NOERR ) call handle_err (ncstat, subnam)

!! Define a one-dimensional variable called 'wapmat'
      ncstat=nf_def_var(monncid,'wapmat',NF_FLOAT,1,timedim,wapmat_id)
      if ( ncstat.ne.NF_NOERR ) call handle_err (ncstat, subnam)
      ncstat = nf_put_att_text(monncid, wapmat_id, 'units', 3, 'm/s')
      if ( ncstat.ne.NF_NOERR ) call handle_err (ncstat, subnam)
      ncstat = nf_put_att_text(monncid, wapmat_id, 'long_name',
     &            44, 'Absolute atmospheric Ekman velocity (p-grid)')
      if ( ncstat.ne.NF_NOERR ) call handle_err (ncstat, subnam)

!! Define a two-dimensional variable called 'entmat'
      edims(1) = lamdim
      edims(2) = timedim
      ncstat = nf_def_var(monncid,'entmat',NF_FLOAT,2,edims,entmat_id)
      if ( ncstat.ne.NF_NOERR ) call handle_err (ncstat, subnam)
      ncstat = nf_put_att_text(monncid, entmat_id, 'units', 3, 'm/s')
      if ( ncstat.ne.NF_NOERR ) call handle_err (ncstat, subnam)
      ncstat = nf_put_att_text(monncid, entmat_id, 'long_name',
     &            31, 'Average atmospheric entrainment')
      if ( ncstat.ne.NF_NOERR ) call handle_err (ncstat, subnam)

!! Define a two-dimensional variable called 'enamat'
      edims(1) = lamdim
      edims(2) = timedim
      ncstat = nf_def_var(monncid,'enamat',NF_FLOAT,2,edims,enamat_id)
      if ( ncstat.ne.NF_NOERR ) call handle_err (ncstat, subnam)
      ncstat = nf_put_att_text(monncid, enamat_id, 'units', 3, 'm/s')
      if ( ncstat.ne.NF_NOERR ) call handle_err (ncstat, subnam)
      ncstat = nf_put_att_text(monncid, enamat_id, 'long_name',
     &            32, 'Absolute atmospheric entrainment')
      if ( ncstat.ne.NF_NOERR ) call handle_err (ncstat, subnam)

!! Define a two-dimensional variable called 'etamat'
      edims(1) = lamdim
      edims(2) = timedim
      ncstat = nf_def_var(monncid,'etamat',NF_FLOAT,2,edims,etamat_id)
      if ( ncstat.ne.NF_NOERR ) call handle_err (ncstat, subnam)
      ncstat = nf_put_att_text(monncid, etamat_id, 'units', 1, 'm')
      if ( ncstat.ne.NF_NOERR ) call handle_err (ncstat, subnam)
      ncstat = nf_put_att_text(monncid, etamat_id, 'long_name',
     &            36, 'Average atmospheric interface height')
      if ( ncstat.ne.NF_NOERR ) call handle_err (ncstat, subnam)

!! Define a two-dimensional variable called 'ddtkeat'
      edims(1) = ladim
      edims(2) = timedim
      ncstat = nf_def_var(monncid,'ddtkeat',NF_FLOAT,2,edims,ddtkeat_id)
      if ( ncstat.ne.NF_NOERR ) call handle_err (ncstat, subnam)
      ncstat = nf_put_att_text(monncid, ddtkeat_id, 'units',5,'W/m^2')
      if ( ncstat.ne.NF_NOERR ) call handle_err (ncstat, subnam)
      ncstat = nf_put_att_text(monncid, ddtkeat_id, 'long_name',
     &            32, 'd/dt(average KE in atmos. layer)')
      if ( ncstat.ne.NF_NOERR ) call handle_err (ncstat, subnam)

!! Define a two-dimensional variable called 'ddtpeat'
      edims(1) = lamdim
      edims(2) = timedim
      ncstat = nf_def_var(monncid,'ddtpeat',NF_FLOAT,2,edims,ddtpeat_id)
      if ( ncstat.ne.NF_NOERR ) call handle_err (ncstat, subnam)
      ncstat = nf_put_att_text(monncid, ddtpeat_id, 'units',5,'W/m^2')
      if ( ncstat.ne.NF_NOERR ) call handle_err (ncstat, subnam)
      ncstat = nf_put_att_text(monncid, ddtpeat_id, 'long_name',
     &            36, 'd/dt(average PE at atmos. interface)')
      if ( ncstat.ne.NF_NOERR ) call handle_err (ncstat, subnam)

!! Define a two-dimensional variable called 'et2mat'
      edims(1) = lamdim
      edims(2) = timedim
      ncstat = nf_def_var(monncid,'et2mat',NF_FLOAT,2,edims,et2mat_id)
      if ( ncstat.ne.NF_NOERR ) call handle_err (ncstat, subnam)
      ncstat = nf_put_att_text(monncid, et2mat_id, 'units',5,'W/m^2')
      if ( ncstat.ne.NF_NOERR ) call handle_err (ncstat, subnam)
      ncstat = nf_put_att_text(monncid, et2mat_id, 'long_name',
     &            25, 'Average atmospheric eta^2')
      if ( ncstat.ne.NF_NOERR ) call handle_err (ncstat, subnam)

!! Define a two-dimensional variable called 'pkenatat'
      edims(1) = lamdim
      edims(2) = timedim
      ncstat = nf_def_var(monncid,'pkenat',NF_FLOAT,2,edims,pkenat_id)
      if ( ncstat.ne.NF_NOERR ) call handle_err (ncstat, subnam)
      ncstat = nf_put_att_text(monncid, pkenat_id, 'units',5,'W/m^2')
      if ( ncstat.ne.NF_NOERR ) call handle_err (ncstat, subnam)
      ncstat = nf_put_att_text(monncid, pkenat_id, 'long_name',
     &            25, 'Average atmospheric eta*e')
      if ( ncstat.ne.NF_NOERR ) call handle_err (ncstat, subnam)

!! Define a one-dimensional variable called 'utauat'
      ncstat=nf_def_var(monncid,'utauat',NF_FLOAT,1,timedim,utauat_id)
      if ( ncstat.ne.NF_NOERR ) call handle_err (ncstat, subnam)
      ncstat = nf_put_att_text(monncid, utauat_id, 'units',5,'W/m^2')
      if ( ncstat.ne.NF_NOERR ) call handle_err (ncstat, subnam)
      ncstat = nf_put_att_text(monncid, utauat_id, 'long_name',
     &            25, 'Average atmospheric p1*we')
      if ( ncstat.ne.NF_NOERR ) call handle_err (ncstat, subnam)

!! Define a two-dimensional variable called 'atstpos'
      edims(1) = ladim
      edims(2) = timedim
      ncstat = nf_def_var(monncid,'atstpos',NF_INT,2,edims,atstpos_id)
      if ( ncstat.ne.NF_NOERR ) call handle_err (ncstat, subnam)
      ncstat=nf_put_att_text(monncid,atstpos_id,'units',10,'gridsquare')
      if ( ncstat.ne.NF_NOERR ) call handle_err (ncstat, subnam)
      ncstat = nf_put_att_text(monncid,atstpos_id,'long_name',53,
     &        'Position of maximum atmospheric velocity (gridsquare)')
      if ( ncstat.ne.NF_NOERR ) call handle_err (ncstat, subnam)

!! Define a two-dimensional variable called 'atstval'
      edims(1) = ladim
      edims(2) = timedim
      ncstat=nf_def_var(monncid,'atstval',NF_FLOAT,2,edims,atstval_id)
      if ( ncstat.ne.NF_NOERR ) call handle_err (ncstat, subnam)
      ncstat = nf_put_att_text(monncid, atstval_id, 'units', 3, 'm/s')
      if ( ncstat.ne.NF_NOERR ) call handle_err (ncstat, subnam)
      ncstat = nf_put_att_text(monncid, atstval_id, 'long_name',
     &            34, 'Maximum atmospheric velocity (m/s)')
      if ( ncstat.ne.NF_NOERR ) call handle_err (ncstat, subnam)

!! Define a two-dimensional variable called 'pavgat'
      edims(1) = ladim
      edims(2) = timedim
      ncstat = nf_def_var(monncid,'pavgat',NF_FLOAT,2,edims,pavgat_id)
      if ( ncstat.ne.NF_NOERR ) call handle_err (ncstat, subnam)
      ncstat = nf_put_att_text(monncid,pavgat_id,'units',8,'m^2 s^-2')
      if ( ncstat.ne.NF_NOERR ) call handle_err (ncstat, subnam)
      ncstat = nf_put_att_text(monncid,pavgat_id,'long_name',
     &            36, 'Average atmospheric dynamic pressure')
      if ( ncstat.ne.NF_NOERR ) call handle_err (ncstat, subnam)

!! Define a two-dimensional variable called 'qavgat'
      edims(1) = ladim
      edims(2) = timedim
      ncstat = nf_def_var(monncid,'qavgat',NF_FLOAT,2,edims,qavgat_id)
      if ( ncstat.ne.NF_NOERR ) call handle_err (ncstat, subnam)
      ncstat = nf_put_att_text(monncid, qavgat_id, 'units',4,'s^-1')
      if ( ncstat.ne.NF_NOERR ) call handle_err (ncstat, subnam)
      ncstat = nf_put_att_text(monncid, qavgat_id, 'long_name',
     &            29, 'Average atmospheric vorticity')
      if ( ncstat.ne.NF_NOERR ) call handle_err (ncstat, subnam)

!! Define a two-dimensional variable called 'ah4dat'
      edims(1) = ladim
      edims(2) = timedim
      ncstat = nf_def_var(monncid,'ah4dat',NF_FLOAT,2,edims,ah4dat_id)
      if ( ncstat.ne.NF_NOERR ) call handle_err (ncstat, subnam)
      ncstat = nf_put_att_text(monncid, ah4dat_id, 'units',5,'W/m^2')
      if ( ncstat.ne.NF_NOERR ) call handle_err (ncstat, subnam)
      ncstat = nf_put_att_text(monncid, ah4dat_id, 'long_name',
     &            31, 'Average atmospheric dissipation')
      if ( ncstat.ne.NF_NOERR ) call handle_err (ncstat, subnam)

!! Define a two-dimensional variable called 'kealat'
      edims(1) = ladim
      edims(2) = timedim
      ncstat = nf_def_var(monncid,'kealat',NF_FLOAT,2,edims,kealat_id)
      if ( ncstat.ne.NF_NOERR ) call handle_err (ncstat, subnam)
      ncstat = nf_put_att_text(monncid, kealat_id, 'units',5,'J/m^2')
      if ( ncstat.ne.NF_NOERR ) call handle_err (ncstat, subnam)
      ncstat = nf_put_att_text(monncid, kealat_id, 'long_name',
     &            40, 'Average atmospheric layer kinetic energy')
      if ( ncstat.ne.NF_NOERR ) call handle_err (ncstat, subnam)

!! Define a one-dimensional variable called 'tmlmat'
      ncstat=nf_def_var(monncid,'tmlmat',NF_FLOAT,1,timedim,tmlmat_id)
      if ( ncstat.ne.NF_NOERR ) call handle_err (ncstat, subnam)
      ncstat = nf_put_att_text(monncid, tmlmat_id, 'units', 1, 'K')
      if ( ncstat.ne.NF_NOERR ) call handle_err (ncstat, subnam)
      ncstat = nf_put_att_text(monncid, tmlmat_id, 'long_name',
     &            43, 'Average atmospheric mixed layer temperature')
      if ( ncstat.ne.NF_NOERR ) call handle_err (ncstat, subnam)

!! Define one-dimensional variables called 'astmin' and 'astmax'
      ncstat=nf_def_var(monncid,'astmin',NF_FLOAT,1,timedim,astmin_id)
      if ( ncstat.ne.NF_NOERR ) call handle_err (ncstat, subnam)
      ncstat = nf_put_att_text(monncid, astmin_id, 'units', 1, 'K')
      if ( ncstat.ne.NF_NOERR ) call handle_err (ncstat, subnam)
      ncstat = nf_put_att_text(monncid, astmin_id, 'long_name',
     &            40, 'Min. atmos. mixed layer temperature')
      if ( ncstat.ne.NF_NOERR ) call handle_err (ncstat, subnam)
      ncstat=nf_def_var(monncid,'astmax',NF_FLOAT,1,timedim,astmax_id)
      if ( ncstat.ne.NF_NOERR ) call handle_err (ncstat, subnam)
      ncstat = nf_put_att_text(monncid, astmax_id, 'units', 1, 'K')
      if ( ncstat.ne.NF_NOERR ) call handle_err (ncstat, subnam)
      ncstat = nf_put_att_text(monncid, astmax_id, 'long_name',
     &            40, 'Max. atmospheric mixed layer temperature')
      if ( ncstat.ne.NF_NOERR ) call handle_err (ncstat, subnam)

!! Define a one-dimensional variable called 'hmlmat'
      ncstat=nf_def_var(monncid,'hmlmat',NF_FLOAT,1,timedim,hmlmat_id)
      if ( ncstat.ne.NF_NOERR ) call handle_err (ncstat, subnam)
      ncstat = nf_put_att_text(monncid, hmlmat_id, 'units', 1, 'm')
      if ( ncstat.ne.NF_NOERR ) call handle_err (ncstat, subnam)
      ncstat = nf_put_att_text(monncid, hmlmat_id, 'long_name',
     &            41, 'Average atmospheric mixed layer thickness')
      if ( ncstat.ne.NF_NOERR ) call handle_err (ncstat, subnam)

!! Define a one-dimensional variable called 'hcmlat'
      ncstat=nf_def_var(monncid,'hcmlat',NF_FLOAT,1,timedim,hcmlat_id)
      if ( ncstat.ne.NF_NOERR ) call handle_err (ncstat, subnam)
      ncstat = nf_put_att_text(monncid, hcmlat_id, 'units',5,'K.m^3')
      if ( ncstat.ne.NF_NOERR ) call handle_err (ncstat, subnam)
      ncstat = nf_put_att_text(monncid, hcmlat_id, 'long_name',
     &            45, 'Total heat content of atmospheric mixed layer')
      if ( ncstat.ne.NF_NOERR ) call handle_err (ncstat, subnam)

!! Define a one-dimensional variable called 'tmaooc'
      ncstat=nf_def_var(monncid,'tmaooc',NF_FLOAT,1,timedim,tmaooc_id)
      if ( ncstat.ne.NF_NOERR ) call handle_err (ncstat, subnam)
      ncstat = nf_put_att_text(monncid, tmaooc_id, 'units', 1, 'K')
      if ( ncstat.ne.NF_NOERR ) call handle_err (ncstat, subnam)
      ncstat = nf_put_att_text(monncid, tmaooc_id, 'long_name', 51,
     &           'Mean atmospheric mixed layer temperature over ocean')
      if ( ncstat.ne.NF_NOERR ) call handle_err (ncstat, subnam)

!! Define a one-dimensional variable called 'olrtop'
      ncstat = nf_def_var(monncid,'olrtop',NF_FLOAT,1,timedim,olrtop_id)
      if ( ncstat.ne.NF_NOERR ) call handle_err (ncstat, subnam)
      ncstat = nf_put_att_text(monncid, olrtop_id, 'units',5,'W/m^2')
      if ( ncstat.ne.NF_NOERR ) call handle_err (ncstat, subnam)
      ncstat = nf_put_att_text(monncid, olrtop_id, 'long_name',
     &            28, 'Outgoing longwave radiation')
      if ( ncstat.ne.NF_NOERR ) call handle_err (ncstat, subnam)

!! Define one-dimensional variables called 'umminat' and 'ummaxat'
      ncstat=nf_def_var(monncid,'umminat',NF_FLOAT,1,timedim,umminat_id)
      if ( ncstat.ne.NF_NOERR ) call handle_err (ncstat, subnam)
      ncstat = nf_put_att_text(monncid, umminat_id, 'units', 3, 'm/s')
      if ( ncstat.ne.NF_NOERR ) call handle_err (ncstat, subnam)
      ncstat = nf_put_att_text(monncid, umminat_id, 'long_name',
     &            35, 'Minimum u in atmosphere mixed layer')
      if ( ncstat.ne.NF_NOERR ) call handle_err (ncstat, subnam)
      ncstat=nf_def_var(monncid,'ummaxat',NF_FLOAT,1,timedim,ummaxat_id)
      if ( ncstat.ne.NF_NOERR ) call handle_err (ncstat, subnam)
      ncstat = nf_put_att_text(monncid, ummaxat_id, 'units', 3, 'm/s')
      if ( ncstat.ne.NF_NOERR ) call handle_err (ncstat, subnam)
      ncstat = nf_put_att_text(monncid, ummaxat_id, 'long_name',
     &            35, 'Maximum u in atmosphere mixed layer')
      if ( ncstat.ne.NF_NOERR ) call handle_err (ncstat, subnam)

!! Define one-dimensional variables called 'vmminat' and 'vmmaxat'
      ncstat=nf_def_var(monncid,'vmminatat',NF_FLOAT,1,timedim,vmminat_id)
      if ( ncstat.ne.NF_NOERR ) call handle_err (ncstat, subnam)
      ncstat = nf_put_att_text(monncid, vmminat_id, 'units', 3, 'm/s')
      if ( ncstat.ne.NF_NOERR ) call handle_err (ncstat, subnam)
      ncstat = nf_put_att_text(monncid, vmminat_id, 'long_name',
     &            35, 'Minimum v in atmosphere mixed layer')
      if ( ncstat.ne.NF_NOERR ) call handle_err (ncstat, subnam)
      ncstat=nf_def_var(monncid,'vmmaxat',NF_FLOAT,1,timedim,vmmaxat_id)
      if ( ncstat.ne.NF_NOERR ) call handle_err (ncstat, subnam)
      ncstat = nf_put_att_text(monncid, vmmaxat_id, 'units', 3, 'm/s')
      if ( ncstat.ne.NF_NOERR ) call handle_err (ncstat, subnam)
      ncstat = nf_put_att_text(monncid, vmmaxat_id, 'long_name',
     &            35, 'Maximum v in atmosphere mixed layer')
      if ( ncstat.ne.NF_NOERR ) call handle_err (ncstat, subnam)

!! Define two-dimensional variables called 'ugminat' and 'ugmaxat'
      edims(1) = ladim
      edims(2) = timedim
      ncstat = nf_def_var(monncid,'ugminat',NF_FLOAT,2,edims,ugminat_id)
      if ( ncstat.ne.NF_NOERR ) call handle_err (ncstat, subnam)
      ncstat = nf_put_att_text(monncid, ugminat_id, 'units',3,'m/s')
      if ( ncstat.ne.NF_NOERR ) call handle_err (ncstat, subnam)
      ncstat = nf_put_att_text(monncid, ugminat_id, 'long_name',
     &            33, 'Minimum u in atmosphere Q-G layer')
      if ( ncstat.ne.NF_NOERR ) call handle_err (ncstat, subnam)
      edims(1) = ladim
      edims(2) = timedim
      ncstat = nf_def_var(monncid,'ugmaxat',NF_FLOAT,2,edims,ugmaxat_id)
      if ( ncstat.ne.NF_NOERR ) call handle_err (ncstat, subnam)
      ncstat = nf_put_att_text(monncid, ugmaxat_id, 'units',3,'m/s')
      if ( ncstat.ne.NF_NOERR ) call handle_err (ncstat, subnam)
      ncstat = nf_put_att_text(monncid, ugmaxat_id, 'long_name',
     &            33, 'Maximum u in atmosphere Q-G layer')
      if ( ncstat.ne.NF_NOERR ) call handle_err (ncstat, subnam)

!! Define two-dimensional variables called 'vgminat' and 'vgmaxat'
      edims(1) = ladim
      edims(2) = timedim
      ncstat = nf_def_var(monncid,'vgminat',NF_FLOAT,2,edims,vgminat_id)
      if ( ncstat.ne.NF_NOERR ) call handle_err (ncstat, subnam)
      ncstat = nf_put_att_text(monncid, vgminat_id, 'units',3,'m/s')
      if ( ncstat.ne.NF_NOERR ) call handle_err (ncstat, subnam)
      ncstat = nf_put_att_text(monncid, vgminat_id, 'long_name',
     &            33, 'Minimum v in atmosphere Q-G layer')
      if ( ncstat.ne.NF_NOERR ) call handle_err (ncstat, subnam)
      edims(1) = ladim
      edims(2) = timedim
      ncstat = nf_def_var(monncid,'vgmaxat',NF_FLOAT,2,edims,vgmaxat_id)
      if ( ncstat.ne.NF_NOERR ) call handle_err (ncstat, subnam)
      ncstat = nf_put_att_text(monncid, vgmaxat_id, 'units',3,'m/s')
      if ( ncstat.ne.NF_NOERR ) call handle_err (ncstat, subnam)
      ncstat = nf_put_att_text(monncid, vgmaxat_id, 'long_name',
     &            33, 'Maximum v in atmosphere Q-G layer')
      if ( ncstat.ne.NF_NOERR ) call handle_err (ncstat, subnam)

!! Define a one-dimensional variable called 'cnmlat'
      ncstat=nf_def_var(monncid,'cnmlat',NF_FLOAT,1,timedim,cnmlat_id)
      if ( ncstat.ne.NF_NOERR ) call handle_err (ncstat, subnam)
      ncstat = nf_put_att_text(monncid, cnmlat_id, 'units', 6, 'nondim')
      if ( ncstat.ne.NF_NOERR ) call handle_err (ncstat, subnam)
      ncstat = nf_put_att_text(monncid, cnmlat_id, 'long_name', 45,
     &           'Max. Courant number in atmosphere mixed layer')
      if ( ncstat.ne.NF_NOERR ) call handle_err (ncstat, subnam)

!! Define a two-dimensional variable called 'cnqgat'
      edims(1) = ladim
      edims(2) = timedim
      ncstat = nf_def_var(monncid,'cnqgat',NF_FLOAT,2,edims,cnqgat_id)
      if ( ncstat.ne.NF_NOERR ) call handle_err (ncstat, subnam)
      ncstat = nf_put_att_text(monncid, cnqgat_id, 'units',6,'nondim')
      if ( ncstat.ne.NF_NOERR ) call handle_err (ncstat, subnam)
      ncstat = nf_put_att_text(monncid, cnqgat_id, 'long_name', 43,
     &           'Max. Courant number in atmosphere Q-G layer')
      if ( ncstat.ne.NF_NOERR ) call handle_err (ncstat, subnam)
#  endif /* not ocean_only */

#  ifndef atmos_only
*     Ocean diagnostics from monnc_comp
*     ---------------------------------
!! Define a one-dimensional variable called 'wetmoc'
      ncstat=nf_def_var(monncid,'wetmoc',NF_FLOAT,1,timedim,wetmoc_id)
      if ( ncstat.ne.NF_NOERR ) call handle_err (ncstat, subnam)
      ncstat = nf_put_att_text(monncid, wetmoc_id, 'units', 3, 'm/s')
      if ( ncstat.ne.NF_NOERR ) call handle_err (ncstat, subnam)
      ncstat = nf_put_att_text(monncid, wetmoc_id, 'long_name',
     &            39, 'Average oceanic Ekman velocity (T-grid)')
      if ( ncstat.ne.NF_NOERR ) call handle_err (ncstat, subnam)

!! Define a one-dimensional variable called 'watmoc'
      ncstat=nf_def_var(monncid,'watmoc',NF_FLOAT,1,timedim,watmoc_id)
      if ( ncstat.ne.NF_NOERR ) call handle_err (ncstat, subnam)
      ncstat = nf_put_att_text(monncid, watmoc_id, 'units', 3, 'm/s')
      if ( ncstat.ne.NF_NOERR ) call handle_err (ncstat, subnam)
      ncstat = nf_put_att_text(monncid, watmoc_id, 'long_name',
     &            40, 'Absolute oceanic Ekman velocity (T-grid)')
      if ( ncstat.ne.NF_NOERR ) call handle_err (ncstat, subnam)

!! Define a one-dimensional variable called 'wepmoc'
      ncstat=nf_def_var(monncid,'wepmoc',NF_FLOAT,1,timedim,wepmoc_id)
      if ( ncstat.ne.NF_NOERR ) call handle_err (ncstat, subnam)
      ncstat = nf_put_att_text(monncid, wepmoc_id, 'units', 3, 'm/s')
      if ( ncstat.ne.NF_NOERR ) call handle_err (ncstat, subnam)
      ncstat = nf_put_att_text(monncid, wepmoc_id, 'long_name',
     &            39, 'Average oceanic Ekman velocity (p-grid)')
      if ( ncstat.ne.NF_NOERR ) call handle_err (ncstat, subnam)

!! Define a one-dimensional variable called 'wapmoc'
      ncstat=nf_def_var(monncid,'wapmoc',NF_FLOAT,1,timedim,wapmoc_id)
      if ( ncstat.ne.NF_NOERR ) call handle_err (ncstat, subnam)
      ncstat = nf_put_att_text(monncid, wapmoc_id, 'units', 3, 'm/s')
      if ( ncstat.ne.NF_NOERR ) call handle_err (ncstat, subnam)
      ncstat = nf_put_att_text(monncid, wapmoc_id, 'long_name',
     &            40, 'Absolute oceanic Ekman velocity (p-grid)')
      if ( ncstat.ne.NF_NOERR ) call handle_err (ncstat, subnam)

!! Define a one-dimensional variable called 'entmoc'
      ncstat=nf_def_var(monncid,'entmoc',NF_FLOAT,1,timedim,entmoc_id)
      if ( ncstat.ne.NF_NOERR ) call handle_err (ncstat, subnam)
      ncstat = nf_put_att_text(monncid, entmoc_id, 'units', 3, 'm/s')
      if ( ncstat.ne.NF_NOERR ) call handle_err (ncstat, subnam)
      ncstat = nf_put_att_text(monncid, entmoc_id, 'long_name',
     &            27, 'Average oceanic entrainment')
      if ( ncstat.ne.NF_NOERR ) call handle_err (ncstat, subnam)

!! Define a one-dimensional variable called 'enamoc'
      ncstat=nf_def_var(monncid,'enamoc',NF_FLOAT,1,timedim,enamoc_id)
      if ( ncstat.ne.NF_NOERR ) call handle_err (ncstat, subnam)
      ncstat = nf_put_att_text(monncid, enamoc_id, 'units', 3, 'm/s')
      if ( ncstat.ne.NF_NOERR ) call handle_err (ncstat, subnam)
      ncstat = nf_put_att_text(monncid, enamoc_id, 'long_name',
     &            28, 'Absolute oceanic entrainment')
      if ( ncstat.ne.NF_NOERR ) call handle_err (ncstat, subnam)

!! Define a two-dimensional variable called 'pkenoc'
      ncstat=nf_def_var(monncid,'pkenoc',NF_FLOAT,1,timedim,pkenoc_id)
      if ( ncstat.ne.NF_NOERR ) call handle_err (ncstat, subnam)
      ncstat = nf_put_att_text(monncid, pkenoc_id, 'units',5,'W/m^2')
      if ( ncstat.ne.NF_NOERR ) call handle_err (ncstat, subnam)
      ncstat = nf_put_att_text(monncid, pkenoc_id, 'long_name',
     &            21, 'Average oceanic eta*e')
      if ( ncstat.ne.NF_NOERR ) call handle_err (ncstat, subnam)

!! Define a two-dimensional variable called 'etamoc'
      edims(1) = lomdim
      edims(2) = timedim
      ncstat = nf_def_var(monncid,'etamoc',NF_FLOAT,2,edims,etamoc_id)
      if ( ncstat.ne.NF_NOERR ) call handle_err (ncstat, subnam)
      ncstat = nf_put_att_text(monncid, etamoc_id, 'units', 1, 'm')
      if ( ncstat.ne.NF_NOERR ) call handle_err (ncstat, subnam)
      ncstat = nf_put_att_text(monncid, etamoc_id, 'long_name',
     &            32, 'Average oceanic interface height')
      if ( ncstat.ne.NF_NOERR ) call handle_err (ncstat, subnam)

!! Define a two-dimensional variable called 'ddtkeoc'
      edims(1) = lodim
      edims(2) = timedim
      ncstat = nf_def_var(monncid,'ddtkeoc',NF_FLOAT,2,edims,ddtkeoc_id)
      if ( ncstat.ne.NF_NOERR ) call handle_err (ncstat, subnam)
      ncstat = nf_put_att_text(monncid, ddtkeoc_id, 'units',5,'W/m^2')
      if ( ncstat.ne.NF_NOERR ) call handle_err (ncstat, subnam)
      ncstat = nf_put_att_text(monncid, ddtkeoc_id, 'long_name',
     &            31, 'd/dt(average KE in ocean layer)')
      if ( ncstat.ne.NF_NOERR ) call handle_err (ncstat, subnam)

!! Define a two-dimensional variable called 'ddtpeoc'
      edims(1) = lomdim
      edims(2) = timedim
      ncstat = nf_def_var(monncid,'ddtpeoc',NF_FLOAT,2,edims,ddtpeoc_id)
      if ( ncstat.ne.NF_NOERR ) call handle_err (ncstat, subnam)
      ncstat = nf_put_att_text(monncid, ddtpeoc_id, 'units',5,'W/m^2')
      if ( ncstat.ne.NF_NOERR ) call handle_err (ncstat, subnam)
      ncstat = nf_put_att_text(monncid, ddtpeoc_id, 'long_name',
     &            35, 'd/dt(average PE at ocean interface)')
      if ( ncstat.ne.NF_NOERR ) call handle_err (ncstat, subnam)

!! Define a two-dimensional variable called 'et2moc'
      edims(1) = lomdim
      edims(2) = timedim
      ncstat = nf_def_var(monncid,'et2moc',NF_FLOAT,2,edims,et2moc_id)
      if ( ncstat.ne.NF_NOERR ) call handle_err (ncstat, subnam)
      ncstat = nf_put_att_text(monncid, et2moc_id, 'units',5,'W/m^2')
      if ( ncstat.ne.NF_NOERR ) call handle_err (ncstat, subnam)
      ncstat = nf_put_att_text(monncid, et2moc_id, 'long_name',
     &            21, 'Average oceanic eta^2')
      if ( ncstat.ne.NF_NOERR ) call handle_err (ncstat, subnam)

!! Define a one-dimensional variable called 'utauoc'
      ncstat=nf_def_var(monncid,'utauoc',NF_FLOAT,1,timedim,utauoc_id)
      if ( ncstat.ne.NF_NOERR ) call handle_err (ncstat, subnam)
      ncstat = nf_put_att_text(monncid, utauoc_id, 'units',5,'W/m^2')
      if ( ncstat.ne.NF_NOERR ) call handle_err (ncstat, subnam)
      ncstat = nf_put_att_text(monncid, utauoc_id, 'long_name',
     &            21, 'Average oceanic p1*we')
      if ( ncstat.ne.NF_NOERR ) call handle_err (ncstat, subnam)

!! Define a two-dimensional variable called 'ocjpos'
      edims(1) = lodim
      edims(2) = timedim
      ncstat = nf_def_var(monncid,'ocjpos',NF_INT,2,edims,ocjpos_id)
      if ( ncstat.ne.NF_NOERR ) call handle_err (ncstat, subnam)
      ncstat=nf_put_att_text(monncid,ocjpos_id,'units',10,'gridsquare')
      if ( ncstat.ne.NF_NOERR ) call handle_err (ncstat, subnam)
      ncstat = nf_put_att_text(monncid,ocjpos_id,'long_name',51,
     &        'Position of maximum ocean jet velocity (gridsquare)')
      if ( ncstat.ne.NF_NOERR ) call handle_err (ncstat, subnam)

!! Define a two-dimensional variable called 'ocjval'
      edims(1) = lodim
      edims(2) = timedim
      ncstat = nf_def_var(monncid,'ocjval',NF_FLOAT,2,edims,ocjval_id)
      if ( ncstat.ne.NF_NOERR ) call handle_err (ncstat, subnam)
      ncstat = nf_put_att_text(monncid, ocjval_id, 'units', 3, 'm/s')
      if ( ncstat.ne.NF_NOERR ) call handle_err (ncstat, subnam)
      ncstat = nf_put_att_text(monncid, ocjval_id, 'long_name',
     &            38, 'Maximum ocean zonal jet velocity (m/s)')
      if ( ncstat.ne.NF_NOERR ) call handle_err (ncstat, subnam)

!! Define a two-dimensional variable called 'pavgoc'
      edims(1) = lodim
      edims(2) = timedim
      ncstat = nf_def_var(monncid,'pavgoc',NF_FLOAT,2,edims,pavgoc_id)
      if ( ncstat.ne.NF_NOERR ) call handle_err (ncstat, subnam)
      ncstat = nf_put_att_text(monncid,pavgoc_id,'units',8,'m^2 s^-2')
      if ( ncstat.ne.NF_NOERR ) call handle_err (ncstat, subnam)
      ncstat = nf_put_att_text(monncid,pavgoc_id,'long_name',
     &            32, 'Average oceanic dynamic pressure')
      if ( ncstat.ne.NF_NOERR ) call handle_err (ncstat, subnam)

!! Define a two-dimensional variable called 'qavgoc'
      edims(1) = lodim
      edims(2) = timedim
      ncstat = nf_def_var(monncid,'qavgoc',NF_FLOAT,2,edims,qavgoc_id)
      if ( ncstat.ne.NF_NOERR ) call handle_err (ncstat, subnam)
      ncstat = nf_put_att_text(monncid, qavgoc_id, 'units', 4, 's^-1')
      if ( ncstat.ne.NF_NOERR ) call handle_err (ncstat, subnam)
      ncstat = nf_put_att_text(monncid, qavgoc_id, 'long_name',
     &            25, 'Average oceanic vorticity')
      if ( ncstat.ne.NF_NOERR ) call handle_err (ncstat, subnam)

!! Define a two-dimensional variable called 'ah2doc'
      edims(1) = lodim
      edims(2) = timedim
      ncstat = nf_def_var(monncid,'ah2doc',NF_FLOAT,2,edims,ah2doc_id)
      if ( ncstat.ne.NF_NOERR ) call handle_err (ncstat, subnam)
      ncstat = nf_put_att_text(monncid, ah2doc_id, 'units',5,'W/m^2')
      if ( ncstat.ne.NF_NOERR ) call handle_err (ncstat, subnam)
      ncstat = nf_put_att_text(monncid, ah2doc_id, 'long_name',
     &            39, 'Average oceanic dissipation (2nd order)')
      if ( ncstat.ne.NF_NOERR ) call handle_err (ncstat, subnam)

!! Define a two-dimensional variable called 'ah4doc'
      edims(1) = lodim
      edims(2) = timedim
      ncstat = nf_def_var(monncid,'ah4doc',NF_FLOAT,2,edims,ah4doc_id)
      if ( ncstat.ne.NF_NOERR ) call handle_err (ncstat, subnam)
      ncstat = nf_put_att_text(monncid, ah4doc_id, 'units',5,'W/m^2')
      if ( ncstat.ne.NF_NOERR ) call handle_err (ncstat, subnam)
      ncstat = nf_put_att_text(monncid, ah4doc_id, 'long_name',
     &            39, 'Average oceanic dissipation (4th order)')
      if ( ncstat.ne.NF_NOERR ) call handle_err (ncstat, subnam)

!! Define a two-dimensional variable called 'kealoc'
      edims(1) = lodim
      edims(2) = timedim
      ncstat = nf_def_var(monncid,'kealoc',NF_FLOAT,2,edims,kealoc_id)
      if ( ncstat.ne.NF_NOERR ) call handle_err (ncstat, subnam)
      ncstat = nf_put_att_text(monncid, kealoc_id, 'units',5,'J/m^2')
      if ( ncstat.ne.NF_NOERR ) call handle_err (ncstat, subnam)
      ncstat = nf_put_att_text(monncid, kealoc_id, 'long_name',
     &            36, 'Average oceanic layer kinetic energy')
      if ( ncstat.ne.NF_NOERR ) call handle_err (ncstat, subnam)

!! Define a two-dimensional variable called 'osfmin'
      edims(1) = lodim
      edims(2) = timedim
      ncstat = nf_def_var(monncid,'osfmin',NF_FLOAT,2,edims,osfmin_id)
      if ( ncstat.ne.NF_NOERR ) call handle_err (ncstat, subnam)
      ncstat = nf_put_att_text(monncid, osfmin_id, 'units', 2, 'Sv')
      if ( ncstat.ne.NF_NOERR ) call handle_err (ncstat, subnam)
      ncstat = nf_put_att_text(monncid, osfmin_id, 'long_name',
     &            22, 'Minimum streamfunction')
      if ( ncstat.ne.NF_NOERR ) call handle_err (ncstat, subnam)

!! Define a two-dimensional variable called 'osfmax'
      edims(1) = lodim
      edims(2) = timedim
      ncstat = nf_def_var(monncid,'osfmax',NF_FLOAT,2,edims,osfmax_id)
      if ( ncstat.ne.NF_NOERR ) call handle_err (ncstat, subnam)
      ncstat = nf_put_att_text(monncid, osfmax_id, 'units', 2, 'Sv')
      if ( ncstat.ne.NF_NOERR ) call handle_err (ncstat, subnam)
      ncstat = nf_put_att_text(monncid, osfmax_id, 'long_name',
     &            22, 'Maximum streamfunction')
      if ( ncstat.ne.NF_NOERR ) call handle_err (ncstat, subnam)

!! Define a two-dimensional variable called 'occirc'
      edims(1) = lodim
      edims(2) = timedim
      ncstat = nf_def_var(monncid,'occirc',NF_FLOAT,2,edims,occirc_id)
      if ( ncstat.ne.NF_NOERR ) call handle_err (ncstat, subnam)
      ncstat = nf_put_att_text(monncid, occirc_id, 'units', 2, 'Sv')
      if ( ncstat.ne.NF_NOERR ) call handle_err (ncstat, subnam)
      ncstat = nf_put_att_text(monncid, occirc_id, 'long_name',
     &            42, 'Geostrophic zonal transport in ocean layer')
      if ( ncstat.ne.NF_NOERR ) call handle_err (ncstat, subnam)

!! Define a one-dimensional variable called 'btdgoc'
      ncstat=nf_def_var(monncid,'btdgoc',NF_FLOAT,1,timedim,btdgoc_id)
      if ( ncstat.ne.NF_NOERR ) call handle_err (ncstat, subnam)
      ncstat = nf_put_att_text(monncid, btdgoc_id, 'units',5,'W/m^2')
      if ( ncstat.ne.NF_NOERR ) call handle_err (ncstat, subnam)
      ncstat = nf_put_att_text(monncid, btdgoc_id, 'long_name',
     &            39, 'Average oceanic bottom drag dissipation')
      if ( ncstat.ne.NF_NOERR ) call handle_err (ncstat, subnam)

!! Define a one-dimensional variable called 'tmlmoc'
      ncstat=nf_def_var(monncid,'tmlmoc',NF_FLOAT,1,timedim,tmlmoc_id)
      if ( ncstat.ne.NF_NOERR ) call handle_err (ncstat, subnam)
      ncstat = nf_put_att_text(monncid, tmlmoc_id, 'units', 1, 'K')
      if ( ncstat.ne.NF_NOERR ) call handle_err (ncstat, subnam)
      ncstat = nf_put_att_text(monncid, tmlmoc_id, 'long_name',
     &            39, 'Average oceanic mixed layer temperature')
      if ( ncstat.ne.NF_NOERR ) call handle_err (ncstat, subnam)

!! Define one-dimensional variables called 'sstmin' and 'sstmax'
      ncstat=nf_def_var(monncid,'sstmin',NF_FLOAT,1,timedim,sstmin_id)
      if ( ncstat.ne.NF_NOERR ) call handle_err (ncstat, subnam)
      ncstat = nf_put_att_text(monncid, sstmin_id, 'units', 1, 'K')
      if ( ncstat.ne.NF_NOERR ) call handle_err (ncstat, subnam)
      ncstat = nf_put_att_text(monncid, sstmin_id, 'long_name',
     &            36, 'Min. oceanic mixed layer temperature')
      if ( ncstat.ne.NF_NOERR ) call handle_err (ncstat, subnam)
      ncstat=nf_def_var(monncid,'sstmax',NF_FLOAT,1,timedim,sstmax_id)
      if ( ncstat.ne.NF_NOERR ) call handle_err (ncstat, subnam)
      ncstat = nf_put_att_text(monncid, sstmax_id, 'units', 1, 'K')
      if ( ncstat.ne.NF_NOERR ) call handle_err (ncstat, subnam)
      ncstat = nf_put_att_text(monncid, sstmax_id, 'long_name',
     &            36, 'Max. oceanic mixed layer temperature')
      if ( ncstat.ne.NF_NOERR ) call handle_err (ncstat, subnam)

!! Define a one-dimensional variable called 'hfmloc'
      ncstat=nf_def_var(monncid,'hfmloc',NF_FLOAT,1,timedim,hfmloc_id)
      if ( ncstat.ne.NF_NOERR ) call handle_err (ncstat, subnam)
      ncstat = nf_put_att_text(monncid, hfmloc_id, 'units',5,'W/m^2')
      if ( ncstat.ne.NF_NOERR ) call handle_err (ncstat, subnam)
      ncstat = nf_put_att_text(monncid, hfmloc_id, 'long_name',
     &            40, 'Heat flux at bottom of ocean mixed layer')
      if ( ncstat.ne.NF_NOERR ) call handle_err (ncstat, subnam)

!! Define a one-dimensional variable called 'occtot'
      ncstat = nf_def_var(monncid,'occtot',NF_FLOAT,1,timedim,occtot_id)
      if ( ncstat.ne.NF_NOERR ) call handle_err (ncstat, subnam)
      ncstat = nf_put_att_text(monncid, occtot_id, 'units', 2, 'Sv')
      if ( ncstat.ne.NF_NOERR ) call handle_err (ncstat, subnam)
      ncstat = nf_put_att_text(monncid, occtot_id, 'long_name',
     &            27, 'Total zonal ocean transport')
      if ( ncstat.ne.NF_NOERR ) call handle_err (ncstat, subnam)

!! Define one-dimensional variables called 'umminoc' and 'ummaxoc'
      ncstat=nf_def_var(monncid,'umminoc',NF_FLOAT,1,timedim,umminoc_id)
      if ( ncstat.ne.NF_NOERR ) call handle_err (ncstat, subnam)
      ncstat = nf_put_att_text(monncid, umminoc_id, 'units', 3, 'm/s')
      if ( ncstat.ne.NF_NOERR ) call handle_err (ncstat, subnam)
      ncstat = nf_put_att_text(monncid, umminoc_id, 'long_name',
     &            30, 'Minimum u in ocean mixed layer')
      if ( ncstat.ne.NF_NOERR ) call handle_err (ncstat, subnam)
      ncstat=nf_def_var(monncid,'ummaxoc',NF_FLOAT,1,timedim,ummaxoc_id)
      if ( ncstat.ne.NF_NOERR ) call handle_err (ncstat, subnam)
      ncstat = nf_put_att_text(monncid, ummaxoc_id, 'units', 3, 'm/s')
      if ( ncstat.ne.NF_NOERR ) call handle_err (ncstat, subnam)
      ncstat = nf_put_att_text(monncid, ummaxoc_id, 'long_name',
     &            30, 'Maximum u in ocean mixed layer')
      if ( ncstat.ne.NF_NOERR ) call handle_err (ncstat, subnam)

!! Define one-dimensional variables called 'vmminoc' and 'vmmaxoc'
      ncstat=nf_def_var(monncid,'vmminoc',NF_FLOAT,1,timedim,vmminoc_id)
      if ( ncstat.ne.NF_NOERR ) call handle_err (ncstat, subnam)
      ncstat = nf_put_att_text(monncid, vmminoc_id, 'units', 3, 'm/s')
      if ( ncstat.ne.NF_NOERR ) call handle_err (ncstat, subnam)
      ncstat = nf_put_att_text(monncid, vmminoc_id, 'long_name',
     &            30, 'Minimum v in ocean mixed layer')
      if ( ncstat.ne.NF_NOERR ) call handle_err (ncstat, subnam)
      ncstat=nf_def_var(monncid,'vmmaxoc',NF_FLOAT,1,timedim,vmmaxoc_id)
      if ( ncstat.ne.NF_NOERR ) call handle_err (ncstat, subnam)
      ncstat = nf_put_att_text(monncid, vmmaxoc_id, 'units', 3, 'm/s')
      if ( ncstat.ne.NF_NOERR ) call handle_err (ncstat, subnam)
      ncstat = nf_put_att_text(monncid, vmmaxoc_id, 'long_name',
     &            30, 'Maximum v in ocean mixed layer')
      if ( ncstat.ne.NF_NOERR ) call handle_err (ncstat, subnam)

!! Define two-dimensional variables called 'ugminoc' and 'ugmaxoc'
      edims(1) = lodim
      edims(2) = timedim
      ncstat = nf_def_var(monncid,'ugminoc',NF_FLOAT,2,edims,ugminoc_id)
      if ( ncstat.ne.NF_NOERR ) call handle_err (ncstat, subnam)
      ncstat = nf_put_att_text(monncid, ugminoc_id, 'units',3,'m/s')
      if ( ncstat.ne.NF_NOERR ) call handle_err (ncstat, subnam)
      ncstat = nf_put_att_text(monncid, ugminoc_id, 'long_name',
     &            28, 'Minimum u in ocean Q-G layer')
      if ( ncstat.ne.NF_NOERR ) call handle_err (ncstat, subnam)
      edims(1) = lodim
      edims(2) = timedim
      ncstat = nf_def_var(monncid,'ugmaxoc',NF_FLOAT,2,edims,ugmaxoc_id)
      if ( ncstat.ne.NF_NOERR ) call handle_err (ncstat, subnam)
      ncstat = nf_put_att_text(monncid, ugmaxoc_id, 'units',3,'m/s')
      if ( ncstat.ne.NF_NOERR ) call handle_err (ncstat, subnam)
      ncstat = nf_put_att_text(monncid, ugmaxoc_id, 'long_name',
     &            28, 'Maximum u in ocean Q-G layer')
      if ( ncstat.ne.NF_NOERR ) call handle_err (ncstat, subnam)

!! Define two-dimensional variables called 'vgminoc' and 'vgmaxoc'
      edims(1) = lodim
      edims(2) = timedim
      ncstat = nf_def_var(monncid,'vgminoc',NF_FLOAT,2,edims,vgminoc_id)
      if ( ncstat.ne.NF_NOERR ) call handle_err (ncstat, subnam)
      ncstat = nf_put_att_text(monncid, vgminoc_id, 'units',3,'m/s')
      if ( ncstat.ne.NF_NOERR ) call handle_err (ncstat, subnam)
      ncstat = nf_put_att_text(monncid, vgminoc_id, 'long_name',
     &            28, 'Minimum v in ocean Q-G layer')
      if ( ncstat.ne.NF_NOERR ) call handle_err (ncstat, subnam)
      edims(1) = lodim
      edims(2) = timedim
      ncstat = nf_def_var(monncid,'vgmaxoc',NF_FLOAT,2,edims,vgmaxoc_id)
      if ( ncstat.ne.NF_NOERR ) call handle_err (ncstat, subnam)
      ncstat = nf_put_att_text(monncid, vgmaxoc_id, 'units',3,'m/s')
      if ( ncstat.ne.NF_NOERR ) call handle_err (ncstat, subnam)
      ncstat = nf_put_att_text(monncid, vgmaxoc_id, 'long_name',
     &            28, 'Maximum v in ocean Q-G layer')
      if ( ncstat.ne.NF_NOERR ) call handle_err (ncstat, subnam)

!! Define a one-dimensional variable called 'cnmloc'
      ncstat=nf_def_var(monncid,'cnmloc',NF_FLOAT,1,timedim,cnmloc_id)
      if ( ncstat.ne.NF_NOERR ) call handle_err (ncstat, subnam)
      ncstat = nf_put_att_text(monncid, cnmloc_id, 'units', 6, 'nondim')
      if ( ncstat.ne.NF_NOERR ) call handle_err (ncstat, subnam)
      ncstat = nf_put_att_text(monncid, cnmloc_id, 'long_name', 40,
     &           'Max. Courant number in ocean mixed layer')
      if ( ncstat.ne.NF_NOERR ) call handle_err (ncstat, subnam)

!! Define a two-dimensional variable called 'cnqgoc'
      edims(1) = lodim
      edims(2) = timedim
      ncstat = nf_def_var(monncid,'cnqgoc',NF_FLOAT,2,edims,cnqgoc_id)
      if ( ncstat.ne.NF_NOERR ) call handle_err (ncstat, subnam)
      ncstat = nf_put_att_text(monncid, cnqgoc_id, 'units',6,'nondim')
      if ( ncstat.ne.NF_NOERR ) call handle_err (ncstat, subnam)
      ncstat = nf_put_att_text(monncid, cnqgoc_id, 'long_name', 38,
     &           'Max. Courant number in ocean Q-G layer')
      if ( ncstat.ne.NF_NOERR ) call handle_err (ncstat, subnam)
#  endif /* not atmos_only */

#  ifndef atmos_only
*     Ocean diagnostics from oml
*     --------------------------
!! Define a one-dimensional variable called 'cfraoc'
      ncstat=nf_def_var(monncid,'cfraoc',NF_FLOAT,1,timedim,cfraoc_id)
      if ( ncstat.ne.NF_NOERR ) call handle_err (ncstat, subnam)
      ncstat = nf_put_att_text(monncid, cfraoc_id, 'units', 1, ' ')
      if ( ncstat.ne.NF_NOERR ) call handle_err (ncstat, subnam)
      ncstat = nf_put_att_text(monncid, cfraoc_id, 'long_name',
     &            28, 'Fraction of ocean convecting')
      if ( ncstat.ne.NF_NOERR ) call handle_err (ncstat, subnam)

!! Define a one-dimensional variable called 'centoc'
      ncstat=nf_def_var(monncid,'centoc',NF_FLOAT,1,timedim,centoc_id)
      if ( ncstat.ne.NF_NOERR ) call handle_err (ncstat, subnam)
      ncstat = nf_put_att_text(monncid, centoc_id, 'units', 3, 'm/s')
      if ( ncstat.ne.NF_NOERR ) call handle_err (ncstat, subnam)
      ncstat = nf_put_att_text(monncid, centoc_id, 'long_name',
     &            44, 'Sum of oceanic entrainment due to convection')
      if ( ncstat.ne.NF_NOERR ) call handle_err (ncstat, subnam)
#  endif /* not atmos_only */

#  ifndef atmos_only
*     Ocean m. l. boundary diagnostics from omladf
*     --------------------------------------------
!! Define a one-dimensional variable called 'ttmads'
      ncstat=nf_def_var(monncid,'ttmads',NF_FLOAT,1,timedim,ttmads_id)
      if ( ncstat.ne.NF_NOERR ) call handle_err (ncstat, subnam)
      ncstat = nf_put_att_text(monncid, ttmads_id, 'units', 3, 'K/s')
      if ( ncstat.ne.NF_NOERR ) call handle_err (ncstat, subnam)
      ncstat = nf_put_att_text(monncid, ttmads_id, 'long_name', 48,
     &         'Average advective heat flux at Southern boundary')
      if ( ncstat.ne.NF_NOERR ) call handle_err (ncstat, subnam)

!! Define a one-dimensional variable called 'vfmads'
      ncstat=nf_def_var(monncid,'vfmads',NF_FLOAT,1,timedim,vfmads_id)
      if ( ncstat.ne.NF_NOERR ) call handle_err (ncstat, subnam)
      ncstat = nf_put_att_text(monncid, vfmads_id, 'units', 3, 'm/s')
      if ( ncstat.ne.NF_NOERR ) call handle_err (ncstat, subnam)
      ncstat = nf_put_att_text(monncid, vfmads_id, 'long_name',
     &            43, 'Average velocity at Southern ocean boundary')
      if ( ncstat.ne.NF_NOERR ) call handle_err (ncstat, subnam)

!! Define a one-dimensional variable called 'ttmdfs'
      ncstat=nf_def_var(monncid,'ttmdfs',NF_FLOAT,1,timedim,ttmdfs_id)
      if ( ncstat.ne.NF_NOERR ) call handle_err (ncstat, subnam)
      ncstat = nf_put_att_text(monncid, ttmdfs_id, 'units', 3, 'K/s')
      if ( ncstat.ne.NF_NOERR ) call handle_err (ncstat, subnam)
      ncstat = nf_put_att_text(monncid, ttmdfs_id, 'long_name', 48,
     &         'Average diffusive heat flux at Southern boundary')
      if ( ncstat.ne.NF_NOERR ) call handle_err (ncstat, subnam)

!! Define a one-dimensional variable called 'ttmadn'
      ncstat=nf_def_var(monncid,'ttmadn',NF_FLOAT,1,timedim,ttmadn_id)
      if ( ncstat.ne.NF_NOERR ) call handle_err (ncstat, subnam)
      ncstat = nf_put_att_text(monncid, ttmadn_id, 'units', 3, 'K/s')
      if ( ncstat.ne.NF_NOERR ) call handle_err (ncstat, subnam)
      ncstat = nf_put_att_text(monncid, ttmadn_id, 'long_name', 48,
     &         'Average advective heat flux at Northern boundary')
      if ( ncstat.ne.NF_NOERR ) call handle_err (ncstat, subnam)

!! Define a one-dimensional variable called 'vfmadn'
      ncstat=nf_def_var(monncid,'vfmadn',NF_FLOAT,1,timedim,vfmadn_id)
      if ( ncstat.ne.NF_NOERR ) call handle_err (ncstat, subnam)
      ncstat = nf_put_att_text(monncid, vfmadn_id, 'units', 3, 'm/s')
      if ( ncstat.ne.NF_NOERR ) call handle_err (ncstat, subnam)
      ncstat = nf_put_att_text(monncid, vfmadn_id, 'long_name',
     &            43, 'Average velocity at Northern ocean boundary')
      if ( ncstat.ne.NF_NOERR ) call handle_err (ncstat, subnam)

!! Define a one-dimensional variable called 'ttmdfn'
      ncstat=nf_def_var(monncid,'ttmdfn',NF_FLOAT,1,timedim,ttmdfn_id)
      if ( ncstat.ne.NF_NOERR ) call handle_err (ncstat, subnam)
      ncstat = nf_put_att_text(monncid, ttmdfn_id, 'units', 3, 'K/s')
      if ( ncstat.ne.NF_NOERR ) call handle_err (ncstat, subnam)
      ncstat = nf_put_att_text(monncid, ttmdfn_id, 'long_name', 48,
     &         'Average diffusive heat flux at Northern boundary')
      if ( ncstat.ne.NF_NOERR ) call handle_err (ncstat, subnam)
#  endif /* not atmos_only */

#  ifndef ocean_only
*     Flux diagnostics from xforc
*     ---------------------------
!! Define a one-dimensional variable called 'slhfav'
      ncstat=nf_def_var(monncid,'slhfav',NF_FLOAT,1,timedim,slhfav_id)
      if ( ncstat.ne.NF_NOERR ) call handle_err (ncstat, subnam)
      ncstat = nf_put_att_text(monncid, slhfav_id, 'units',5,'W/m^2')
      if ( ncstat.ne.NF_NOERR ) call handle_err (ncstat, subnam)
      ncstat = nf_put_att_text(monncid, slhfav_id, 'long_name', 45,
     &         'Average oceanic sensible and latent heat flux')
      if ( ncstat.ne.NF_NOERR ) call handle_err (ncstat, subnam)

!! Define a one-dimensional variable called 'oradav'
      ncstat=nf_def_var(monncid,'oradav',NF_FLOAT,1,timedim,oradav_id)
      if ( ncstat.ne.NF_NOERR ) call handle_err (ncstat, subnam)
      ncstat = nf_put_att_text(monncid, oradav_id, 'units',5,'W/m^2')
      if ( ncstat.ne.NF_NOERR ) call handle_err (ncstat, subnam)
      ncstat = nf_put_att_text(monncid, oradav_id, 'long_name',
     &               35, 'Average oceanic radiative heat flux')
      if ( ncstat.ne.NF_NOERR ) call handle_err (ncstat, subnam)

!! Define a one-dimensional variable called 'arocav'
      ncstat=nf_def_var(monncid,'arocav',NF_FLOAT,1,timedim,arocav_id)
      if ( ncstat.ne.NF_NOERR ) call handle_err (ncstat, subnam)
      ncstat = nf_put_att_text(monncid, arocav_id, 'units',5,'W/m^2')
      if ( ncstat.ne.NF_NOERR ) call handle_err (ncstat, subnam)
      ncstat = nf_put_att_text(monncid, arocav_id, 'long_name', 52,
     &         'Average atmospheric radiative heat flux (over ocean)')
      if ( ncstat.ne.NF_NOERR ) call handle_err (ncstat, subnam)

!! Define a one-dimensional variable called 'arlaav'
      ncstat=nf_def_var(monncid,'arlaav',NF_FLOAT,1,timedim,arlaav_id)
      if ( ncstat.ne.NF_NOERR ) call handle_err (ncstat, subnam)
      ncstat = nf_put_att_text(monncid, arlaav_id, 'units',5,'W/m^2')
      if ( ncstat.ne.NF_NOERR ) call handle_err (ncstat, subnam)
      ncstat = nf_put_att_text(monncid, arlaav_id, 'long_name', 51,
     &         'Average atmospheric radiative heat flux (over land)')
      if ( ncstat.ne.NF_NOERR ) call handle_err (ncstat, subnam)
#  endif /* not ocean_only */

#  ifndef ocean_only
*     Atmosphere diagnostics from aml
*     -------------------------------
!! Define a one-dimensional variable called 'cfraat'
      ncstat=nf_def_var(monncid,'cfraat',NF_FLOAT,1,timedim,cfraat_id)
      if ( ncstat.ne.NF_NOERR ) call handle_err (ncstat, subnam)
      ncstat = nf_put_att_text(monncid, cfraat_id, 'units', 1, ' ')
      if ( ncstat.ne.NF_NOERR ) call handle_err (ncstat, subnam)
      ncstat = nf_put_att_text(monncid, cfraat_id, 'long_name',
     &            33, 'Fraction of atmosphere convecting')
      if ( ncstat.ne.NF_NOERR ) call handle_err (ncstat, subnam)

!! Define a one-dimensional variable called 'centat'
      ncstat=nf_def_var(monncid,'centat',NF_FLOAT,1,timedim,centat_id)
      if ( ncstat.ne.NF_NOERR ) call handle_err (ncstat, subnam)
      ncstat = nf_put_att_text(monncid, centat_id, 'units',3,'m/s')
      if ( ncstat.ne.NF_NOERR ) call handle_err (ncstat, subnam)
      ncstat = nf_put_att_text(monncid, centat_id, 'long_name', 48,
     &         'Sum of atmospheric entrainment due to convection')
      if ( ncstat.ne.NF_NOERR ) call handle_err (ncstat, subnam)
#  endif /* not ocean_only */

#  ifndef ocean_only
*     Consistency diagnostics from atinvq
*     -----------------------------------
!! Define a two-dimensional variable called 'ermasa'
      edims(1) = lamdim
      edims(2) = timedim
      ncstat = nf_def_var(monncid,'ermasa',NF_FLOAT,2,edims,ermasa_id)
      if ( ncstat.ne.NF_NOERR ) call handle_err (ncstat, subnam)
      ncstat = nf_put_att_text(monncid, ermasa_id, 'units',7,'m^4/s^2')
      if ( ncstat.ne.NF_NOERR ) call handle_err (ncstat, subnam)
      ncstat = nf_put_att_text(monncid, ermasa_id, 'long_name',
     &            29, 'Sum of atmospheric mass error')
      if ( ncstat.ne.NF_NOERR ) call handle_err (ncstat, subnam)

!! Define a two-dimensional variable called 'emfrat'
      edims(1) = lamdim
      edims(2) = timedim
      ncstat = nf_def_var(monncid,'emfrat',NF_FLOAT,2,edims,emfrat_id)
      if ( ncstat.ne.NF_NOERR ) call handle_err (ncstat, subnam)
      ncstat = nf_put_att_text(monncid, emfrat_id, 'units', 1, ' ')
      if ( ncstat.ne.NF_NOERR ) call handle_err (ncstat, subnam)
      ncstat = nf_put_att_text(monncid, emfrat_id, 'long_name',
     &            33, 'Fractional atmospheric mass error')
      if ( ncstat.ne.NF_NOERR ) call handle_err (ncstat, subnam)
#  endif /* not ocean_only */

#  ifndef atmos_only
*     Consistency diagnostics from ocinvq
*     -----------------------------------
#    ifdef cyclic_ocean
!! Define a two-dimensional variable called 'ermaso'
      edims(1) = lomdim
      edims(2) = timedim
      ncstat = nf_def_var(monncid,'ermaso',NF_FLOAT,2,edims,ermaso_id)
      if ( ncstat.ne.NF_NOERR ) call handle_err (ncstat, subnam)
      ncstat = nf_put_att_text(monncid, ermaso_id, 'units',7,'m^4/s^2')
      if ( ncstat.ne.NF_NOERR ) call handle_err (ncstat, subnam)
      ncstat = nf_put_att_text(monncid, ermaso_id, 'long_name',
     &            25, 'Sum of oceanic mass error')
      if ( ncstat.ne.NF_NOERR ) call handle_err (ncstat, subnam)

!! Define a two-dimensional variable called 'emfroc'
      edims(1) = lomdim
      edims(2) = timedim
      ncstat = nf_def_var(monncid,'emfroc',NF_FLOAT,2,edims,emfroc_id)
      if ( ncstat.ne.NF_NOERR ) call handle_err (ncstat, subnam)
      ncstat = nf_put_att_text(monncid, emfroc_id, 'units', 1, ' ')
      if ( ncstat.ne.NF_NOERR ) call handle_err (ncstat, subnam)
      ncstat = nf_put_att_text(monncid, emfroc_id, 'long_name',
     &            29, 'Fractional oceanic mass error')
      if ( ncstat.ne.NF_NOERR ) call handle_err (ncstat, subnam)
#    endif
#  endif /* not atmos_only */

!! Leave definition mode: entering data mode.
      ncstat = nf_enddef(monncid)
      if ( ncstat.ne.NF_NOERR ) then
        print *,' monnc_init: enddef problem for monitoring file'
        call handle_err (ncstat, subnam)
      endif

!! Convert height into km and store in 'zo'
#  ifndef atmos_only
      tomp(1) = 0.5d-3*hoc(1)
      do k=2,nlo
        tomp(k) = tomp(k-1) + 0.5d-3*( hoc(k-1) + hoc(k) )
      enddo
**    write(*,'(a,9f12.5)') '  zo (km): ',tomp
      ncstat = nf_put_var_double(monncid, lo_id, tomp)
      if ( ncstat.ne.NF_NOERR ) call handle_err (ncstat, subnam)
!! Convert height into km and store in 'zom'
      tomp(1) = 1.0d-3*hoc(1)
      do k=2,nlo-1
        tomp(k) = tomp(k-1) + 1.0d-3*hoc(k)
      enddo
**    write(*,'(a,9f12.5)') '  zom (km): ',tomp
      ncstat = nf_put_var_double(monncid, lom_id, tomp)
      if ( ncstat.ne.NF_NOERR ) call handle_err (ncstat, subnam)
#  endif /* not atmos_only */

!! Convert height into km and store in 'za'
#  ifndef ocean_only
      tamp(1) = 0.5d-3*hat(1)
      do k=2,nla
        tamp(k) = tamp(k-1) + 0.5d-3*( hat(k-1) + hat(k) )
      enddo
**    write(*,'(a,9f12.5)') '  za (km): ',tamp
      ncstat = nf_put_var_double(monncid, la_id, tamp)
      if ( ncstat.ne.NF_NOERR ) call handle_err (ncstat, subnam)
!! Convert height into km and store in 'zam'
      tamp(1) = 1.0d-3*hat(1)
      do k=2,nla-1
        tamp(k) = tamp(k-1) + 1.0d-3*hat(k)
      enddo
**    write(*,'(a,9f12.5)') '  zam (km): ',tamp
      ncstat = nf_put_var_double(monncid, lam_id, tamp)
      if ( ncstat.ne.NF_NOERR ) call handle_err (ncstat, subnam)
#  endif /* not ocean_only */

#endif

      END SUBROUTINE monnc_init

c***********************************************************************

      SUBROUTINE monnc_out

*     Modules
      USE parameters, ONLY : nla, nlo
      USE timinfo, ONLY : ntdone, nocmon, tyrs
      USE monitor
      USE nc_subs, ONLY : handle_err

      IMPLICIT NONE

#ifdef use_netcdf
      INCLUDE 'netcdf.inc'
*
*     Subroutine arguments

*     Local parameters
      character (len=*), parameter :: subnam = 'monnc_out'

!!  Netcdf variables used locally
      integer ncstat, startt, countt, start(2), count(2)

!!  Store current time as part of 'time' vector
      startt = ntdone/nocmon + 1
      countt = 1
      ncstat=nf_put_vara_double(monncid, time_id, startt, countt, tyrs)
      if ( ncstat.ne.NF_NOERR ) call handle_err (ncstat, subnam)

#  ifndef ocean_only
*     Store atmosphere diagnostics from monnc_comp
*     --------------------------------------------
!!  Store wetmat in 'wetmat'
      ncstat=nf_put_vara_double(monncid,wetmat_id,startt,countt,wetmat)
      if ( ncstat.ne.NF_NOERR ) call handle_err (ncstat, subnam)

!!  Store watmat in 'watmat'
      ncstat=nf_put_vara_double(monncid,watmat_id,startt,countt,watmat)
      if ( ncstat.ne.NF_NOERR ) call handle_err (ncstat, subnam)

!!  Store wepmat in 'wepmat'
      ncstat=nf_put_vara_double(monncid,wepmat_id,startt,countt,wepmat)
      if ( ncstat.ne.NF_NOERR ) call handle_err (ncstat, subnam)

!!  Store wapmat in 'wapmat'
      ncstat=nf_put_vara_double(monncid,wapmat_id,startt,countt,wapmat)
      if ( ncstat.ne.NF_NOERR ) call handle_err (ncstat, subnam)

!!  Store  entmat in 'entmat'
      start(1) = 1
      start(2) = ntdone/nocmon + 1
      count(1) = nla-1
      count(2) = 1
      ncstat=nf_put_vara_double(monncid,entmat_id,start,count,entmat)
      if ( ncstat.ne.NF_NOERR ) call handle_err (ncstat, subnam)

!!  Store enamat in 'enamat'
      start(1) = 1
      start(2) = ntdone/nocmon + 1
      count(1) = nla-1
      count(2) = 1
      ncstat=nf_put_vara_double(monncid,enamat_id,start,count,enamat)
      if ( ncstat.ne.NF_NOERR ) call handle_err (ncstat, subnam)

!!  Store etamat in 'etamat'
      start(1) = 1
      start(2) = ntdone/nocmon + 1
      count(1) = nla-1
      count(2) = 1
      ncstat=nf_put_vara_double(monncid,etamat_id,start,count,etamat)
      if ( ncstat.ne.NF_NOERR ) call handle_err (ncstat, subnam)

!!  Store ddtkeat in 'ddtkeat'
      start(1) = 1
      start(2) = ntdone/nocmon + 1
      count(1) = nla
      count(2) = 1
      ncstat=nf_put_vara_double(monncid,ddtkeat_id,start,count,ddtkeat)
      if ( ncstat.ne.NF_NOERR ) call handle_err (ncstat, subnam)

!!  Store ddtpeat in 'ddtpeat'
      start(1) = 1
      start(2) = ntdone/nocmon + 1
      count(1) = nla-1
      count(2) = 1
      ncstat=nf_put_vara_double(monncid,ddtpeat_id,start,count,ddtpeat)
      if ( ncstat.ne.NF_NOERR ) call handle_err (ncstat, subnam)

!!  Store et2mat in 'et2mat'
      start(1) = 1
      start(2) = ntdone/nocmon + 1
      count(1) = nla-1
      count(2) = 1
      ncstat=nf_put_vara_double(monncid,et2mat_id,start,count,et2mat)
      if ( ncstat.ne.NF_NOERR ) call handle_err (ncstat, subnam)

!!  Store pkenat in 'pkenat'
      start(1) = 1
      start(2) = ntdone/nocmon + 1
      count(1) = nla-1
      count(2) = 1
      ncstat=nf_put_vara_double(monncid,pkenat_id,start,count,pkenat)
      if ( ncstat.ne.NF_NOERR ) call handle_err (ncstat, subnam)

!!  Store utauat in 'utauat'
      ncstat=nf_put_vara_double(monncid,utauat_id,startt,countt,utauat)
      if ( ncstat.ne.NF_NOERR ) call handle_err (ncstat, subnam)

!!  Store atstpos in 'atstpos'
      start(1) = 1
      start(2) = ntdone/nocmon + 1
      count(1) = nla
      count(2) = 1
      ncstat=nf_put_vara_int(monncid,atstpos_id,start,count,atstpos)
      if ( ncstat.ne.NF_NOERR ) call handle_err (ncstat, subnam)

!!  Store atstval in 'atstval'
      start(1) = 1
      start(2) = ntdone/nocmon + 1
      count(1) = nla
      count(2) = 1
      ncstat=nf_put_vara_double(monncid,atstval_id,start,count,atstval)
      if ( ncstat.ne.NF_NOERR ) call handle_err (ncstat, subnam)

!!  Store pavgat in 'pavgat'
      start(1) = 1
      start(2) = ntdone/nocmon + 1
      count(1) = nla
      count(2) = 1
      ncstat=nf_put_vara_double(monncid,pavgat_id,start,count,pavgat)
      if ( ncstat.ne.NF_NOERR ) call handle_err (ncstat, subnam)

!!  Store qavgat in 'qavgat'
      start(1) = 1
      start(2) = ntdone/nocmon + 1
      count(1) = nla
      count(2) = 1
      ncstat=nf_put_vara_double(monncid,qavgat_id,start,count,qavgat)
      if ( ncstat.ne.NF_NOERR ) call handle_err (ncstat, subnam)

!!  Store ah4dat in 'ah4dat'
      start(1) = 1
      start(2) = ntdone/nocmon + 1
      count(1) = nla
      count(2) = 1
      ncstat=nf_put_vara_double(monncid,ah4dat_id,start,count,ah4dat)
      if ( ncstat.ne.NF_NOERR ) call handle_err (ncstat, subnam)

!!  Store kealat in 'kealat'
      start(1) = 1
      start(2) = ntdone/nocmon + 1
      count(1) = nla
      count(2) = 1
      ncstat=nf_put_vara_double(monncid,kealat_id,start,count,kealat)
      if ( ncstat.ne.NF_NOERR ) call handle_err (ncstat, subnam)

!!  Store astmin and astmax in 'astmin' and 'astmax'
      ncstat=nf_put_vara_double(monncid,astmin_id,startt,countt,astmin)
      if ( ncstat.ne.NF_NOERR ) call handle_err (ncstat, subnam)
      ncstat=nf_put_vara_double(monncid,astmax_id,startt,countt,astmax)
      if ( ncstat.ne.NF_NOERR ) call handle_err (ncstat, subnam)

!!  Store tmlmat in 'tmlmat'
      ncstat=nf_put_vara_double(monncid,tmlmat_id,startt,countt,tmlmat)
      if ( ncstat.ne.NF_NOERR ) call handle_err (ncstat, subnam)

!!  Store hmlmat in 'hmlmat'
      ncstat=nf_put_vara_double(monncid,hmlmat_id,startt,countt,hmlmat)
      if ( ncstat.ne.NF_NOERR ) call handle_err (ncstat, subnam)

!!  Store hcmlat in 'hcmlat'
      ncstat=nf_put_vara_double(monncid,hcmlat_id,startt,countt,hcmlat)
      if ( ncstat.ne.NF_NOERR ) call handle_err (ncstat, subnam)

!!  Store tmaooc in 'tmaooc'
      ncstat=nf_put_vara_double(monncid,tmaooc_id,startt,countt,tmaooc)
      if ( ncstat.ne.NF_NOERR ) call handle_err (ncstat, subnam)

!!  Store olrtop in 'olrtop'
      ncstat=nf_put_vara_double(monncid,olrtop_id,startt,countt,olrtop)
      if ( ncstat.ne.NF_NOERR ) call handle_err (ncstat, subnam)

!!  Store umminat, ummaxat, vmminat, vmmaxat in corresponding fields
      ncstat = nf_put_vara_double (monncid, umminat_id,
     &                             startt, countt, umminat)
      if ( ncstat.ne.NF_NOERR ) call handle_err (ncstat, subnam)
      ncstat = nf_put_vara_double (monncid, ummaxat_id,
     &                             startt, countt, ummaxat)
      if ( ncstat.ne.NF_NOERR ) call handle_err (ncstat, subnam)
      ncstat = nf_put_vara_double (monncid, vmminat_id,
     &                             startt, countt, vmminat)
      if ( ncstat.ne.NF_NOERR ) call handle_err (ncstat, subnam)
      ncstat = nf_put_vara_double (monncid, vmmaxat_id,
     &                             startt, countt, vmmaxat)
      if ( ncstat.ne.NF_NOERR ) call handle_err (ncstat, subnam)

!!  Store ugminat, ugmaxat, vgminat, vgmaxat in corresponding fields
      start(1) = 1
      start(2) = ntdone/nocmon + 1
      count(1) = nla
      count(2) = 1
      ncstat=nf_put_vara_double(monncid,ugminat_id,start,count,ugminat)
      if ( ncstat.ne.NF_NOERR ) call handle_err (ncstat, subnam)
      ncstat=nf_put_vara_double(monncid,ugmaxat_id,start,count,ugmaxat)
      if ( ncstat.ne.NF_NOERR ) call handle_err (ncstat, subnam)
      ncstat=nf_put_vara_double(monncid,vgminat_id,start,count,vgminat)
      if ( ncstat.ne.NF_NOERR ) call handle_err (ncstat, subnam)
      ncstat=nf_put_vara_double(monncid,vgmaxat_id,start,count,vgmaxat)
      if ( ncstat.ne.NF_NOERR ) call handle_err (ncstat, subnam)

!!  Store cnmlat in 'cnmlat'
      ncstat=nf_put_vara_double(monncid,cnmlat_id,startt,countt,cnmlat)
      if ( ncstat.ne.NF_NOERR ) call handle_err (ncstat, subnam)

!!  Store cnqgat in 'cnqgat'
      start(1) = 1
      start(2) = ntdone/nocmon + 1
      count(1) = nla
      count(2) = 1
      ncstat=nf_put_vara_double(monncid,cnqgat_id,start,count,cnqgat)
      if ( ncstat.ne.NF_NOERR ) call handle_err (ncstat, subnam)
#  endif /* not ocean_only */

#  ifndef atmos_only
*     Store ocean diagnostics from monnc_comp
*     ---------------------------------------
!!  Store wetmoc in 'wetmoc'
      ncstat=nf_put_vara_double(monncid,wetmoc_id,startt,countt,wetmoc)
      if ( ncstat.ne.NF_NOERR ) call handle_err (ncstat, subnam)

!!  Store watmoc in 'watmoc'
      ncstat=nf_put_vara_double(monncid,watmoc_id,startt,countt,watmoc)
      if ( ncstat.ne.NF_NOERR ) call handle_err (ncstat, subnam)

!!  Store wepmoc in 'wepmoc'
      ncstat=nf_put_vara_double(monncid,wepmoc_id,startt,countt,wepmoc)
      if ( ncstat.ne.NF_NOERR ) call handle_err (ncstat, subnam)

!!  Store wapmoc in 'wapmoc'
      ncstat=nf_put_vara_double(monncid,wapmoc_id,startt,countt,wapmoc)
      if ( ncstat.ne.NF_NOERR ) call handle_err (ncstat, subnam)

!!  Store entmoc in 'entmoc'
      ncstat=nf_put_vara_double(monncid,entmoc_id,startt,countt,entmoc)
      if ( ncstat.ne.NF_NOERR ) call handle_err (ncstat, subnam)

!!  Store enamoc in 'enamoc'
      ncstat=nf_put_vara_double(monncid,enamoc_id,startt,countt,enamoc)
      if ( ncstat.ne.NF_NOERR ) call handle_err (ncstat, subnam)

!!  Store pkenoc in 'pkenoc'
      ncstat=nf_put_vara_double(monncid,pkenoc_id,startt,countt,pkenoc)
      if ( ncstat.ne.NF_NOERR ) call handle_err (ncstat, subnam)

!!  Store etamoc in 'etamoc'
      start(1) = 1
      start(2) = ntdone/nocmon + 1
      count(1) = nlo-1
      count(2) = 1
      ncstat=nf_put_vara_double(monncid,etamoc_id,start,count,etamoc)
      if ( ncstat.ne.NF_NOERR ) call handle_err (ncstat, subnam)

!!  Store ddtkeoc in 'ddtkeoc'
      start(1) = 1
      start(2) = ntdone/nocmon + 1
      count(1) = nlo
      count(2) = 1
      ncstat=nf_put_vara_double(monncid,ddtkeoc_id,start,count,ddtkeoc)
      if ( ncstat.ne.NF_NOERR ) call handle_err (ncstat, subnam)

!!  Store ddtpeoc in 'ddtpeoc'
      start(1) = 1
      start(2) = ntdone/nocmon + 1
      count(1) = nlo-1
      count(2) = 1
      ncstat=nf_put_vara_double(monncid,ddtpeoc_id,start,count,ddtpeoc)
      if ( ncstat.ne.NF_NOERR ) call handle_err (ncstat, subnam)

!!  Store et2moc in 'et2moc'
      start(1) = 1
      start(2) = ntdone/nocmon + 1
      count(1) = nlo-1
      count(2) = 1
      ncstat=nf_put_vara_double(monncid,et2moc_id,start,count,et2moc)
      if ( ncstat.ne.NF_NOERR ) call handle_err (ncstat, subnam)

!!  Store utauoc in 'utauoc'
      ncstat=nf_put_vara_double(monncid,utauoc_id,startt,countt,utauoc)
      if ( ncstat.ne.NF_NOERR ) call handle_err (ncstat, subnam)

!!  Store ocjpos in 'ocjpos'
      start(1) = 1
      start(2) = ntdone/nocmon + 1
      count(1) = nlo
      count(2) = 1
      ncstat = nf_put_vara_int(monncid,ocjpos_id,start,count,ocjpos)
      if ( ncstat.ne.NF_NOERR ) call handle_err (ncstat, subnam)

!!  Store ocjval in 'ocjval'
      start(1) = 1
      start(2) = ntdone/nocmon + 1
      count(1) = nlo
      count(2) = 1
      ncstat = nf_put_vara_double(monncid,ocjval_id,start,count,ocjval)
      if ( ncstat.ne.NF_NOERR ) call handle_err (ncstat, subnam)

!!  Store pavgoc in 'pavgoc'
      start(1) = 1
      start(2) = ntdone/nocmon + 1
      count(1) = nlo
      count(2) = 1
      ncstat=nf_put_vara_double(monncid,pavgoc_id,start,count,pavgoc)
      if ( ncstat.ne.NF_NOERR ) call handle_err (ncstat, subnam)

!!  Store qavgoc in 'qavgoc'
      start(1) = 1
      start(2) = ntdone/nocmon + 1
      count(1) = nlo
      count(2) = 1
      ncstat=nf_put_vara_double(monncid,qavgoc_id,start,count,qavgoc)
      if ( ncstat.ne.NF_NOERR ) call handle_err (ncstat, subnam)

!!  Store ah2doc in 'ah2doc'
      start(1) = 1
      start(2) = ntdone/nocmon + 1
      count(1) = nlo
      count(2) = 1
      ncstat=nf_put_vara_double(monncid,ah2doc_id,start,count,ah2doc)
      if ( ncstat.ne.NF_NOERR ) call handle_err (ncstat, subnam)

!!  Store ah4doc in 'ah4doc'
      start(1) = 1
      start(2) = ntdone/nocmon + 1
      count(1) = nlo
      count(2) = 1
      ncstat=nf_put_vara_double(monncid,ah4doc_id,start,count,ah4doc)
      if ( ncstat.ne.NF_NOERR ) call handle_err (ncstat, subnam)

!!  Store kealoc in 'kealoc'
      start(1) = 1
      start(2) = ntdone/nocmon + 1
      count(1) = nlo
      count(2) = 1
      ncstat=nf_put_vara_double(monncid,kealoc_id,start,count,kealoc)
      if ( ncstat.ne.NF_NOERR ) call handle_err (ncstat, subnam)

!!  Store osfmin in 'osfmin'
      start(1) = 1
      start(2) = ntdone/nocmon + 1
      count(1) = nlo
      count(2) = 1
      ncstat = nf_put_vara_double(monncid,osfmin_id,start,count,osfmin)
      if ( ncstat.ne.NF_NOERR ) call handle_err (ncstat, subnam)

!!  Store osfmax in 'osfmax'
      start(1) = 1
      start(2) = ntdone/nocmon + 1
      count(1) = nlo
      count(2) = 1
      ncstat = nf_put_vara_double(monncid,osfmax_id,start,count,osfmax)
      if ( ncstat.ne.NF_NOERR ) call handle_err (ncstat, subnam)

!!  Store occirc in 'occirc'
      start(1) = 1
      start(2) = ntdone/nocmon + 1
      count(1) = nlo
      count(2) = 1
      ncstat = nf_put_vara_double(monncid,occirc_id,start,count,occirc)
      if ( ncstat.ne.NF_NOERR ) call handle_err (ncstat, subnam)

!!  Store btdgoc in 'btdgoc'
      ncstat=nf_put_vara_double(monncid,btdgoc_id,startt,countt,btdgoc)
      if ( ncstat.ne.NF_NOERR ) call handle_err (ncstat, subnam)

!!  Store sstmin and sstmax in 'sstmin' and 'sstmax'
      ncstat=nf_put_vara_double(monncid,sstmin_id,startt,countt,sstmin)
      if ( ncstat.ne.NF_NOERR ) call handle_err (ncstat, subnam)
      ncstat=nf_put_vara_double(monncid,sstmax_id,startt,countt,sstmax)
      if ( ncstat.ne.NF_NOERR ) call handle_err (ncstat, subnam)

!!  Store tmlmoc in 'tmlmoc'
      ncstat=nf_put_vara_double(monncid,tmlmoc_id,startt,countt,tmlmoc)
      if ( ncstat.ne.NF_NOERR ) call handle_err (ncstat, subnam)

!!  Store hfmloc in 'hfmloc'
      ncstat=nf_put_vara_double(monncid,hfmloc_id,startt,countt,hfmloc)
      if ( ncstat.ne.NF_NOERR ) call handle_err (ncstat, subnam)

!!  Store occtot in 'occtot'
      ncstat=nf_put_vara_double(monncid,occtot_id,startt,countt,occtot)
      if ( ncstat.ne.NF_NOERR ) call handle_err (ncstat, subnam)

!!  Store umminoc, ummaxoc, vmminoc, vmmaxoc in corresponding fields
      ncstat = nf_put_vara_double (monncid, umminoc_id,
     &                             startt, countt, umminoc)
      if ( ncstat.ne.NF_NOERR ) call handle_err (ncstat, subnam)
      ncstat = nf_put_vara_double (monncid, ummaxoc_id,
     &                             startt, countt, ummaxoc)
      if ( ncstat.ne.NF_NOERR ) call handle_err (ncstat, subnam)
      ncstat = nf_put_vara_double (monncid, vmminoc_id,
     &                             startt, countt, vmminoc)
      if ( ncstat.ne.NF_NOERR ) call handle_err (ncstat, subnam)
      ncstat = nf_put_vara_double (monncid, vmmaxoc_id,
     &                             startt, countt, vmmaxoc)
      if ( ncstat.ne.NF_NOERR ) call handle_err (ncstat, subnam)

!!  Store ugminoc, ugmaxoc, vgminoc, vgmaxoc in corresponding fields
      start(1) = 1
      start(2) = ntdone/nocmon + 1
      count(1) = nla
      count(2) = 1
      ncstat=nf_put_vara_double(monncid,ugminoc_id,start,count,ugminoc)
      if ( ncstat.ne.NF_NOERR ) call handle_err (ncstat, subnam)
      ncstat=nf_put_vara_double(monncid,ugmaxoc_id,start,count,ugmaxoc)
      if ( ncstat.ne.NF_NOERR ) call handle_err (ncstat, subnam)
      ncstat=nf_put_vara_double(monncid,vgminoc_id,start,count,vgminoc)
      if ( ncstat.ne.NF_NOERR ) call handle_err (ncstat, subnam)
      ncstat=nf_put_vara_double(monncid,vgmaxoc_id,start,count,vgmaxoc)
      if ( ncstat.ne.NF_NOERR ) call handle_err (ncstat, subnam)

!!  Store cnmloc in 'cnmloc'
      ncstat=nf_put_vara_double(monncid,cnmloc_id,startt,countt,cnmloc)
      if ( ncstat.ne.NF_NOERR ) call handle_err (ncstat, subnam)

!!  Store cnqgoc in 'cnqgoc'
      start(1) = 1
      start(2) = ntdone/nocmon + 1
      count(1) = nlo
      count(2) = 1
      ncstat=nf_put_vara_double(monncid,cnqgoc_id,start,count,cnqgoc)
      if ( ncstat.ne.NF_NOERR ) call handle_err (ncstat, subnam)

*     Store ocean diagnostics from oml
*     --------------------------------
!!  Store cfraoc in 'cfraoc'
      ncstat=nf_put_vara_double(monncid,cfraoc_id,startt,countt,cfraoc)
      if ( ncstat.ne.NF_NOERR ) call handle_err (ncstat, subnam)

!!  Store centoc in 'centoc'
      ncstat=nf_put_vara_double(monncid,centoc_id,startt,countt,centoc)
      if ( ncstat.ne.NF_NOERR ) call handle_err (ncstat, subnam)

*     Store ocean m. l. boundary diagnostics from omladf
*     --------------------------------------------------
!!  Store ttmads in 'ttmads'
      ncstat=nf_put_vara_double(monncid,ttmads_id,startt,countt,ttmads)
      if ( ncstat.ne.NF_NOERR ) call handle_err (ncstat, subnam)

!!  Store vfmads in 'vfmads'
      ncstat=nf_put_vara_double(monncid,vfmads_id,startt,countt,vfmads)
      if ( ncstat.ne.NF_NOERR ) call handle_err (ncstat, subnam)

!!  Store ttmdfs in 'ttmdfs'
      ncstat=nf_put_vara_double(monncid,ttmdfs_id,startt,countt,ttmdfs)
      if ( ncstat.ne.NF_NOERR ) call handle_err (ncstat, subnam)

!!  Store ttmadn in 'ttmadn'
      ncstat=nf_put_vara_double(monncid,ttmadn_id,startt,countt,ttmadn)
      if ( ncstat.ne.NF_NOERR ) call handle_err (ncstat, subnam)

!!  Store vfmadn in 'vfmadn'
      ncstat=nf_put_vara_double(monncid,vfmadn_id,startt,countt,vfmadn)
      if ( ncstat.ne.NF_NOERR ) call handle_err (ncstat, subnam)

!!  Store ttmdfn in 'ttmdfn'
      ncstat=nf_put_vara_double(monncid,ttmdfn_id,startt,countt,ttmdfn)
      if ( ncstat.ne.NF_NOERR ) call handle_err (ncstat, subnam)
#  endif /* not atmos_only */

#  ifndef ocean_only
*     Store flux diagnostics from xforc
*     ---------------------------------
!!  Store slhfav in 'slhfav'
      ncstat=nf_put_vara_double(monncid,slhfav_id,startt,countt,slhfav)
      if ( ncstat.ne.NF_NOERR ) call handle_err (ncstat, subnam)

!!  Store oradav in 'oradav'
      ncstat=nf_put_vara_double(monncid,oradav_id,startt,countt,oradav)
      if ( ncstat.ne.NF_NOERR ) call handle_err (ncstat, subnam)

!!  Store arocav in 'arocav'
      ncstat=nf_put_vara_double(monncid,arocav_id,startt,countt,arocav)
      if ( ncstat.ne.NF_NOERR ) call handle_err (ncstat, subnam)

!!  Store arlaav in 'arlaav'
      ncstat=nf_put_vara_double(monncid,arlaav_id,startt,countt,arlaav)
      if ( ncstat.ne.NF_NOERR ) call handle_err (ncstat, subnam)

*     Store atmosphere diagnostics from aml
*     -------------------------------------
!!  Store cfraat in 'cfraat'
      ncstat=nf_put_vara_double(monncid,cfraat_id,startt,countt,cfraat)
      if ( ncstat.ne.NF_NOERR ) call handle_err (ncstat, subnam)

!!  Store centat in 'centat'
      ncstat=nf_put_vara_double(monncid,centat_id,startt,countt,centat)
      if ( ncstat.ne.NF_NOERR ) call handle_err (ncstat, subnam)
#  endif /* not ocean_only */

#  ifndef ocean_only
*     Store consistency diagnostics from atinvq
*     -----------------------------------------
!!  Store ermasa in 'ermasa'
      start(1) = 1
      start(2) = ntdone/nocmon + 1
      count(1) = nla-1
      count(2) = 1
      ncstat=nf_put_vara_double(monncid,ermasa_id,start,count,ermasa)
      if ( ncstat.ne.NF_NOERR ) call handle_err (ncstat, subnam)

!!  Store emfrat in 'emfrat'
      start(1) = 1
      start(2) = ntdone/nocmon + 1
      count(1) = nla-1
      count(2) = 1
      ncstat=nf_put_vara_double(monncid,emfrat_id,start,count,emfrat)
      if ( ncstat.ne.NF_NOERR ) call handle_err (ncstat, subnam)
#  endif /* not ocean_only */

#  ifndef atmos_only
*     Store consistency diagnostics from ocinvq
*     -----------------------------------------
#    ifdef cyclic_ocean
!!  Store ermaso in 'ermaso'
      start(1) = 1
      start(2) = ntdone/nocmon + 1
      count(1) = nlo-1
      count(2) = 1
      ncstat=nf_put_vara_double(monncid,ermaso_id,start,count,ermaso)
      if ( ncstat.ne.NF_NOERR ) call handle_err (ncstat, subnam)

!!  Store emfroc in 'emfroc'
      start(1) = 1
      start(2) = ntdone/nocmon + 1
      count(1) = nlo-1
      count(2) = 1
      ncstat=nf_put_vara_double(monncid,emfroc_id,start,count,emfroc)
      if ( ncstat.ne.NF_NOERR ) call handle_err (ncstat, subnam)
#    endif /* cyclic_ocean */
#  endif /* not atmos_only */

#endif /* use_netcdf */

      END SUBROUTINE monnc_out

c***********************************************************************
*
      END MODULE mondiag
*
c***********************************************************************
