c***********************************************************************
c     Q-GCM Version 1.5.0 : last modified 07/06/2013
c
c    08/16/2020 (Sergey Kravtsov)
c
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
c
c
c    CONSTANT ATMOSPHERIC MIXED-LAYER THICKNESS VERSION
c
c              (S. Kravtsov, May 18 - August 16, 2020)
c
c    The version with a.m.l. emissivity zm<1 (August 16, 2020)
c
c***********************************************************************
c
c
c    MOIST MODEL INVOKED BY ACTIVATING THE moist FLAG 
c    (ORIGINAL DRY MODEL BY DEFAULT)
c
c    S. Kravtsov (10/09/2020)
c
c***********************************************************************
c
      MODULE amlsubs

*     Contains subroutines aml and amladf which step the atmospheric
*     mixed layer evolution equations (7.16) and (7.17), and update
*     the entrainment between atmospheric layers given by (7.18),
*     with possible convective correction given by (7.19).

*     Modules

      IMPLICIT NONE

      PRIVATE


#ifndef ocean_only
      PUBLIC  :: aml
      PRIVATE :: amladf

# ifdef moist

      PRIVATE :: imoadf,trapin

# endif

      CONTAINS
#endif

c***********************************************************************
*
#ifndef ocean_only

      SUBROUTINE aml
*
*     Timestep atmospheric mixed layer height - equation (7.16)
*     and mixed layer temperature - equation (7.17).
*     Also compute entrainment between atmospheric layers
*     - equation (7.18), plus convective correction (7.19).

*     Modules
      USE parameters, ONLY : nxta, nyta, nla, nxpa, nypa, atnorm,
     &                       nx1, ny1, nxaooc, nyaooc
      USE atconst, ONLY : dxa, dya, gpat, tat, dqat,rhoat,
     &                    tdta, xc1ast, xcexp, dtopat,
     &                    hat, hta, ytarel,tmbara,Lv,
     &                    rh,rhc,rhm,rhmc,
     &                    rho,rhom
      USE athomog, ONLY : xan, enisat, eninat
      USE atstate, ONLY : pa, pam, entat, entat1, wekta, 
     &                    precip_i, hum, humm
      USE intrfac, ONLY : ast, astm, ssta, fnetat, hmixa,
     &                    hmixam, hmat, hmamin,
     &                    evapa,precip,xfa0
      USE radiate, ONLY : rrcpat, aface, dface,
     &                    aface1, dface1, eface,  f2,
     &                    afacel, aface1l, dfacel, zm,
     &                    epsfac,gamma,
     &                    Flambar,Fmebar,Febar
      USE monitor, ONLY : cfraat, centat
      USE intsubs, ONLY : xintp
      USE xfosubs, ONLY : fsprim,esat

      IMPLICIT NONE


*     Subroutine arguments
*
*     Local variables
      integer i,j,l

      double precision hmainv,
     &                 entfac,dt1,dt2,
     &                 tmrhs(nxta,nyta),hmrhs(nxta,nyta),cfrasm,centsm,
     &              dhdiab,hnew,dhnew,dhfix,dtfix,trhtot,astnew,xfaent,
     &                 xfaent1,dtanew,conena,xfa(nxta,nyta),
     &                 ensums,ensumn,dtainv,xfasum,t_aux,
     &                 fsp,astmold(nxta,nyta),xfa1(nxta,nyta),
     &                 eta1,eta2,add_latent,xxf,xxf1,xfa2(nxta,nyta)


#  ifdef moist
      double precision Gc,g,RR,RV,eeps,Gtilde,tauqi
*
* tauqi is the inverse relaxation time scale for specific humidity (s^-1)
*
      integer nz
      parameter (nz=101,
     &           Gc=6.5d-3,g=9.82d0,RR=287.0d0,RV=461.d0,eeps=RR/RV,
     &           Gtilde=g/(RR*Gc),
     &           tauqi=0.d0/(24.d0*3600.d0*2.d0))

      double precision delz,zz,TT,TTs,fact,es,pd(nz),pp,t_aux1,h_ref,
     &                 dhfix1,hnew1,hnew2,hnew3,fact0,es0
      double precision humrhs(nxta,nyta,nla)
      double precision vmask(nyta),pmask(nyta),ppmask,tpi,yy8
      double precision A,B,C,D,q0,qp,qm,zp,zpp
      double precision dqw,add_latentm,add_latentm1
#  endif
*

        dt1=tat(2)-tat(1)
        dt2=tat(3)-tat(2) !Both e1 and e2 <0
        dtainv=1.d0/(dt1+f2*dt2)

# ifdef moist

      tpi=8.d0*atan(1.d0) !pi
      do j=1,nyta
        yy8=dble(j-1)
c         vmask(j)=0.5d0*(1.d0-tanh(dble(j-15)))
c         vmask(j)=1.5d0*cos(tpi*yy8/25.d0)*exp(-yy8**2/100.d0);
c         pmask(j)=0.5d0*(1.d0+tanh(0.5d0*(yy8-20.d0)))
        pmask(j)=1.d0
        vmask(j)=1.d0
      enddo

c
c    Make humidity gradient specification interactive
c
      
c      do l=1,nla-1
c         dqw=0.d0
c
c!$OMP PARALLEL DEFAULT (NONE)
c!$OMP&           PRIVATE (i,j)
c!$OMP&           SHARED  (l,dqw,humm)
c
c!$OMP   DO  SCHEDULE (STATIC)
c!$OMP&   REDUCTION(+:dqw)
c        do j=1,nyta
c          do i=1,nxta
c             dqw = dqw+humm(i,j,l+1)-humm(i,j,l)
c          enddo
c        enddo 
c!$OMP   END DO
c
c!$OMP END PARALLEL
c        dqat(l)=dqw*atnorm
c      
c      enddo
c
# endif


      hmainv = 1.0d0/hmat
      entfac = 1.0d0/( tdta*(tat(2)-tat(1)) )

*     Entrainment factors: aface(l),  cface and dface such that
*     e(1) = Sum(l)[ aface(l)*eta(l) ]  + cface*aD + dface*aTm'
*     We are assuming all entrainment is across interface 1.

c----------------------------------------------------------------------------------
c----------------------------------------------------------------------------------
c
c     The above formulas are modified here to allow for the entraintment through
c     interface 2; also, the a.m.l. thickness is assumed to be constant and
c     topography contributions to anomalous radiative fluxes are neglected,
c     as in the revised radsubs.F. There are SST contributions to radiative fluxes
c     over ocean and differences in radiative coefficients over ocean and land.
c     The entrainment over land has contributions  proportional to Fs  if
c     epsm (here, zm) <1.
c
c     e(2) = Sum(l)[ aface1(l)*eta(l) ] + dface1*aTm'
c     
c     Note: all cface are equal to zero.
c
c
c     (S. Kravtsov, 07/22/2020, 08/10/2020, 08/16/2020)
c-----------------------------------------------------------------------------------
c-----------------------------------------------------------------------------------



*     Initialise rh sides with advective, diffusive
*     and Del-4th terms - first 3 terms in equation
*     (7.16) and first four terms in equation (7.17)
*     ----------------------------------------------

*****************************************************************************************
*
* In amladf, hmrhs is now humidity advection/diffusion
* (if "moist" flag is activated) or, otherwise, zero) 
* (S. Kravtsov, 10/09/2020)
*
*****************************************************************************************
      call amladf (tmrhs, hmrhs, pa(1,1,1))

#  ifdef moist
*
*    Compute moisture advection/diffusion in interior layers
*
      do l=1,nla 
        call imoadf(humrhs(1,1,l),pa(1,1,l),hum(1,1,l),humm(1,1,l))
      enddo

#  endif
*     Add forcing term and diabatic effect - the 5th and
*     6th terms in equation (7.17), then timestep ast
*     Also do last term in equation (7.16) and step that
*     --------------------------------------------------
*     MONITORING - extra section for cfraat, centat
      cfrasm = 0.0d0
      centsm = 0.0d0

# ifdef moist
      add_latentm=0.d0
      add_latentm1=0.d0
# endif

!$OMP PARALLEL DEFAULT (NONE)
!$OMP&         PRIVATE (i,j,l,fsp)
!$OMP&         SHARED  (astm,hmixam,hmat,tat,tdta,hmrhs,
!$OMP&               hmamin,tmrhs,rrcpat,fnetat,hmainv,wekta,dface,
!$OMP&               dface1,xcexp,xc1ast,entfac,hmixa,xfa,xfa1,cfrasm,
!$OMP&               centsm,ast,aface,aface1,pam,entat,entat1,
!$OMP&               dtopat,dtainv,xfasum,f2,zm,ytarel,ssta,gpat,hta,
!$OMP&               hat,eface,astmold,dfacel,epsfac,afacel,aface1l)
#  ifdef moist
!$OMP&         SHARED  (evapa,precip,tmbara,gamma,rhoat,Lv,vmask,pmask,
!$OMP&                  humrhs,hum,humm,precip_i,dqat,rhm,rhmc,rh,rhc,
!$OMP&          rho,rhom,add_latentm,add_latentm1,Flambar,Fmebar,Febar,
!$OMP&          xfa0,xfa2)
#  endif

!$OMP DO SCHEDULE (STATIC)
!$OMP&   PRIVATE (dhdiab,hnew,dhnew,dhfix,dtfix,trhtot,eta1,eta2,
!$OMP&            astnew,xfaent,xfaent1,dtanew,conena,t_aux,add_latent)
#  ifdef moist
!$OMP&   PRIVATE (t_aux1,dhfix1,h_ref,delz,zz,TT,TTs,
!$OMP&               es,fact,es0,fact0,pd,pp,ppmask,
!$OMP&               A,B,C,D,q0,qp,qm,zp,zpp)
!$OMP&   REDUCTION(+:add_latentm)
#  endif
!$OMP&   REDUCTION(+:cfrasm) REDUCTION(-:centsm)

******************************************************************************************
*
* Below, comment out hmix step (S. Kravtsov, 05/18/2020)
*
******************************************************************************************
      do j=1,nyta
*       radiative forcing perturbation (for computing entrainment over land)
*
       fsp = fsprim( ytarel(j) )

# ifdef moist

       ppmask=pmask(j)

# endif

        do i=1,nxta

*******************************************************************************************
*
*  Mixed-layer entrainment contributions to ast tendency (use dhdiab workspace, subtract from trhtot)
*  and layer 1/2 entrainment (use dhfix workspace, subtract from xfaent)
*      (S. Kravtsov, 05/18/2020)
*
* Use reference temperature t_aux in mixed-layer entrainment formulas
*
*      (S. Kravtsov, 07/22/2020)
*
*******************************************************************************************


c
c     Atmospheric interface and temperature anomaly 
c

          eta1=0.25d0*(pam(i,j,1)-pam(i,j,2)
     &          +pam(i+1,j,1)-pam(i+1,j,2)
     &          +pam(i,j+1,1)-pam(i,j+1,2)
     &          +pam(i+1,j+1,1)-pam(i+1,j+1,2))/gpat(1)
          eta2=0.25d0*(pam(i,j,2)-pam(i,j,3)
     &         +pam(i+1,j,2)-pam(i+1,j,3)
     &         +pam(i,j+1,2)-pam(i,j+1,3)
     &         +pam(i+1,j+1,2)-pam(i+1,j+1,3))/gpat(2)

          t_aux=-(eta1*(tat(2)-tat(1))+eta2*(tat(3)-tat(2)))/hta

*
*      Step mixed-layer humidity equation, figure out 
*      the latent heat contribution to a.m.l.
*      thermal forcing (already included in the dry model)
*
*
      add_latent=0.d0

#  ifdef moist
c
c     Add the climatology to get the total temperature
c
          t_aux1=t_aux+tmbara
c
c         Compute reference specific humidity just above the a.m.l.
c         (assuming relative humidity of 0.7)
c          h_ref=0.7d-3*esat(t_aux1)

          
          TTs=t_aux1+Gc*hmat
          zz = hmat
          TT=TTs-Gc*zz
          fact0=(TTs/TT)**Gtilde !pressure factor
          es0=esat(TT)  !saturation water vapor pressure above a.m.l.

          zz = hmat+0.5d0*hat(1)
          TT=TTs-Gc*zz
          fact=(TTs/TT)**Gtilde !pressure factor
          es=esat(TT)  !saturation water vapor pressure at the lower-layer
                       !mid-level

          h_ref=humm(i,j,1)*es0*fact0/(es*fact) !assumed constant rh 
                                                ! in the lower layer
c
c    Compute entrainment contribution to the mixed-layer specific humidity
c
          if (wekta(i,j).le.0 ) then
            dhdiab = hmainv*wekta(i,j)
     &                     *( h_ref - hmixam(i,j) )
            dhfix1 = 0.d0
          else
            dhdiab = 0.d0
            dhfix1 = wekta(i,j)*(hmixam(i,j) - h_ref)
     &      /hat(1)
c     &                     *Lv*rhoat
          endif
*
*         Compute saturation specific humidity in the mixed layer
*         to implement relaxation to climatology
*
c         TTs=astm(i,j)+tmbara-0.5d0*(gamma-Gc)*hmat
c         zz=0.5d0*hmat
c         TT=TTs-Gc*zz
c         fact0=(TTs/TT)**Gtilde !pressure factor
c         es0=esat(TT)  !saturation water vapor pressure at the mixed-layer 
c                       !mid-level


*         Step mixed-layer moisture equation w/o precip: 
            
        hnew = hmixam(i,j) + tdta*(hmrhs(i,j)
     &        + hmainv*evapa(i,j)
     &        - hmainv*wekta(i,j)*hmixam(i,j)
     &        - dhdiab)
c     &     -tauqi*(hmixam(i,j)-rhm*0.001d0*fact0*es0))
c     &   -tauqi*vmask(j)*(hmixam(i,j)-fact0*es0*humm(i,j,1)/(fact*es)))

c        dhnew=
c     &     tauqi*vmask(j)*(fact*es*hmixam(i,j)/(fact0*es0)-humm(i,j,1))

c
c    Now compute precipitation and adjust mixed-layer
c    specific humidity accordingly. Use Laine et al. '11/
c    Deremble et al. '12 idea, but with a constant lapse rate
c    temperature profile within the mixed layer
c

c      delz = hmat/dble(nz-1)
c      TTs=astm(i,j)+tmbara-0.5d0*(gamma-Gc)*hmat
c        do l=1,nz
c          zz = dble(l-1)*delz
c          TT=TTs-Gc*zz
c          fact=(TTs/TT)**Gtilde
c          es=esat(TT)
c          pd(l) = max(hnew-0.001d0*fact*es,0.d0)
c        enddo
c      call trapin (pd, nz, delz, pp)
c
c

*--------------------------------------------------------------------------------
*     Precip computation using local quadratic fits to qs
*
      TTs=astm(i,j)+tmbara-0.5d0*(gamma-Gc)*hmat
      zz=0.5d0*hmat
      TT=TTs-Gc*zz
      fact=(TTs/TT)**Gtilde
      es=esat(TT)
      q0=0.001d0*eeps*fact*es
c
      zz=zz+1.d0
      TT=TTs-Gc*zz
      fact=(TTs/TT)**Gtilde
      es=esat(TT)
      qp=0.001d0*eeps*fact*es
c
      zz=zz-2.d0
      TT=TTs-Gc*zz
      fact=(TTs/TT)**Gtilde
      es=esat(TT)
      qm=0.001d0*eeps*fact*es
c
      C=hnew-q0
      B=0.5d0*(qm-qp)
      A=0.5d0*(2.d0*q0-qp-qm)
c
      D=B*B-4.d0*A*C
      pp=0.d0
      if(D.ge.0.d0)then
         if(A.ne.0.d0)then
           zp=0.5d0*(-B+sqrt(D))/A
         else
           zp=-C/B
         endif
         zp=max(zp,-0.5d0*hmat) 
         zpp=0.5d0*hmat
         zp=min(zpp,zp)
         pp=A*zpp*zpp*zpp/3.d0+B*zpp*zpp/2.d0+C*zpp
     &  -(A*zp*zp*zp/3.d0+B*zp*zp/2.d0+C*zp)
      endif
c
      hnew=hnew-pp*hmainv
      precip(i,j)=pp/tdta
      add_latent=rhom*Lv*ppmask*precip(i,j)
      add_latentm=add_latentm+add_latent*atnorm
*
*     End of precip computation with local quadratic fits
*===================================================================================

# endif

          if (wekta(i,j).le.0 ) then
            dhdiab = hmainv*wekta(i,j)
     &                     *( t_aux - astm(i,j) )
            dhfix = 0.d0
          else
            dhdiab = 0.d0
            dhfix = wekta(i,j)*( astm(i,j) - t_aux )
     &                     *dtainv
          endif
c
     
*         Predict new ast - equation (7.17)
          trhtot =  tmrhs(i,j) 
     &             + hmainv*rrcpat*(fnetat(i,j)+add_latent)
     &             - hmainv*wekta(i,j)*astm(i,j)
     &             - dhdiab
          astnew = astm(i,j) + tdta*trhtot

*         Find layer 1/2, 2/3 entrainment at T points 
*         (corrections over ocean are added later)
*
          xfaent = 
     &   + (dface+dfacel)*( xcexp*astm(i,j) + xc1ast(i,j) )
     &   +dhfix
     &   -fsp*(1.d0-zm)*epsfac*dtainv*rrcpat

          xfaent1 =  f2*xfaent
          xfaent = xfaent     +(aface(1)+afacel(1))*eta1
     &                        +(aface(2)+afacel(2))*eta2
          xfaent1 = xfaent1 +(aface1(1)+aface1l(1))*eta1
     &                      +(aface1(2)+aface1l(2))*eta2


********************************************************************************************
*
*   Skip convection parameterization below and allow perturbation 
*   mixed-layer temperature be warmer
*   than layer-1 temperature. Convection is fishy in a perturbation layer model
*   like this. (S. Kravtsov, 07/22/2020)
*
*
*********************************************************************************************

*         Check for convection & if necessary correct layer 1/2
*         entrainment and mixed layer temperature - equation (7.19)
*         dtanew should be >= 0 (stable case)
*         Correction is nonzero only if dtanew < 0

c          t_aux=0.25d0*(tat(2)-tat(1))*(pam(i,j,2)-pam(i,j,1)
c     &                          +pam(i+1,j,2)-pam(i+1,j,1)
c     &                          +pam(i,j+1,2)-pam(i,j+1,1)
c     &                          +pam(i+1,j+1,2)-pam(i+1,j+1,1))/(gpat(1)*hta)
c     &         +0.25d0*(tat(3)-tat(2))*(pam(i,j,3)-pam(i,j,2)
c     &                          +pam(i+1,j,3)-pam(i+1,j,2)
c     &                          +pam(i,j+1,3)-pam(i,j+1,2)
c     &                          +pam(i+1,j+1,3)-pam(i+1,j+1,2))/(gpat(2)*hta)
c          dtanew = tat(1) - astnew
c          dtanew = t_aux - astnew
c          conena = entfac*hmixa(i,j)*min( 0.0d0, dtanew )
c          xfa(i,j) = xfaent - xcexp*conena

          xfa(i,j) = xfaent
          xfa1(i,j)= xfaent1

c          astnew = astnew + min( 0.0d0, dtanew )

*         MONITORING - extra section for cfraat, centat
c          cfrasm = cfrasm + ( 0.5d0 - sign( 0.5d0, dtanew ) )
c          centsm = centsm - conena

          astmold(i,j)=astm(i,j)

          astm(i,j) = ast(i,j)
          ast(i,j) = astnew
#  ifdef moist
          hmixam(i,j) = hmixa(i,j)
          hmixa(i,j) = max(hnew,0.d0)
          humrhs(i,j,1)=humrhs(i,j,1)+dhfix1
c     &    +dhnew 
                                           !added moisture entrainment and convection flux 
                                           !to the lower layer moisture equation
#  endif
        enddo
      enddo
!$OMP END DO
*
* Subtract the basin-mean latent-heat release to take out the climatological release
*

# ifdef moist

!$OMP DO SCHEDULE (STATIC)
        do j=1,nyta
          do i=1,nxta
            ast(i,j)=ast(i,j)-tdta*hmainv*rrcpat*add_latentm
c            ast(i,j)=ast(i,j)-tdta*hmainv*rrcpat*(Flambar-Fmebar)
          enddo
        enddo
!$OMP END DO NOWAIT

# endif

*-------------------------------------------------------------------------------------------
*
*      Now include xfa corrections over ocean,
*      due to a different radiation formulation there:
*
*


!$OMP DO SCHEDULE (STATIC)
!$OMP&   PRIVATE (eta1,eta2)

      do j=ny1,ny1+nyaooc-1
*       Compute radiative forcing perturbation
        fsp = fsprim( ytarel(j) )
        do i=nx1,nx1+nxaooc-1

          eta1=0.25d0*(pam(i,j,1)-pam(i,j,2)
     &          +pam(i+1,j,1)-pam(i+1,j,2)
     &          +pam(i,j+1,1)-pam(i,j+1,2)
     &          +pam(i+1,j+1,1)-pam(i+1,j+1,2))/gpat(1)
          eta2=0.25d0*(pam(i,j,2)-pam(i,j,3)
     &         +pam(i+1,j,2)-pam(i+1,j,3)
     &         +pam(i,j+1,2)-pam(i,j+1,3)
     &         +pam(i+1,j+1,2)-pam(i+1,j+1,3))/gpat(2)

          xfa(i,j) = xfa(i,j)
     &      - dfacel*( xcexp*astmold(i,j) + xc1ast(i,j) )
     &      +fsp*(1.d0-zm)*epsfac*dtainv*rrcpat
     &      + eface*ssta(i-nx1+1,j-ny1+1)
     &      -afacel(1)*eta1-afacel(2)*eta2

          xfa1(i,j) = xfa1(i,j)+f2*
     &      (- dfacel*( xcexp*astmold(i,j) + xc1ast(i,j) )
     &      +fsp*(1.d0-zm)*epsfac*dtainv*rrcpat
     &      + eface*ssta(i-nx1+1,j-ny1+1))
     &      -aface1l(1)*eta1-aface1l(2)*eta2

        enddo
      enddo

!$OMP END DO NOWAIT


*------------------------------------------------------------------------------------
*
*   Step interior moisture equations, compute layer precipitation,
*   modify interior entrainments to adjust for latent heat release
*
*

# ifdef moist


c
c First iteration, with "dry" entrainment
c

!$OMP DO SCHEDULE (STATIC)
!$OMP&   REDUCTION(+:add_latentm1)
!$OMP&   PRIVATE (eta1,eta2,t_aux,t_aux1,hnew1,hnew2,hnew3,
!$OMP&            delz,zz,TT,TTs,es,fact,es0,fact0,
!$OMP&            pd,pp,add_latent,ppmask,
!$OMP&            A,B,C,D,q0,qp,qm,zp,zpp,xxf,xxf1)

      do j=1,nyta
         ppmask=pmask(j)
        do i=1,nxta

c
c     Atmospheric interface and temperature anomaly
c

          eta1=0.25d0*(pam(i,j,1)-pam(i,j,2)
     &          +pam(i+1,j,1)-pam(i+1,j,2)
     &          +pam(i,j+1,1)-pam(i,j+1,2)
     &          +pam(i+1,j+1,1)-pam(i+1,j+1,2))/gpat(1)
          eta2=0.25d0*(pam(i,j,2)-pam(i,j,3)
     &         +pam(i+1,j,2)-pam(i+1,j,3)
     &         +pam(i,j+1,2)-pam(i,j+1,3)
     &         +pam(i+1,j+1,2)-pam(i+1,j+1,3))/gpat(2)

          t_aux=-(eta1*(tat(2)-tat(1))+eta2*(tat(3)-tat(2)))/hta
c
c     Add the climatology to get the total temperature
c
          t_aux1=t_aux+tmbara
c
c  Add entrainment correction due to precip from the previous time step
c
         xxf=xfa(i,j)
         xxf1=xfa1(i,j)
c

*
*         Step lower-layer moisture equation w/o precip:
*

c          TTs=t_aux1+Gc*hmat
c          zz = hmat+0.5d0*hat(1)
c          TT=TTs-Gc*zz
c          fact0=(TTs/TT)**Gtilde !pressure factor
c          es0=esat(TT)  !saturation water vapor pressure at the lower-layer
c                       !mid-level
c
c          zz = hmat+hat(1)+0.5d0*hat(2)
c          TT=TTs-Gc*zz
c          fact=(TTs/TT)**Gtilde !pressure factor
c          es=esat(TT)  !saturation water vapor pressure at the middle-layer
c                       !mid-level



          hnew1 = humm(i,j,1) + tdta*(humrhs(i,j,1)
c     &               + min(xxf,0.d0)*(rho(2)/rho(1))*dqat(1)/hat(1))
     & +min(xxf,0.d0)*(humm(i,j,1)-humm(i,j,2)*rho(2)/rho(1))/hat(1))
c     &            -tauqi*(humm(i,j,1)-rh(1)*0.001d0*fact0*es0)) 
c     &   -tauqi*vmask(j)*(humm(i,j,1)-fact0*es0*humm(i,j,2)/(fact*es)))



c*
c*         compute precip and adjust accordingly
c*
c          pp=max(hnew1-rhc(1)*0.001d0*fact0*es0,0.d0)
c          hnew1=hnew1-pp
c          precip_i(i,j,1)=pp/tdta
c          add_latent=rhoat*Lv*pp/tdta
c          hnew1=max(hnew1,0.d0)

*         Step mid-layer moisture equation w/o precip:

          hnew2 = humm(i,j,2) + tdta*(humrhs(i,j,2)
c     &                + (max(xxf,0.d0)*(rho(1)/rho(2))*dqat(1)
c     &            + min(xxf1,0.d0)*(rho(3)/rho(2))*dqat(2))/hat(2))
     & +(max(xxf,0.d0)*(humm(i,j,1)*rho(1)/rho(2)-humm(i,j,2))
     & +min(xxf1,0.d0)*(humm(i,j,2)-humm(i,j,3)*rho(3)/rho(2)))/hat(2))

c     &            -tauqi*(humm(i,j,2)-rh(2)*0.001d0*fact*es))
c     &   -tauqi*vmask(j)*(humm(i,j,2)-fact*es*humm(i,j,1)/(fact0*es0)))


c          zz = hmat+hat(1)+hat(2)+0.5d0*hat(3)
c          TT=TTs-Gc*zz
c          fact0=(TTs/TT)**Gtilde !pressure factor
c          es0=esat(TT)  !saturation water vapor pressure at the upper-layer
c                       !mid-level

c          hnew2=hnew2+tdta*
c     &  (-tauqi*vmask(j)*(humm(i,j,2)-fact*es*humm(i,j,3)/(fact0*es0)))

c*
c*         compute precip and adjust accordingly
c*
c          pp=max(hnew2-rhc(2)*0.001d0*fact*es,0.d0)
c          hnew2=hnew2-pp
c          precip_i(i,j,2)=pp/tdta
c          add_latent=add_latent+rhoat*Lv*pp/tdta
c          hnew2=max(hnew2,0.d0)

*         Step upper-layer moisture equation w/o precip:

c          zz = hmat+hat(1)+hat(2)+0.5d0*hat(3)
c          TT=TTs-Gc*zz
c          fact=(TTs/TT)**Gtilde !pressure factor
c          es=esat(TT)  !saturation water vapor pressure at the middle-layer
c                       !mid-level

          hnew3 = humm(i,j,3) + tdta*(humrhs(i,j,3)
c     &             + max(xxf1,0.d0)*(rho(2)/rho(3))*dqat(2)/hat(3))
     & + max(xxf1,0.d0)*(humm(i,j,2)*rho(2)/rho(3)-humm(i,j,3))/hat(3))
c     &            -tauqi*(humm(i,j,3)-rh(3)*0.001d0*fact*es))
c     &   -tauqi*vmask(j)*(humm(i,j,3)-fact0*es0*humm(i,j,2)/(fact*es)))

c*
c*         compute precip and adjust accordingly
c*
c          pp=max(hnew3-rhc(3)*0.001d0*fact*es,0.d0)
c          hnew3=hnew3-pp
c          precip_i(i,j,3)=pp/tdta
c          add_latent=add_latent+rhoat*Lv*pp/tdta
c          hnew3=max(hnew3,0.d0)

c
c    Now compute precipitation and adjust
c    humidities accordingly. Use Laine et al. '11/
c    Deremble et al. '12 idea, but with a constant lapse rate
c    temperature profile
c

c      delz = hat(1)/dble(nz-1)
cc      TTs=t_aux1
c        do l=1,nz
c          zz = hmat+dble(l-1)*delz
c          TT=TTs-Gc*zz
c          fact=(TTs/TT)**Gtilde
c          es=esat(TT)
c          pd(l) = max(hnew1-0.001d0*fact*es,0.d0)
c        enddo
c      call trapin (pd, nz, delz, pp)
c      hnew1=hnew1-pp/hat(1)
c      precip_i(i,j,1)=pp/tdta
c      add_latent=rhoat*Lv*pp/tdta
c
c      delz = hat(2)/dble(nz-1)
c        do l=1,nz
c          zz = hmat+hat(1)+dble(l-1)*delz
c          TT=TTs-Gc*zz
c          fact=(TTs/TT)**Gtilde
c          es=esat(TT)
c          pd(l) = max(hnew2-0.001d0*fact*es,0.d0)
c        enddo
c      call trapin (pd, nz, delz, pp)
c      hnew2=hnew2-pp/hat(2)
c      precip_i(i,j,2)=pp/tdta
c      add_latent=add_latent+rhoat*Lv*pp/tdta
c
c      delz = hat(3)/dble(nz-1)
c        do l=1,nz
c          zz = hmat+hat(1)+hat(2)+dble(l-1)*delz
c          TT=TTs-Gc*zz
c          fact=(TTs/TT)**Gtilde
c          es=esat(TT)
c          pd(l) = max(hnew3-0.001d0*fact*es,0.d0)
c        enddo
c      call trapin (pd, nz, delz, pp)
c      hnew3=hnew3-pp/hat(3)
c      precip_i(i,j,3)=pp/tdta
c      add_latent=add_latent+rhoat*Lv*pp/tdta
c
c      add_latent=add_latent*ppmask
c
*--------------------------------------------------------------------------------
*     Precip computation using local quadratic fits to qs, layer 1
*

      TTs=t_aux1-0.5d0*(gamma-Gc)*hmat

      zz=hmat+0.5d0*hat(1)
      TT=TTs-Gc*zz
      fact=(TTs/TT)**Gtilde
      es=esat(TT)
      q0=0.001d0*eeps*fact*es
c
      zz=zz+1.d0
      TT=TTs-Gc*zz
      fact=(TTs/TT)**Gtilde
      es=esat(TT)
      qp=0.001d0*eeps*fact*es
c
      zz=zz-2.d0
      TT=TTs-Gc*zz
      fact=(TTs/TT)**Gtilde
      es=esat(TT)
      qm=0.001d0*eeps*fact*es

      C=hnew1-q0
      B=0.5d0*(qm-qp)
      A=0.5d0*(2.d0*q0-qp-qm)

      D=B*B-4.d0*A*C
      pp=0.d0
      if(D.ge.0.d0)then
         if(A.ne.0.d0)then
           zp=0.5d0*(-B+sqrt(D))/A
         else
           zp=-C/B
         endif
         zp=max(zp,-0.5d0*hat(1))
         zpp=0.5d0*hat(1)
         zp=min(zpp,zp)
         pp=A*zpp*zpp*zpp/3.d0+B*zpp*zpp/2.d0+C*zpp
     &  -(A*zp*zp*zp/3.d0+B*zp*zp/2.d0+C*zp)
      endif

      hnew1=hnew1-pp/hat(1)
      precip_i(i,j,1)=pp/tdta
      add_latent=rho(1)*Lv*pp/tdta

*
*     End of precip computation with local quadratic fits, layer 1
*===================================================================================

*--------------------------------------------------------------------------------
*     Precip computation using local quadratic fits to qs, layer 2
*
      zz=hmat+hat(1)+0.5d0*hat(2)
      TT=TTs-Gc*zz
      fact=(TTs/TT)**Gtilde
      es=esat(TT)
      q0=0.001d0*eeps*fact*es
c
      zz=zz+1.d0
      TT=TTs-Gc*zz
      fact=(TTs/TT)**Gtilde
      es=esat(TT)
      qp=0.001d0*eeps*fact*es
c
      zz=zz-2.d0
      TT=TTs-Gc*zz
      fact=(TTs/TT)**Gtilde
      es=esat(TT)
      qm=0.001d0*eeps*fact*es

      C=hnew2-q0
      B=0.5d0*(qm-qp)
      A=0.5d0*(2.d0*q0-qp-qm)

      D=B*B-4.d0*A*C
      pp=0.d0
      if(D.ge.0.d0)then
         if(A.ne.0.d0)then
           zp=0.5d0*(-B+sqrt(D))/A
         else
           zp=-C/B
         endif
         zp=max(zp,-0.5d0*hat(2))
         zpp=0.5d0*hat(2)
         zp=min(zpp,zp)
         pp=A*zpp*zpp*zpp/3.d0+B*zpp*zpp/2.d0+C*zpp
     &  -(A*zp*zp*zp/3.d0+B*zp*zp/2.d0+C*zp)
      endif

      hnew2=hnew2-pp/hat(2)
      precip_i(i,j,2)=pp/tdta
      add_latent=add_latent+rho(2)*Lv*pp/tdta

*
*     End of precip computation with local quadratic fits, layer 2
*===================================================================================

*--------------------------------------------------------------------------------
*     Precip computation using local quadratic fits to qs, layer 3
*
      zz=hmat+hat(1)+hat(2)+0.5d0*hat(3)
      TT=TTs-Gc*zz
      fact=(TTs/TT)**Gtilde
      es=esat(TT)
      q0=0.001d0*eeps*fact*es
c
      zz=zz+1.d0
      TT=TTs-Gc*zz
      fact=(TTs/TT)**Gtilde
      es=esat(TT)
      qp=0.001d0*eeps*fact*es
c
      zz=zz-2.d0
      TT=TTs-Gc*zz
      fact=(TTs/TT)**Gtilde
      es=esat(TT)
      qm=0.001d0*eeps*fact*es

      C=hnew3-q0
      B=0.5d0*(qm-qp)
      A=0.5d0*(2.d0*q0-qp-qm)

      D=B*B-4.d0*A*C
      pp=0.d0
      if(D.ge.0.d0)then
         if(A.ne.0.d0)then
           zp=0.5d0*(-B+sqrt(D))/A
         else
           zp=-C/B
         endif
         zp=max(zp,-0.5d0*hat(3))
         zpp=0.5d0*hat(3)
         zp=min(zpp,zp)
         pp=A*zpp*zpp*zpp/3.d0+B*zpp*zpp/2.d0+C*zpp
     &  -(A*zp*zp*zp/3.d0+B*zp*zp/2.d0+C*zp)
      endif

      hnew3=hnew3-pp/hat(3)
      precip_i(i,j,3)=pp/tdta
      add_latent=add_latent+rho(3)*Lv*pp/tdta
c      add_latent=Lv*(rho(1)*precip_i(i,j,1)
c     &              +rho(2)*precip_i(i,j,2)
c     &              +rho(3)*precip_i(i,j,3))

*
*     End of precip computation with local quadratic fits, layer 3
*===================================================================================

c       add_latent=add_latent-Fmebar !subtract climatological flux
c      add_latent=add_latent*ppmask
       add_latentm1=add_latentm1+add_latent*atnorm

         hnew1=max(hnew1,0.d0)
         hnew2=max(hnew2,0.d0)
         hnew3=max(hnew3,0.d0)

*
*     update humidities:
*

# ifdef dry_latent

          humm(i,j,1) = hum(i,j,1)
          hum(i,j,1) = hnew1
c
          humm(i,j,2) = hum(i,j,2)
          hum(i,j,2) = hnew2
c
          humm(i,j,3) = hum(i,j,3)
          hum(i,j,3) = hnew3

*
*     Now adjust entrainment rates:
*

      xfa(i,j)=xfa(i,j)+add_latent*rrcpat*dtainv
      xfa1(i,j)=xfa1(i,j)+f2*add_latent*rrcpat*dtainv
# else
      xfa0(i,j)=add_latent*rrcpat*dtainv
# endif

        enddo
      enddo

c      print*,'********************************' 

!$OMP END DO 

!$OMP DO SCHEDULE (STATIC)
      do j=1,nyta
        do i=1,nxta
*
*     Now remove climatological effect:
*

# ifdef dry_latent
      xfa(i,j)=xfa(i,j)-add_latentm1*rrcpat*dtainv
      xfa1(i,j)=xfa1(i,j)-f2*add_latentm1*rrcpat*dtainv
# else
      xfa0(i,j)=xfa0(i,j)-add_latentm1*rrcpat*dtainv
# endif

        enddo
      enddo

!$OMP END DO NOWAIT

c-------------------------------------------------------
c SECOND ITERATION, WITH FIRST-ORDER PRECIP CORRECTION
c-------------------------------------------------------

# ifndef dry_latent

!$OMP SINGLE
       add_latentm1=0.d0
!$OMP END SINGLE


!$OMP DO SCHEDULE (STATIC)
!$OMP&   REDUCTION(+:add_latentm1)
!$OMP&   PRIVATE (eta1,eta2,t_aux,t_aux1,hnew1,hnew2,hnew3,
!$OMP&            delz,zz,TT,TTs,es,fact,es0,fact0,
!$OMP&            pd,pp,add_latent,ppmask,
!$OMP&            A,B,C,D,q0,qp,qm,zp,zpp,xxf,xxf1)

      do j=1,nyta
         ppmask=pmask(j)
        do i=1,nxta

c
c     Atmospheric interface and temperature anomaly
c

          eta1=0.25d0*(pam(i,j,1)-pam(i,j,2)
     &          +pam(i+1,j,1)-pam(i+1,j,2)
     &          +pam(i,j+1,1)-pam(i,j+1,2)
     &          +pam(i+1,j+1,1)-pam(i+1,j+1,2))/gpat(1)
          eta2=0.25d0*(pam(i,j,2)-pam(i,j,3)
     &         +pam(i+1,j,2)-pam(i+1,j,3)
     &         +pam(i,j+1,2)-pam(i,j+1,3)
     &         +pam(i+1,j+1,2)-pam(i+1,j+1,3))/gpat(2)

          t_aux=-(eta1*(tat(2)-tat(1))+eta2*(tat(3)-tat(2)))/hta
c
c     Add the climatology to get the total temperature
c
          t_aux1=t_aux+tmbara
c
c  Add entrainment correction due to precip from the previous iteration 
c
         xxf=xfa(i,j)+xfa0(i,j)
         xxf1=xfa1(i,j)+f2*xfa0(i,j)
c
c         print*,i,j,xfa(i,j),xfa0(i,j)
c
*
*         Step lower-layer moisture equation w/o precip:
*

          hnew1 = humm(i,j,1) + tdta*(humrhs(i,j,1)
     & + min(xxf,0.d0)*(humm(i,j,1)-humm(i,j,2)*rho(2)/rho(1))/hat(1))


*         Step mid-layer moisture equation w/o precip:

          hnew2 = humm(i,j,2) + tdta*(humrhs(i,j,2)
     & + (max(xxf,0.d0)*(humm(i,j,1)*rho(1)/rho(2)-humm(i,j,2))
     & +min(xxf1,0.d0)*(humm(i,j,2)-humm(i,j,3)*rho(3)/rho(2)))/hat(2))


*         Step upper-layer moisture equation w/o precip:

          hnew3 = humm(i,j,3) + tdta*(humrhs(i,j,3)
     & + max(xxf1,0.d0)*(humm(i,j,2)*rho(2)/rho(3)-humm(i,j,3))/hat(3))

c
c    Now compute precipitation and adjust
c    humidities accordingly. Use Laine et al. '11/
c    Deremble et al. '12 idea, but with a constant lapse rate
c    temperature profile
c
c
*--------------------------------------------------------------------------------
*     Precip computation using local quadratic fits to qs, layer 1
*

      TTs=t_aux1-0.5d0*(gamma-Gc)*hmat

      zz=hmat+0.5d0*hat(1)
      TT=TTs-Gc*zz
      fact=(TTs/TT)**Gtilde
      es=esat(TT)
      q0=0.001d0*eeps*fact*es
c
      zz=zz+1.d0
      TT=TTs-Gc*zz
      fact=(TTs/TT)**Gtilde
      es=esat(TT)
      qp=0.001d0*eeps*fact*es
c
      zz=zz-2.d0
      TT=TTs-Gc*zz
      fact=(TTs/TT)**Gtilde
      es=esat(TT)
      qm=0.001d0*eeps*fact*es

      C=hnew1-q0
      B=0.5d0*(qm-qp)
      A=0.5d0*(2.d0*q0-qp-qm)

      D=B*B-4.d0*A*C
      pp=0.d0
      if(D.ge.0.d0)then
         if(A.ne.0.d0)then
           zp=0.5d0*(-B+sqrt(D))/A
         else
           zp=-C/B
         endif
         zp=max(zp,-0.5d0*hat(1))
         zpp=0.5d0*hat(1)
         zp=min(zpp,zp)
         pp=A*zpp*zpp*zpp/3.d0+B*zpp*zpp/2.d0+C*zpp
     &  -(A*zp*zp*zp/3.d0+B*zp*zp/2.d0+C*zp)
      endif

      hnew1=hnew1-pp/hat(1)
      precip_i(i,j,1)=pp/tdta

*
*     End of precip computation with local quadratic fits, layer 1
*===================================================================================

*--------------------------------------------------------------------------------
*     Precip computation using local quadratic fits to qs, layer 2
*
      zz=hmat+hat(1)+0.5d0*hat(2)
      TT=TTs-Gc*zz
      fact=(TTs/TT)**Gtilde
      es=esat(TT)
      q0=0.001d0*eeps*fact*es
c
      zz=zz+1.d0
      TT=TTs-Gc*zz
      fact=(TTs/TT)**Gtilde
      es=esat(TT)
      qp=0.001d0*eeps*fact*es
c
      zz=zz-2.d0
      TT=TTs-Gc*zz
      fact=(TTs/TT)**Gtilde
      es=esat(TT)
      qm=0.001d0*eeps*fact*es

      C=hnew2-q0
      B=0.5d0*(qm-qp)
      A=0.5d0*(2.d0*q0-qp-qm)

      D=B*B-4.d0*A*C
      pp=0.d0
      if(D.ge.0.d0)then
         if(A.ne.0.d0)then
           zp=0.5d0*(-B+sqrt(D))/A
         else
           zp=-C/B
         endif
         zp=max(zp,-0.5d0*hat(2))
         zpp=0.5d0*hat(2)
         zp=min(zpp,zp)
         pp=A*zpp*zpp*zpp/3.d0+B*zpp*zpp/2.d0+C*zpp
     &  -(A*zp*zp*zp/3.d0+B*zp*zp/2.d0+C*zp)
      endif

      hnew2=hnew2-pp/hat(2)
      precip_i(i,j,2)=pp/tdta

*
*     End of precip computation with local quadratic fits, layer 2
*===================================================================================

*--------------------------------------------------------------------------------
*     Precip computation using local quadratic fits to qs, layer 3
*
      zz=hmat+hat(1)+hat(2)+0.5d0*hat(3)
      TT=TTs-Gc*zz
      fact=(TTs/TT)**Gtilde
      es=esat(TT)
      q0=0.001d0*eeps*fact*es
c
      zz=zz+1.d0
      TT=TTs-Gc*zz
      fact=(TTs/TT)**Gtilde
      es=esat(TT)
      qp=0.001d0*eeps*fact*es
c
      zz=zz-2.d0
      TT=TTs-Gc*zz
      fact=(TTs/TT)**Gtilde
      es=esat(TT)
      qm=0.001d0*eeps*fact*es

      C=hnew3-q0
      B=0.5d0*(qm-qp)
      A=0.5d0*(2.d0*q0-qp-qm)

      D=B*B-4.d0*A*C
      pp=0.d0
      if(D.ge.0.d0)then
         if(A.ne.0.d0)then
           zp=0.5d0*(-B+sqrt(D))/A
         else
           zp=-C/B
         endif
         zp=max(zp,-0.5d0*hat(3))
         zpp=0.5d0*hat(3)
         zp=min(zpp,zp)
         pp=A*zpp*zpp*zpp/3.d0+B*zpp*zpp/2.d0+C*zpp
     &  -(A*zp*zp*zp/3.d0+B*zp*zp/2.d0+C*zp)
      endif

      hnew3=hnew3-pp/hat(3)
      precip_i(i,j,3)=pp/tdta
      add_latent=Lv*(rho(1)*precip_i(i,j,1)
     &              +rho(2)*precip_i(i,j,2)
     &              +rho(3)*precip_i(i,j,3))

*
*     End of precip computation with local quadratic fits, layer 3
*===================================================================================

       add_latentm1=add_latentm1+add_latent*atnorm

         hnew1=max(hnew1,0.d0)
         hnew2=max(hnew2,0.d0)
         hnew3=max(hnew3,0.d0)

*
*     update humidities:
*

          humm(i,j,1) = hum(i,j,1)
          hum(i,j,1) = hnew1

          humm(i,j,2) = hum(i,j,2)
          hum(i,j,2) = hnew2

          humm(i,j,3) = hum(i,j,3)
          hum(i,j,3) = hnew3

*
*     Now adjust entrainment rates:
*

      xfa(i,j)=xfa(i,j)+add_latent*rrcpat*dtainv
      xfa1(i,j)=xfa1(i,j)+f2*add_latent*rrcpat*dtainv
      xfa2(i,j)=add_latent*rrcpat*dtainv

        enddo
      enddo

c      print*,'********************************'

!$OMP END DO

!$OMP DO SCHEDULE (STATIC)
      do j=1,nyta
        do i=1,nxta
*
*     Now remove climatological effect:
*

      xfa(i,j)=xfa(i,j)-add_latentm1*rrcpat*dtainv
      xfa1(i,j)=xfa1(i,j)-f2*add_latentm1*rrcpat*dtainv
      xfa2(i,j)=xfa2(i,j)-add_latentm1*rrcpat*dtainv

        enddo
      enddo

!$OMP END DO NOWAIT

# endif  !for dry_latent

# endif


*     Average entrainment (computed in xfa at T points) onto p points
*     ---------------------------------------------------------------
*     Internal points
!$OMP DO SCHEDULE (STATIC)
      do j=2,nypa-1
*       Cyclicity condition
        entat(  1 ,j) = 0.25d0*(  xfa(nxta,j-1) + xfa(1,j-1)
     &                          + xfa(nxta, j ) + xfa(1, j ) )
        do i=2,nxpa-1
          entat(i,j) = 0.25d0*(  xfa(i-1,j-1) + xfa(i,j-1)
     &                         + xfa(i-1, j ) + xfa(i, j ) )
        enddo
        entat(nxpa,j) = entat(1,j)
      enddo
!$OMP END DO NOWAIT

*     Do edge and corner values to conserve area integral
!$OMP SINGLE
      entat(  1 ,  1 ) = 0.5d0*( xfa(nxta,  1 ) + xfa(1,  1 ) )
      entat(  1 ,nypa) = 0.5d0*( xfa(nxta,nyta) + xfa(1,nyta) )
      do i=2,nxpa-1
        entat(i,  1 ) = 0.5d0*( xfa(i-1,  1 ) + xfa(i,  1 ) )
        entat(i,nypa) = 0.5d0*( xfa(i-1,nyta) + xfa(i,nyta) )
      enddo
      entat(nxpa,  1 ) = entat(  1 ,  1 )
      entat(nxpa,nypa) = entat(  1 ,nypa)
!$OMP END SINGLE


!$OMP DO SCHEDULE (STATIC)
      do j=2,nypa-1
*       Cyclicity condition
        entat1(  1 ,j) = 0.25d0*(  xfa1(nxta,j-1) + xfa1(1,j-1)
     &                          + xfa1(nxta, j ) + xfa1(1, j ) )
        do i=2,nxpa-1
          entat1(i,j) = 0.25d0*(  xfa1(i-1,j-1) + xfa1(i,j-1)
     &                         + xfa1(i-1, j ) + xfa1(i, j ) )
        enddo
        entat1(nxpa,j) = entat1(1,j)
      enddo
!$OMP END DO NOWAIT

*     Do edge and corner values to conserve area integral
!$OMP SINGLE
      entat1(  1 ,  1 ) = 0.5d0*( xfa1(nxta,  1 ) + xfa1(1,  1 ) )
      entat1(  1 ,nypa) = 0.5d0*( xfa1(nxta,nyta) + xfa1(1,nyta) )
      do i=2,nxpa-1
        entat1(i,  1 ) = 0.5d0*( xfa1(i-1,  1 ) + xfa1(i,  1 ) )
        entat1(i,nypa) = 0.5d0*( xfa1(i-1,nyta) + xfa1(i,nyta) )
      enddo
      entat1(nxpa,  1 ) = entat1(  1 ,  1 )
      entat1(nxpa,nypa) = entat1(  1 ,nypa)
!$OMP END SINGLE


!$OMP END PARALLEL


*     MONITORING - extra section for cfraat, centat
      cfraat = cfrasm*atnorm
      centat = centsm*dxa*dya

*     Compute area integral of entrainment between layers 1 and 2
*     and also between layers 2 and 3! (SK: 07/22/2020)
      call xintp (xan(1), entat, nxpa, nypa)
      xan(1) = xan(1)*dxa*dya
      call xintp (xan(2), entat1, nxpa, nypa)
      xan(2) = xan(2)*dxa*dya

*     Compute boundary integrals of entrainment
*     -----------------------------------------
      ensums = 0.5d0*entat(1,  1 )
      ensumn = 0.5d0*entat(1,nypa)
      do i=2,nxpa-1
        ensums = ensums + entat(i,  1 )
        ensumn = ensumn + entat(i,nypa)
      enddo
      ensums = ensums + 0.5d0*entat(nxpa,  1 )
      ensumn = ensumn + 0.5d0*entat(nxpa,nypa)
      enisat(1) = dxa*ensums
      eninat(1) = dxa*ensumn

      ensums = 0.5d0*entat1(1,  1 )
      ensumn = 0.5d0*entat1(1,nypa)
      do i=2,nxpa-1
        ensums = ensums + entat1(i,  1 )
        ensumn = ensumn + entat1(i,nypa)
      enddo
      ensums = ensums + 0.5d0*entat1(nxpa,  1 )
      ensumn = ensumn + 0.5d0*entat1(nxpa,nypa)
      enisat(2) = dxa*ensums
      eninat(2) = dxa*ensumn

c      print*,'aml ends here'

      END SUBROUTINE aml

#endif /* not ocean_only */
c
c***********************************************************************
c
#ifndef ocean_only

      SUBROUTINE amladf (tmrhs, hmrhs, pa)
*
*     Computes (in tmrhs and hmrhs) the advective and diffusive
*     contributions to the RHS of the evolution equations for the
*     atmospheric mixed layer temperature and thickness.
*     Uses time-lagged mixed layer quantities.
*     pa is the pressure in the lowest atmospheric layer (layer 1).
*     Temperature equation contains both Del-sqd and Del-4th terms.
*     Diffusion implemented using dummy points west & east in del2t.
*     Apply no-flux boundary condition to Del-4th diffusion
*     across solid boundaries, equivalent to T'_{nnn} = 0
*     Version with N.E. corner bugfix included
*     (thanks to Tim Leslie for spotting that!).

c
c
c
c This is a version for the model with the constant atmospheric mixed-layer
c thickness. The thickness variable in this version is replaced with
c the mixed-layer specific humidity (hmixa,hmixam). The advective and
c diffusive moisture terms are only computed if the flag "moist" is activated.
c BCs for humidity mimic those for temperature: open boundary conditions for
c advection and no-flux conditions for diffusion are used. The hyper-viscosity
c in humidity equation is the same as for temperature. Note that different
c advection velocities are used in humidity and temperature equations when
c the temperature-dependent mixed-layer winds are applied (flag temp_fdbck).
c In particular, the temperature-dependent wind component is excluded  from the
c temperature equations, since its temperature advection is identically zero.
c
c   Sergey Kravtsov (10/09/2020)
c
c=======================================================================



*     Modules
      USE parameters, ONLY : nxta, nyta, nxpa, nypa
      USE atconst, ONLY : hdxam1, dxam2, rdxaf0
      USE intrfac, ONLY : ast, astm, uekat, vekat, 
     &                    uekat1, vekat1, hmixa,
     &                    hmixam, at2d, at4d

      IMPLICIT NONE

*     Subroutine arguments
      double precision, INTENT(OUT) :: tmrhs(nxta,nyta),hmrhs(nxta,nyta)
      double precision, INTENT(IN) :: pa(nxpa,nypa)
*
*     Local variables
      integer i,j
      double precision um,tm,hm,up,tp,hp,xadvt,xadvh,vm,vp,yadvt,
     &                 yadvh,d2tfac,d4tfac,hmdfac,del2t(0:nxta+1,nyta),
     &                 del2h(0:nxta+1,nyta),um1,up1,vm1,vp1

      d2tfac = at2d*dxam2
      d4tfac = at4d*dxam2**2
c      hmdfac = ahmd*dxam2
*
*     C-grid advection scheme, second order accurate
*     The Ekman components uekat1, vekat1 of the horizontal
*     flow have been previously computed in xforc

*     Inner points + meridional boundaries
*     ------------------------------------
!$OMP PARALLEL DEFAULT (NONE)
!$OMP&         PRIVATE (i,j,um,vm,tm,hm,up,vp,tp,hp,xadvt,xadvh,yadvt,
!$OMP&             yadvh,um1,vm1,up1,vp1)
!$OMP&         SHARED  (rdxaf0,uekat1,vekat1,pa,hdxam1,ast,hmixa,
!$OMP&                  tmrhs,del2t,del2h,astm,hmrhs,hmdfac,hmixam,
!$OMP&                  uekat,vekat)

!$OMP DO SCHEDULE (STATIC)
      do j=2,nyta-1

*       Western boundary (i=1)
*       ----------------
*       Zonally cyclic atmosphere
        um = -rdxaf0*( pa(1,j+1) - pa(1,j) ) + uekat1(1,j)
        tm = ast(1,j) + ast(nxta,j)
# ifdef moist
        um1 = -rdxaf0*( pa(1,j+1) - pa(1,j) ) + uekat(1,j)
        hm = hmixa(1,j) + hmixa(nxta,j)
# endif
        up = -rdxaf0*( pa(2,j+1) - pa(2,j) ) + uekat1(2,j)
        tp = ast(1,j) + ast(2,j)
# ifdef moist
        up1 = -rdxaf0*( pa(2,j+1) - pa(2,j) ) + uekat(2,j)
        hp = hmixa(1,j) + hmixa(2,j)
# endif
        xadvt = hdxam1*( up*tp - um*tm )
# ifdef moist
        xadvh = hdxam1*( up1*hp - um1*hm )
# endif
        vm =  rdxaf0*( pa(2, j ) - pa(1, j ) ) + vekat1(1, j )
        vp =  rdxaf0*( pa(2,j+1) - pa(1,j+1) ) + vekat1(1,j+1)
        yadvt = hdxam1*( vp*( ast(1,j+1)+ast(1, j ) )
     &                  -vm*( ast(1, j )+ast(1,j-1) ) )
# ifdef moist
        vm1 =  rdxaf0*( pa(2, j ) - pa(1, j ) ) + vekat(1, j )
        vp1 =  rdxaf0*( pa(2,j+1) - pa(1,j+1) ) + vekat(1,j+1)
        yadvh = hdxam1*( vp1*( hmixa(1,j+1)+hmixa(1, j ) )
     &                  -vm1*( hmixa(1, j )+hmixa(1,j-1) ) )
# endif
        tmrhs( 1 ,j) = -( xadvt + yadvt )
        del2t( 1 ,j) =   astm( 1 ,j-1) + astm(nxta,j) + astm(2,j)
     &                 + astm( 1 ,j+1) - 4.0d0*astm(1,j)
# ifdef moist
        hmrhs( 1 ,j) = -( xadvh + yadvh )
        del2h( 1 ,j) =   hmixam( 1 ,j-1) + hmixam(nxta,j) + hmixam(2,j)
     &                 + hmixam( 1 ,j+1) - 4.0d0*hmixam(1,j) 
# endif

*       Inner points (standard case)
*       ----------------------------
        do i=2,nxta-1
          um = up
          tm = tp
# ifdef moist
          um1=up1
          hm = hp
# endif
          up = -rdxaf0*( pa(i+1,j+1) - pa(i+1,j) ) + uekat1(i+1,j)
          tp = ast(i,j) + ast(i+1,j)
# ifdef moist
          up1 = -rdxaf0*( pa(i+1,j+1) - pa(i+1,j) ) + uekat(i+1,j)
          hp = hmixa(i,j) + hmixa(i+1,j)
# endif
          xadvt = hdxam1*( up*tp - um*tm )
# ifdef moist
          xadvh = hdxam1*( up1*hp - um1*hm )
# endif
          vm =  rdxaf0*( pa(i+1, j ) - pa(i, j ) ) + vekat1(i, j )
          vp =  rdxaf0*( pa(i+1,j+1) - pa(i,j+1) ) + vekat1(i,j+1)
          yadvt = hdxam1*( vp*(ast(i,j+1)+ast(i, j ))
     &                    -vm*(ast(i, j )+ast(i,j-1)) )
# ifdef moist
          vm1 =  rdxaf0*( pa(i+1, j ) - pa(i, j ) ) + vekat(i, j )
          vp1 =  rdxaf0*( pa(i+1,j+1) - pa(i,j+1) ) + vekat(i,j+1)
          yadvh = hdxam1*( vp1*(hmixa(i,j+1)+hmixa(i, j ))
     &                    -vm1*(hmixa(i, j )+hmixa(i,j-1)) )
# endif
          tmrhs(i,j) = -( xadvt + yadvt )
          del2t(i,j) =   astm( i ,j-1) + astm(i-1,j) + astm(i+1,j)
     &                 + astm( i ,j+1) - 4.0d0*astm(i,j)
# ifdef moist
          hmrhs(i,j) = -( xadvh + yadvh ) 
          del2h(i,j) =  hmixam( i ,j-1) + hmixam(i-1,j) + hmixam(i+1,j)
     &                + hmixam( i ,j+1) - 4.0d0*hmixam(i,j)
# endif
        enddo

*       Eastern boundary (i=nxta)
*       ----------------
*       Zonally cyclic atmosphere
        um = up
        tm = tp
# ifdef moist
        um1 = up1
        hm = hp
# endif
        up = -rdxaf0*( pa(nxpa,j+1) - pa(nxpa,j) ) + uekat1(nxpa,j)
        tp = ast(1,j) + ast(nxta,j)
# ifdef moist
        up1 = -rdxaf0*( pa(nxpa,j+1) - pa(nxpa,j) ) + uekat(nxpa,j)
        hp = hmixa(1,j) + hmixa(nxta,j)
# endif
        xadvt = hdxam1*( up*tp - um*tm )
# ifdef moist
        xadvh = hdxam1*( up1*hp - um1*hm )
# endif
        vm =  rdxaf0*( pa(nxpa, j ) - pa(nxta, j ) ) + vekat1(nxta, j )
        vp =  rdxaf0*( pa(nxpa,j+1) - pa(nxta,j+1) ) + vekat1(nxta,j+1)
        yadvt = hdxam1*( vp*( ast(nxta,j+1)+ast(nxta, j ) )
     &                  -vm*( ast(nxta, j )+ast(nxta,j-1) ) )
# ifdef moist
        vm1 = rdxaf0*( pa(nxpa, j ) - pa(nxta, j ) ) + vekat(nxta, j )
        vp1 = rdxaf0*( pa(nxpa,j+1) - pa(nxta,j+1) ) + vekat(nxta,j+1)
        yadvh = hdxam1*( vp1*( hmixa(nxta,j+1)+hmixa(nxta, j ) )
     &                  -vm1*( hmixa(nxta, j )+hmixa(nxta,j-1) ) )
# endif
        tmrhs(nxta,j) = -( xadvt + yadvt )
        del2t(nxta,j) =   astm(nxta,j-1) + astm(nxta-1,j) + astm(1,j)
     &                  + astm(nxta,j+1) - 4.0d0*astm(nxta,j)
# ifdef moist
        hmrhs(nxta,j) = -( xadvh + yadvh ) 
        del2h(nxta,j)= hmixam(nxta,j-1) + hmixam(nxta-1,j) +hmixam(1,j)
     &              + hmixam(nxta,j+1) - 4.0d0*hmixam(nxta,j) 
# endif

*       Fill dummy points of del2t array
*       --------------------------------
*       Zonally cyclic atmosphere
        del2t(   0  ,j) = del2t(nxta,j)
        del2t(nxta+1,j) = del2t(  1 ,j)
# ifdef moist
        del2h(   0  ,j) = del2h(nxta,j)
        del2h(nxta+1,j) = del2h(  1 ,j)
# endif

      enddo
!$OMP END DO NOWAIT

*     Zonal boundaries
*     ================
*     (N.B. can't do in recurrence form as above because //ising on i)
*     No normal heat flux in temperature equation
*     Normal flux allowed in thickness equation;
*     thickness set to hmat outside boundary.
*     N.B. P constant along boundaries, so Px vanishes.
*     Inner points (standard case)
*     ----------------------------
!$OMP DO SCHEDULE (STATIC)
      do i=2,nxta-1

*       Southern boundary (j=1)
        um = -rdxaf0*( pa( i ,2) - pa( i ,1) ) + uekat1( i ,1)
        up = -rdxaf0*( pa(i+1,2) - pa(i+1,1) ) + uekat1(i+1,1)
        xadvt = hdxam1*( up*( ast(i+1,1)+ast( i ,1) )
     &                  -um*( ast( i ,1)+ast(i-1,1) ) )
# ifdef moist
        um1 = -rdxaf0*( pa( i ,2) - pa( i ,1) ) + uekat( i ,1)
        up1 = -rdxaf0*( pa(i+1,2) - pa(i+1,1) ) + uekat(i+1,1)
        xadvh = hdxam1*( up1*( hmixa(i+1,1)+hmixa( i ,1) )
     &                  -um1*( hmixa( i ,1)+hmixa(i-1,1) ) )
# endif
        vm =  vekat1(i,1)
        vp =  rdxaf0*( pa(i+1,2) - pa(i,2) ) + vekat1(i,2)
        yadvt = hdxam1*vp*(ast(i,2)+ast(i, 1 ))
     &        -hdxam1*vm*(ast(i,1)+ast(i, 1 ))
# ifdef moist
        vm1 =  vekat(i,1)
        vp1 =  rdxaf0*( pa(i+1,2) - pa(i,2) ) + vekat(i,2)
        yadvh = hdxam1*( vp1*(hmixa(i,2)+hmixa(i,1))
     &                  -vm1*(hmixa(i,1)+hmixa(i,1)) )
# endif
        tmrhs(i,  1  ) = -( xadvt + yadvt )
        del2t(i,  1  ) =   astm(i-1,1) + astm(i+1,1)
     &                   + astm(i,2) - 3.0d0*astm(i,1)
# ifdef moist
        hmrhs(i,  1  ) = -( xadvh + yadvh )
        del2h(i,  1  ) =   hmixam(i-1,1) + hmixam(i+1,1)
     &                   + hmixam(i,2) - 3.0d0*hmixam(i,1)
# endif

*       Northern boundary (j=nyta)
        um = -rdxaf0*( pa( i ,nypa) - pa( i ,nyta) ) + uekat1( i ,nyta)
        up = -rdxaf0*( pa(i+1,nypa) - pa(i+1,nyta) ) + uekat1(i+1,nyta)
        xadvt = hdxam1*( up*( ast(i+1,nyta) + ast( i ,nyta) )
     &                  -um*( ast( i ,nyta) + ast(i-1,nyta) ) )
# ifdef moist
        um1 = -rdxaf0*( pa( i ,nypa) - pa( i ,nyta) ) + uekat( i ,nyta)
        up1 = -rdxaf0*( pa(i+1,nypa) - pa(i+1,nyta) ) + uekat(i+1,nyta)
        xadvh = hdxam1*( up1*( hmixa(i+1,nyta)+hmixa( i ,nyta) )
     &                  -um1*( hmixa( i ,nyta)+hmixa(i-1,nyta) ) )
# endif
        vm =  rdxaf0*( pa(i+1,nyta) - pa(i,nyta) ) + vekat1(i,nyta)
        vp =  vekat1(i,nyta+1)
        yadvt = hdxam1*( -vm*( ast(i, nyta )+ast(i,nyta-1) ) )
     &        +hdxam1*( vp*( ast(i, nyta )+ast(i,nyta) ) )
# ifdef moist
        vm1 =  rdxaf0*( pa(i+1,nyta) - pa(i,nyta) ) + vekat(i,nyta)
        vp1 =  vekat(i,nyta+1)
        yadvh = hdxam1*( vp1*( hmixa(i,nyta)+hmixa(i,nyta) )
     &                  -vm1*( hmixa(i,nyta)+hmixa(i,nyta-1) ) )
# endif
        tmrhs(i,nyta) = -( xadvt + yadvt )
        del2t(i,nyta) =   astm(i,nyta-1) + astm(i-1,nyta)
     &                  + astm(i+1,nyta) - 3.0d0*astm(i,nyta)
# ifdef moist
        hmrhs(i,nyta) = -( xadvh + yadvh )
        del2h(i,nyta) =   hmixam(i,nyta-1) + hmixam(i-1,nyta)
     &                  + hmixam(i+1,nyta) - 3.0d0*hmixam(i,nyta)
# endif

      enddo
!$OMP END DO NOWAIT

!$OMP END PARALLEL

*     Corner points, in memory order
*     ==============================
*     No normal flux for temperature and humidity
*     Periodic in other direction.
c-----------------------------------------------------------------
c
c
c    Open boundary option added (with advection of heat through the 
c    boundaries allowed), here and above
c
c     Sergey Kravtsov, 07/29/2020
c
c-----------------------------------------------------------------

*     Southwestern (i=1, j=1)
*     ------------
      um = -rdxaf0*( pa(1,2) - pa(1,1) ) + uekat1(1,1)
      up = -rdxaf0*( pa(2,2) - pa(2,1) ) + uekat1(2,1)
      xadvt = hdxam1*( up*( ast(2,1)+ast(  1 ,1) )
     &                -um*( ast(1,1)+ast(nxta,1) ) )
# ifdef moist
      um1 = -rdxaf0*( pa(1,2) - pa(1,1) ) + uekat(1,1)
      up1 = -rdxaf0*( pa(2,2) - pa(2,1) ) + uekat(2,1)
      xadvh = hdxam1*( up1*( hmixa(2,1)+hmixa(  1 ,1) )
     &                -um1*( hmixa(1,1)+hmixa(nxta,1) ) )
# endif
      vm =  vekat1(1,1)
      vp =  rdxaf0*( pa(2,2) - pa(1,2) ) + vekat1(1,2)
      yadvt = hdxam1*vp*( ast(1,2)+ast(1,1) )
     &       -hdxam1*vm*( ast(1,1)+ast(1,1) )
# ifdef moist
      vm1 =  vekat(1,1)
      vp1 =  rdxaf0*( pa(2,2) - pa(1,2) ) + vekat(1,2)
      yadvh = hdxam1*( vp1*( hmixa(1,2)+hmixa(1,1) )
     &                -vm1*( hmixa(1,1)+hmixa(1,1) ) )
# endif
      tmrhs( 1 , 1 ) = -( xadvt + yadvt )
      del2t( 1 , 1 ) =   astm(nxta,1) + astm(2,1)
     &                 + astm(1,2) - 3.0d0*astm(1,1)
      del2t(nxta+1,1) = del2t(1,1)
# ifdef moist
      hmrhs( 1 , 1 ) = -( xadvh + yadvh )
      del2h( 1 , 1 ) =   hmixam(nxta,1) + hmixa(2,1)
     &                 + hmixam(1,2) - 3.0d0*hmixam(1,1)
      del2h(nxta+1,1) = del2h(1,1)
# endif

*     Southeastern (i=nxta, j=1)
*     ------------
      um = -rdxaf0*( pa(nxta,2) - pa(nxta,1) ) + uekat1(nxta,1)
      up = -rdxaf0*( pa(nxpa,2) - pa(nxpa,1) ) + uekat1(nxpa,1)
      xadvt = hdxam1*( up*( ast(  1 ,1)+ast( nxta ,1) )
     &                -um*( ast(nxta,1)+ast(nxta-1,1) ) )
# ifdef moist
      um1 = -rdxaf0*( pa(nxta,2) - pa(nxta,1) ) + uekat(nxta,1)
      up1 = -rdxaf0*( pa(nxpa,2) - pa(nxpa,1) ) + uekat(nxpa,1)
      xadvh = hdxam1*( up1*( hmixa(  1 ,1)+hmixa( nxta ,1) )
     &                -um1*( hmixa(nxta,1)+hmixa(nxta-1,1) ) )
# endif
      vm =  vekat1(nxta,1)
      vp =  rdxaf0*( pa(nxpa,2) - pa(nxta,2) ) + vekat1(nxta,2)
      yadvt = hdxam1*vp*( ast(nxta,2)+ast(nxta, 1 ) )
     &       -hdxam1*vm*( ast(nxta,1)+ast(nxta, 1 ) )
# ifdef moist
      vm1 =  vekat(nxta,1)
      vp1 =  rdxaf0*( pa(nxpa,2) - pa(nxta,2) ) + vekat(nxta,2)
      yadvh = hdxam1*( vp1*(hmixa(nxta,2)+hmixa(nxta,1))
     &                -vm1*(hmixa(nxta,1)+hmixa(nxta,1)) )
# endif
      tmrhs(nxta,1) = -( xadvt + yadvt )
      del2t(nxta,1) =   astm(nxta-1,1) + astm(1,1)
     &                + astm(nxta,2) - 3.0d0*astm(nxta,1)
      del2t(0,1) = del2t(nxta,1)
# ifdef moist
      hmrhs(nxta,1) = -( xadvh + yadvh )
      del2h(nxta,1) =   hmixam(nxta-1,1) + hmixam(1,1)
     &                + hmixam(nxta,2) - 3.0d0*hmixam(nxta,1)
      del2h(0,1) = del2h(nxta,1)
#  endif

*     Northwestern (i=1, j=nyta)
*     ------------
      um = -rdxaf0*( pa(1,nypa) - pa(1,nyta) ) + uekat1(1,nyta)
      up = -rdxaf0*( pa(2,nypa) - pa(2,nyta) ) + uekat1(2,nyta)
      xadvt = hdxam1*( up*( ast( 2 ,nyta)+ast(  1 ,nyta) )
     &                -um*( ast( 1 ,nyta)+ast(nxta,nyta) ) )
# ifdef moist
      um1 = -rdxaf0*( pa(1,nypa) - pa(1,nyta) ) + uekat(1,nyta)
      up1 = -rdxaf0*( pa(2,nypa) - pa(2,nyta) ) + uekat(2,nyta)
      xadvh = hdxam1*( up1*( hmixa( 2 ,nyta)+hmixa(  1 ,nyta) )
     &                -um1*( hmixa( 1 ,nyta)+hmixa(nxta,nyta) ) )
# endif
      vm =  rdxaf0*( pa(2,nyta) - pa(1,nyta) ) + vekat1(1,nyta)
      vp =  vekat1(1,nyta+1)
      yadvt = hdxam1*(-vm*( ast(1, nyta)+ast(1,nyta-1) ) )
     &       +hdxam1*( vp*( ast(1, nyta)+ast(1,nyta) ) )
# ifdef moist
      vm1 =  rdxaf0*( pa(2,nyta) - pa(1,nyta) ) + vekat(1,nyta)
      vp1 =  vekat(1,nyta+1)
      yadvh = hdxam1*( vp1*( hmixa(1,nyta)+hmixa(1,nyta) )
     &                -vm1*( hmixa(1,nyta)+hmixa(1,nyta-1) ) )
# endif
      tmrhs( 1 ,nyta) = -( xadvt + yadvt )
      del2t( 1 ,nyta) =   astm(1,nyta-1) + astm(nxta,nyta)
     &                  + astm(2,nyta) - 3.0d0*astm(1,nyta)
      del2t(nxta+1,nyta) = del2t( 1 ,nyta)
# ifdef moist
      hmrhs( 1 ,nyta) = -( xadvh + yadvh )
      del2h( 1 ,nyta) =   hmixam(1,nyta-1) + hmixam(nxta,nyta)
     &                  + hmixam(2,nyta) - 3.0d0*hmixam(1,nyta)
      del2h(nxta+1,nyta) = del2h( 1 ,nyta)
# endif

*     Northeastern (i=nxta, j=nyta)
*     ------------
      um = -rdxaf0*( pa(nxta,nypa) - pa(nxta,nyta) ) + uekat1(nxta,nyta)
      up = -rdxaf0*( pa(nxpa,nypa) - pa(nxpa,nyta) ) + uekat1(nxpa,nyta)
      xadvt = hdxam1*( up*( ast(  1 ,nyta)+ast( nxta ,nyta) )
     &                -um*( ast(nxta,nyta)+ast(nxta-1,nyta) ) )
# ifdef moist
      um1 = -rdxaf0*( pa(nxta,nypa) - pa(nxta,nyta) ) + uekat(nxta,nyta)
      up1 = -rdxaf0*( pa(nxpa,nypa) - pa(nxpa,nyta) ) + uekat(nxpa,nyta)
      xadvh = hdxam1*( up1*( hmixa(  1 ,nyta)+hmixa( nxta ,nyta) )
     &                -um1*( hmixa(nxta,nyta)+hmixa(nxta-1,nyta) ) )
# endif
      vm =  rdxaf0*( pa(nxpa,nyta) - pa(nxta,nyta) ) + vekat1(nxta,nyta)
      vp =  vekat1(nxta,nyta+1 )
      yadvt = hdxam1*( -vm*( ast(nxta,nyta)+ast(nxta,nyta-1) ) )
     &       +hdxam1*(  vp*( ast(nxta,nyta)+ast(nxta,nyta) ) )
# ifdef moist
      vm1 =  rdxaf0*( pa(nxpa,nyta) - pa(nxta,nyta) ) + vekat(nxta,nyta)
      vp1 =  vekat(nxta,nyta+1 )
      yadvh = hdxam1*( vp1*( hmixa(nxta,nyta)+hmixa(nxta,nyta) )
     &                -vm1*( hmixa(nxta,nyta)+hmixa(nxta,nyta-1)) )
# endif
      tmrhs(nxta,nyta) = -( xadvt + yadvt )
      del2t(nxta,nyta) =   astm(nxta,nyta-1) + astm(nxta-1,nyta)
     &                   + astm(1,nyta) - 3.0d0*astm(nxta,nyta)
      del2t(0,nyta) = del2t(nxta,nyta)
# ifdef moist
      hmrhs(nxta,nyta) = -( xadvh + yadvh )
      del2h(nxta,nyta) =   hmixam(nxta,nyta-1) + hmixam(nxta-1,nyta)
     &                   + hmixam(1,nyta) - 3.0d0*hmixam(nxta,nyta)
      del2h(0,nyta) = del2h(nxta,nyta)
#  endif

*     Add Del-sqd and Del-4th terms to temperature evolution term
*     and Del-4th term to moisture equation
*     ===========================================================

!$OMP PARALLEL DEFAULT (NONE)
!$OMP&         PRIVATE (i,j)
!$OMP&         SHARED  (tmrhs,d2tfac,del2t,d4tfac)
# ifdef moist
!$OMP&         SHARED  (hmrhs,del2h)
#endif

*     Inner points in y; complete x range
*     -----------------------------------
!$OMP DO SCHEDULE (STATIC)
      do j=2,nyta-1
        do i=1,nxta
          tmrhs(i,j) = tmrhs( i ,j) + d2tfac*del2t( i ,j) - d4tfac*
     &                 (  del2t( i ,j-1) + del2t(i-1,j) + del2t(i+1,j)
     &                  + del2t( i ,j+1) - 4.0d0*del2t(i,j) )
# ifdef moist
          hmrhs(i,j) = hmrhs( i ,j) + d2tfac*del2h( i ,j) - d4tfac*
     &                 (  del2h( i ,j-1) + del2h(i-1,j) + del2h(i+1,j)
     &                  + del2h( i ,j+1) - 4.0d0*del2h(i,j) )
# endif
        enddo
      enddo
!$OMP END DO NOWAIT

*     Zonal boundaries
*     ----------------
!$OMP DO SCHEDULE (STATIC)
      do i=1,nxta
*       Southern boundary (j=1)
        tmrhs(i,  1 ) = tmrhs(i,  1 ) + d2tfac*del2t(i,  1 )
     &                  - d4tfac*(  del2t(i-1,1) + del2t(i+1,1)
     &                            + del2t(i,2) - 3.0d0*del2t(i,1) )
*       Northern boundary (j=nyta)
        tmrhs(i,nyta) = tmrhs(i,nyta) + d2tfac*del2t(i,nyta)
     &                  - d4tfac*(  del2t(i,nyta-1) + del2t(i-1,nyta)
     &                     + del2t(i+1,nyta) - 3.0d0*del2t(i,nyta) )
# ifdef moist
        hmrhs(i,  1 ) = hmrhs(i,  1 )  + d2tfac*del2h(i,  1 )
     &                  - d4tfac*(  del2h(i-1,1) + del2h(i+1,1)
     &                            + del2h(i,2) - 3.0d0*del2h(i,1) )
        hmrhs(i,nyta) = hmrhs(i,nyta) + d2tfac*del2h(i,nyta)
     &                  - d4tfac*(  del2h(i,nyta-1) + del2h(i-1,nyta)
     &                     + del2h(i+1,nyta) - 3.0d0*del2h(i,nyta) )
# endif
      enddo
!$OMP END DO NOWAIT

!$OMP END PARALLEL

      END SUBROUTINE amladf

#endif /* not ocean_only */
c
c***********************************************************************
c

# ifndef ocean_only
# ifdef moist

      SUBROUTINE imoadf (humrhs, pa, hum, humm)
*
*     Computes (in humrhs) the advective and diffusive
*     contributions to the RHS of the evolution equations for the
*     interior moisture.
*     Uses time-lagged interior moisture.
*     pa is the pressure in a given atmospheric layer (layer 1).
*     hum and humm are here moisture at a given and the previos
*     time step, respectively
*
*     otherwise completely analogous to amladf
*
*    (S. Kravtsov, 10/27/2020)
*

c
c=======================================================================



*     Modules
      USE parameters, ONLY : nxta, nyta, nxpa, nypa
      USE atconst, ONLY : hdxam1, dxam2, rdxaf0
      USE intrfac, ONLY : at2d,at4d

      IMPLICIT NONE

*     Subroutine arguments
      double precision, INTENT(OUT) :: humrhs(nxta,nyta)
      double precision, INTENT(IN) :: pa(nxpa,nypa),
     &                   hum(nxta,nyta),humm(nxta,nyta)
*
*     Local variables
      integer i,j
      double precision um,tm,up,tp,xadvt,vm,vp,yadvt,
     &              d2tfac,d4tfac,del2hum(0:nxta+1,nyta)

      d2tfac = at2d*dxam2
      d4tfac = at4d*dxam2**2

*
*     C-grid advection scheme, second order accurate

*     Inner points + meridional boundaries
*     ------------------------------------
!$OMP PARALLEL DEFAULT (NONE)
!$OMP&         PRIVATE (i,j,um,vm,tm,up,vp,tp,xadvt,yadvt)
!$OMP&         SHARED  (rdxaf0,pa,hum,hdxam1,
!$OMP&                  humrhs,del2hum,humm)

!$OMP DO SCHEDULE (STATIC)
      do j=2,nyta-1

*       Western boundary (i=1)
*       ----------------
*       Zonally cyclic atmosphere
        um = -rdxaf0*( pa(1,j+1) - pa(1,j) )
        tm = hum(1,j) + hum(nxta,j)
        up = -rdxaf0*( pa(2,j+1) - pa(2,j) )
        tp = hum(1,j) + hum(2,j)
        xadvt = hdxam1*( up*tp - um*tm )
        vm =  rdxaf0*( pa(2, j ) - pa(1, j ) )
        vp =  rdxaf0*( pa(2,j+1) - pa(1,j+1) )
        yadvt = hdxam1*( vp*( hum(1,j+1)+hum(1, j ) )
     &                  -vm*( hum(1, j )+hum(1,j-1) ) )
        humrhs( 1 ,j) = -( xadvt + yadvt )
        del2hum( 1 ,j) =   humm( 1 ,j-1) + humm(nxta,j) + humm(2,j)
     &                 + humm( 1 ,j+1) - 4.0d0*humm(1,j)

*       Inner points (standard case)
*       ----------------------------
        do i=2,nxta-1
          um = up
          tm = tp
          up = -rdxaf0*( pa(i+1,j+1) - pa(i+1,j) )
          tp = hum(i,j) + hum(i+1,j)
          xadvt = hdxam1*( up*tp - um*tm )
          vm =  rdxaf0*( pa(i+1, j ) - pa(i, j ) )
          vp =  rdxaf0*( pa(i+1,j+1) - pa(i,j+1) )
          yadvt = hdxam1*( vp*(hum(i,j+1)+hum(i, j ))
     &                    -vm*(hum(i, j )+hum(i,j-1)) )
          humrhs(i,j) = -( xadvt + yadvt )
          del2hum(i,j) =   humm( i ,j-1) + humm(i-1,j) + humm(i+1,j)
     &              + humm( i ,j+1) - 4.0d0*humm(i,j)
        enddo

*       Eastern boundary (i=nxta)
*       ----------------
*       Zonally cyclic atmosphere
        um = up
        tm = tp
        up = -rdxaf0*( pa(nxpa,j+1) - pa(nxpa,j) )
        tp = hum(1,j) + hum(nxta,j)
        xadvt = hdxam1*( up*tp - um*tm )
        vm =  rdxaf0*( pa(nxpa, j ) - pa(nxta, j ) )
        vp =  rdxaf0*( pa(nxpa,j+1) - pa(nxta,j+1) )
        yadvt = hdxam1*( vp*( hum(nxta,j+1)+hum(nxta, j ) )
     &                  -vm*( hum(nxta, j )+hum(nxta,j-1) ) )
        humrhs(nxta,j) = -( xadvt + yadvt )
        del2hum(nxta,j) =   humm(nxta,j-1) + humm(nxta-1,j) + humm(1,j)
     &                  + humm(nxta,j+1) - 4.0d0*humm(nxta,j)

*       Fill dummy points of del2hum array
*       --------------------------------
*       Zonally cyclic atmosphere
        del2hum(   0  ,j) = del2hum(nxta,j)
        del2hum(nxta+1,j) = del2hum(  1 ,j)

      enddo
!$OMP END DO NOWAIT

*     Zonal boundaries
*     ================
*     (N.B. can't do in recurrence form as above because //ising on i)
*     No normal heat flux in temperature equation
*     Normal flux allowed in thickness equation;
*     thickness set to hmat outside boundary.
*     N.B. P constant along boundaries, so Px vanishes.
*     Inner points (standard case)
*     ----------------------------
!$OMP DO SCHEDULE (STATIC)
      do i=2,nxta-1

*       Southern boundary (j=1)
        um = -rdxaf0*( pa( i ,2) - pa( i ,1) )
        up = -rdxaf0*( pa(i+1,2) - pa(i+1,1) )
        xadvt = hdxam1*( up*( hum(i+1,1)+hum( i ,1) )
     &                  -um*( hum( i ,1)+hum(i-1,1) ) )
        vm =  0.d0
        vp =  rdxaf0*( pa(i+1,2) - pa(i,2) )
        yadvt = hdxam1*vp*(hum(i,2)+hum(i, 1 ))
     &        -hdxam1*vm*(hum(i,1)+hum(i, 1 ))
        humrhs(i,  1  ) = -( xadvt + yadvt )
        del2hum(i,  1  ) =   humm(i-1,1) + humm(i+1,1)
     &                   + humm(i,2) - 3.0d0*humm(i,1)

*       Northern boundary (j=nyta)
        um = -rdxaf0*( pa( i ,nypa) - pa( i ,nyta) )
        up = -rdxaf0*( pa(i+1,nypa) - pa(i+1,nyta) )
        xadvt = hdxam1*( up*( hum(i+1,nyta) + hum( i ,nyta) )
     &                  -um*( hum( i ,nyta) + hum(i-1,nyta) ) )
        vm =  rdxaf0*( pa(i+1,nyta) - pa(i,nyta) )
        vp =  0.d0
        yadvt = hdxam1*( -vm*( hum(i, nyta )+hum(i,nyta-1) ) )
     &        +hdxam1*( vp*( hum(i, nyta )+hum(i,nyta) ) )
        humrhs(i,nyta) = -( xadvt + yadvt )
        del2hum(i,nyta) =   humm(i,nyta-1) + humm(i-1,nyta)
     &                  + humm(i+1,nyta) - 3.0d0*humm(i,nyta)

      enddo
!$OMP END DO NOWAIT

!$OMP END PARALLEL

*     Corner points, in memory order
*     ==============================
*     No normal flux in y
*     Periodic in other direction.

*     Southwestern (i=1, j=1)
*     ------------
      um = -rdxaf0*( pa(1,2) - pa(1,1) )
      up = -rdxaf0*( pa(2,2) - pa(2,1) )
      xadvt = hdxam1*( up*( hum(2,1)+hum(  1 ,1) )
     &                -um*( hum(1,1)+hum(nxta,1) ) )
      vm =  0.d0
      vp =  rdxaf0*( pa(2,2) - pa(1,2) )
      yadvt = hdxam1*vp*( hum(1,2)+hum(1,1) )
     &       -hdxam1*vm*( hum(1,1)+hum(1,1) )
      humrhs( 1 , 1 ) = -( xadvt + yadvt )
      del2hum( 1 , 1 ) =   humm(nxta,1) + humm(2,1)
     &                 + humm(1,2) - 3.0d0*humm(1,1)
      del2hum(nxta+1,1) = del2hum(1,1)
*     Southeastern (i=nxta, j=1)
*     ------------
      um = -rdxaf0*( pa(nxta,2) - pa(nxta,1) )
      up = -rdxaf0*( pa(nxpa,2) - pa(nxpa,1) )
      xadvt = hdxam1*( up*( hum(  1 ,1)+hum( nxta ,1) )
     &                -um*( hum(nxta,1)+hum(nxta-1,1) ) )
      vm =  0.d0
      vp =  rdxaf0*( pa(nxpa,2) - pa(nxta,2) )
      yadvt = hdxam1*vp*( hum(nxta,2)+hum(nxta, 1 ) )
     &       -hdxam1*vm*( hum(nxta,1)+hum(nxta, 1 ) )
      humrhs(nxta,1) = -( xadvt + yadvt )
      del2hum(nxta,1) =   humm(nxta-1,1) + humm(1,1)
     &                + humm(nxta,2) - 3.0d0*humm(nxta,1)
      del2hum(0,1) = del2hum(nxta,1)

*     Northwestern (i=1, j=nyta)
*     ------------
      um = -rdxaf0*( pa(1,nypa) - pa(1,nyta) )
      up = -rdxaf0*( pa(2,nypa) - pa(2,nyta) )
      xadvt = hdxam1*( up*( hum( 2 ,nyta)+hum(  1 ,nyta) )
     &                -um*( hum( 1 ,nyta)+hum(nxta,nyta) ) )
      vm =  rdxaf0*( pa(2,nyta) - pa(1,nyta) )
      vp =  0.d0
      yadvt = hdxam1*(-vm*( hum(1, nyta)+hum(1,nyta-1) ) )
     &       +hdxam1*( vp*( hum(1, nyta)+hum(1,nyta) ) )
      humrhs( 1 ,nyta) = -( xadvt + yadvt )
      del2hum( 1 ,nyta) =   humm(1,nyta-1) + humm(nxta,nyta)
     &                  + humm(2,nyta) - 3.0d0*humm(1,nyta)
      del2hum(nxta+1,nyta) = del2hum( 1 ,nyta)

*     Northehumern (i=nxta, j=nyta)
*     ------------
      um = -rdxaf0*( pa(nxta,nypa) - pa(nxta,nyta) )
      up = -rdxaf0*( pa(nxpa,nypa) - pa(nxpa,nyta) )
      xadvt = hdxam1*( up*( hum(  1 ,nyta)+hum( nxta ,nyta) )
     &                -um*( hum(nxta,nyta)+hum(nxta-1,nyta) ) )
      vm =  rdxaf0*( pa(nxpa,nyta) - pa(nxta,nyta) )
      vp =  0.d0
      yadvt = hdxam1*( -vm*( hum(nxta,nyta)+hum(nxta,nyta-1) ) )
     &       +hdxam1*(  vp*( hum(nxta,nyta)+hum(nxta,nyta) ) )
      humrhs(nxta,nyta) = -( xadvt + yadvt )
      del2hum(nxta,nyta) =   humm(nxta,nyta-1) + humm(nxta-1,nyta)
     &                   + humm(1,nyta) - 3.0d0*humm(nxta,nyta)
      del2hum(0,nyta) = del2hum(nxta,nyta)

*     Add Del-4th terms
*     ===========================================================

!$OMP PARALLEL DEFAULT (NONE)
!$OMP&         PRIVATE (i,j)
!$OMP&         SHARED  (humrhs,del2hum,d2tfac,d4tfac)
*     Inner points in y; complete x range
*     -----------------------------------
!$OMP DO SCHEDULE (STATIC)
      do j=2,nyta-1
        do i=1,nxta
          humrhs(i,j) = humrhs( i ,j)+d2tfac*del2hum( i ,j) - d4tfac*
     &       (  del2hum( i ,j-1) + del2hum(i-1,j) + del2hum(i+1,j)
     &             + del2hum( i ,j+1) - 4.0d0*del2hum(i,j) )
        enddo
      enddo
!$OMP END DO NOWAIT

*     Zonal boundaries
*     ----------------
!$OMP DO SCHEDULE (STATIC)
      do i=1,nxta
*       Southern boundary (j=1)
        humrhs(i,  1 ) = humrhs(i,  1 )+ d2tfac*del2hum(i,  1 )
     &             - d4tfac*(  del2hum(i-1,1) + del2hum(i+1,1)
     &                    + del2hum(i,2) - 3.0d0*del2hum(i,1) )
*       Northern boundary (j=nyta)
        humrhs(i,nyta) = humrhs(i,nyta)+ d2tfac*del2hum(i,nyta)
     &        - d4tfac*(  del2hum(i,nyta-1) + del2hum(i-1,nyta)
     &                 + del2hum(i+1,nyta) - 3.0d0*del2hum(i,nyta) )
      enddo
!$OMP END DO NOWAIT

!$OMP END PARALLEL

      END SUBROUTINE imoadf

#endif /* moist */
#endif /* not ocean_only */

c
c***********************************************************************
c

# ifdef moist

      SUBROUTINE trapin (fofz, nz, delz, vinteg)
*
*     Computes the extended trapezoidal rule approximation to
*     the integral of a tabulated function. fofz is a vector
*     of length nz containing the function tabulation (including
*     at the end points), delz is the tabulation interval,
*     and vinteg is the returned value of the integral.
*
*     This version modified to use Kahan summation (a.k.a.
*     compensated summation), in order to improve accuracy.
*     We use the variant of Kahan summation where the correction from the
*     last term is added on, in case any further benefit can be achieved.

      IMPLICIT NONE

*     Subroutine arguments
      integer nz
      double precision fofz(nz),delz,vinteg
*
*     Local variables
      integer i
      double precision sum,corr,yadd,sest
*
      sum = 0.5d0*fofz(1)
      corr = 0.0d0
      do i=2,nz-1
        yadd = fofz(i) - corr
        sest = sum + yadd
        corr = (sest - sum ) - yadd
        sum = sest
      enddo
      yadd = 0.5d0*fofz(nz) - corr
      sest = sum + yadd
*     Compute and apply 'last term' correction
      corr = ( sest - sum ) - yadd
      sum = sest - corr
      vinteg = delz*sum

      END SUBROUTINE trapin

#endif /* moist */
c
c***********************************************************************
c

      END MODULE amlsubs
c
c***********************************************************************
