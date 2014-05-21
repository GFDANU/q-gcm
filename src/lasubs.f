*----------------------------------------------------------------------*
*                     Include LAPACK/BLAS routines                     *
*----------------------------------------------------------------------*
*                                           ! equivalent NAG name
*     For eigmod:
*     -----------
      INCLUDE 'lapack/dgebal.f'             ! DGEBAL = F08NHF
      INCLUDE 'lapack/dgehrd.f'             ! DGEHRD = F08NEF
      INCLUDE 'lapack/dorghr.f'             ! DORGHR = F08NFF
      INCLUDE 'lapack/dhseqr.f'             ! DHSEQR = F08PEF
      INCLUDE 'lapack/dtrevc.f'             ! DTREVC = F08QKF
      INCLUDE 'lapack/dgebak.f'             ! DGEBAK = F08NJF
*
*     For homsol/ocinvq/radiat:
*     -------------------------
      INCLUDE 'lapack/dgetrf.f'             ! DGETRF = F07ADF
      INCLUDE 'lapack/dgetrs.f'             ! DGETRS = F07AEF
      INCLUDE 'lapack/dgerfs.f'             ! DGERFS = F07AHF
*
*     LAPACK routines called at level 1:
*     ----------------------------------
      INCLUDE 'lapack/ilaenv.f'             ! ILAENV =
      INCLUDE 'lapack/dlamch.f'             ! DLAMCH =
      INCLUDE 'lapack/dlahrd.f'             ! DLAHRD = F08NEY
      INCLUDE 'lapack/dlarfb.f'             ! DLARFB = F08AEY
      INCLUDE 'lapack/dgehd2.f'             ! DGEHD2 = F08NEZ
      INCLUDE 'lapack/dorgqr.f'             ! DORGQR = F08AFF
      INCLUDE 'lapack/dlaset.f'             ! DLASET =
      INCLUDE 'lapack/dlahqr.f'             ! DLAHQR = F08PEZ
      INCLUDE 'lapack/dlabad.f'             ! DLABAD =
      INCLUDE 'lapack/dlanhs.f'             ! DLANHS = F06RMF
      INCLUDE 'lapack/dlacpy.f'             ! DLACPY =
      INCLUDE 'lapack/dlapy2.f'             ! DLAPY2 =
      INCLUDE 'lapack/dlarfg.f'             ! DLARFG = F08AEV
      INCLUDE 'lapack/dlarfx.f'             ! DLARFX = F08PEW
      INCLUDE 'lapack/dlaln2.f'             ! DLALN2 = F08QHX
      INCLUDE 'lapack/dgetf2.f'             ! DGETF2 = F07ADZ
      INCLUDE 'lapack/dlaswp.f'             ! DLASWP = F07ADY
      INCLUDE 'lapack/dlacon.f'             ! DLACON =
      INCLUDE 'lapack/dlahr2.f'             ! DLAHR2 =
      INCLUDE 'lapack/dlaqr0.f'             ! DLAQR0 =
      INCLUDE 'lapack/dlacn2.f'             ! DLACN2 =
*
*     BLAS routines called at level 1:
*     --------------------------------
      INCLUDE 'lapack/blas/xerbla.f'        ! XERBLA =
      INCLUDE 'lapack/blas/lsame.f'         ! LSAME  =
      INCLUDE 'lapack/blas/dswap.f'         ! DSWAP  = F06EGF
      INCLUDE 'lapack/blas/idamax.f'        ! IDAMAX = F06JLF
      INCLUDE 'lapack/blas/dscal.f'         ! DSCAL  = F06EDF
      INCLUDE 'lapack/blas/dgemm.f'         ! DGEMM  = F06YAF
      INCLUDE 'lapack/blas/dcopy.f'         ! DCOPY  = F06EFF
      INCLUDE 'lapack/blas/dgemv.f'         ! DGEMV  = F06PAF
      INCLUDE 'lapack/blas/daxpy.f'         ! DAXPY  = F06ECF
      INCLUDE 'lapack/blas/ddot.f'          ! DDOT   = F06EAF
      INCLUDE 'lapack/blas/dtrsm.f'         ! DTRSM  = F06YJF
*
*     LAPACK routines called at level 2:
*     ----------------------------------
      INCLUDE 'lapack/ieeeck.f'             ! IEEECK =
      INCLUDE 'lapack/dorg2r.f'             ! DORG2R = F08AFZ
      INCLUDE 'lapack/dladiv.f'             ! DLADIV =
      INCLUDE 'lapack/dlarft.f'             ! DLARFT = F08AEX
      INCLUDE 'lapack/dlassq.f'             ! DLASSQ =
      INCLUDE 'lapack/dlarf.f'              ! DLARF  = F08AEW
      INCLUDE 'lapack/dlanv2.f'             ! DLANV2 = F08PEY
      INCLUDE 'lapack/iladlc.f'             ! ILADLC =
      INCLUDE 'lapack/iladlr.f'             ! ILADLR =
      INCLUDE 'lapack/iparmq.f'             ! IPARMQ =
      INCLUDE 'lapack/dlaqr3.f'             ! DLAQR3 =
      INCLUDE 'lapack/dlaqr4.f'             ! DLAQR4 =
      INCLUDE 'lapack/dlaqr5.f'             ! DLAQR5 =
*
*     BLAS routines called at level 2:
*     --------------------------------
      INCLUDE 'lapack/blas/drot.f'          ! DROT   = F06EPF
      INCLUDE 'lapack/blas/dtrmv.f'         ! DTRMV  = F06PFF
      INCLUDE 'lapack/blas/dnrm2.f'         ! DNRM2  = F06EJF
      INCLUDE 'lapack/blas/dger.f'          ! DGER   = F06PMF
      INCLUDE 'lapack/blas/dtrmm.f'         ! DTRMM  = F06YFF
      INCLUDE 'lapack/blas/dasum.f'         ! DASUM  = F06EKF
*
*     LAPACK routines called at level 3:
*     ----------------------------------
      INCLUDE 'lapack/dlaqr1.f'             ! DLAQR1 =
      INCLUDE 'lapack/dormhr.f'             ! DORMHR =
      INCLUDE 'lapack/dtrexc.f'             ! DTREXC =
      INCLUDE 'lapack/dlaqr2.f'             ! DLAQR2 =
*
*     LAPACK routines called at level 4:
*     ----------------------------------
      INCLUDE 'lapack/dlaexc.f'             ! DLAEXC =
      INCLUDE 'lapack/dormqr.f'             ! DORMQR =
*
*     LAPACK routines called at level 5:
*     ----------------------------------
      INCLUDE 'lapack/dlange.f'             ! DLANGE =
      INCLUDE 'lapack/dlasy2.f'             ! DLASY2 =
      INCLUDE 'lapack/dlartg.f'             ! DLARTG =
      INCLUDE 'lapack/dorm2r.f'             ! DORM2R =
*
*----------------------------------------------------------------------*
