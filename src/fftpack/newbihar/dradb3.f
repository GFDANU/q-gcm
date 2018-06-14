      subroutine dradb3 (ido,l1,cc,ch,wa1,wa2)

*     Subroutine arguments
      integer ido, l1
      double precision cc(ido,3,l1), ch(ido,l1,3), wa1(*), wa2(*)

*     Local parameters
      double precision TAUR, TAUI
      PARAMETER ( TAUR = -0.5 d0, TAUI = 0.8660254037 84438647 d0 )

*     Local variables
      integer k, idp2, i, ic
      double precision  ci2, ci3, cr2, cr3,
     1  di2, di3, dr2, dr3, ti2, tr2
c
      do k=1,l1
         tr2 = cc(ido,2,k)+cc(ido,2,k)
         cr2 = cc(1,1,k)+TAUR*tr2
         ch(1,k,1) = cc(1,1,k)+tr2
         ci3 = TAUI*(cc(1,3,k)+cc(1,3,k))
         ch(1,k,2) = cr2-ci3
         ch(1,k,3) = cr2+ci3
      enddo
c
      if (ido .eq. 1) return
      idp2 = ido+2
      do k=1,l1
         do i=3,ido,2
            ic = idp2-i
            tr2 = cc(i-1,3,k)+cc(ic-1,2,k)
            cr2 = cc(i-1,1,k)+TAUR*tr2
            ch(i-1,k,1) = cc(i-1,1,k)+tr2
            ti2 = cc(i,3,k)-cc(ic,2,k)
            ci2 = cc(i,1,k)+TAUR*ti2
            ch(i,k,1) = cc(i,1,k)+ti2
            cr3 = TAUI*(cc(i-1,3,k)-cc(ic-1,2,k))
            ci3 = TAUI*(cc(i,3,k)+cc(ic,2,k))
            dr2 = cr2 - ci3
            dr3 = cr2 + ci3
            di2 = ci2 + cr3
            di3 = ci2 - cr3
            ch(i-1,k,2) = wa1(i-2)*dr2-wa1(i-1)*di2
            ch(i,k,2) = wa1(i-2)*di2+wa1(i-1)*dr2
            ch(i-1,k,3) = wa2(i-2)*dr3-wa2(i-1)*di3
            ch(i,k,3) = wa2(i-2)*di3+wa2(i-1)*dr3
         enddo
      enddo
c
      return
      end
