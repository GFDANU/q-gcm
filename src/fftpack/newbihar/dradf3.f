      subroutine dradf3 (ido,l1,cc,ch,wa1,wa2)

*     Subroutine arguments
      integer ido, l1
      double precision cc(ido,l1,3), ch(ido,3,l1), wa1(*), wa2(*)

*     Local parameters
      double precision TAUR, TAUI
      PARAMETER ( TAUR = -0.5 d0, TAUI = 0.8660254037 84438647 d0 )

*     Local variables
      integer k, idp2, i
      double precision ci2, cr2, di2, di3, dr2, dr3, ti2, ti3, tr2, tr3
c
      do k=1,l1
         cr2 = cc(1,k,2)+cc(1,k,3)
         ch(1,1,k) = cc(1,k,1)+cr2
         ch(1,3,k) = TAUI*(cc(1,k,3)-cc(1,k,2))
         ch(ido,2,k) = cc(1,k,1)+TAUR*cr2
      enddo
c
      if (ido .eq. 1) return
      idp2 = ido+2
      do k=1,l1
         do i=3,ido,2
            dr2 = wa1(i-2)*cc(i-1,k,2) + wa1(i-1)*cc( i ,k,2)
            di2 = wa1(i-2)*cc( i ,k,2) - wa1(i-1)*cc(i-1,k,2)
            dr3 = wa2(i-2)*cc(i-1,k,3) + wa2(i-1)*cc( i ,k,3)
            di3 = wa2(i-2)*cc( i ,k,3) - wa2(i-1)*cc(i-1,k,3)
            cr2 = dr2 + dr3
            ci2 = di2 + di3
            ch(i-1,1,k) = cc(i-1,k,1) + cr2
            ch( i ,1,k) = cc( i ,k,1) + ci2
            tr2 = cc(i-1,k,1) + TAUR*cr2
            ti2 = cc( i ,k,1) + TAUR*ci2
            tr3 = TAUI*(di2-di3)
            ti3 = TAUI*(dr3-dr2)
            ch(i-1,3,k) = tr2+tr3
            ch(idp2-i-1,2,k) = tr2-tr3
            ch(i,3,k) = ti2+ti3
            ch(idp2-i,2,k) = ti3-ti2
         enddo
      enddo
c
      return
      end
