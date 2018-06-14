      subroutine dradb5 (ido,l1,cc,ch,wa1,wa2,wa3,wa4)

*     Subroutine arguments
      integer ido, l1
      double precision cc(ido,5,l1), ch(ido,l1,5), wa1(*), wa2(*),
     1  wa3(*), wa4(*)

*     Local parameters
      double precision TR11, TI11, TR12, TI12
      PARAMETER ( TR11 =  0.3090169943 74947424 d0,
     &            TI11 =  0.9510565162 95153572 d0,
     &            TR12 = -0.8090169943 74947424 d0,
     &            TI12 =  0.5877852522 92473129 d0 )

*     Local variables
      double precision ci2, ci3, ci4, ci5, cr2, cr3, cr4, cr5,
     2  di2, di3, di4, di5, dr2, dr3, dr4, dr5, ti2, ti3,
     3  ti4, ti5, tr2, tr3, tr4, tr5
      integer k, idp2, i, ic
c
      do k=1,l1
         ti5 = cc(1,3,k)+cc(1,3,k)
         ti4 = cc(1,5,k)+cc(1,5,k)
         tr2 = cc(ido,2,k)+cc(ido,2,k)
         tr3 = cc(ido,4,k)+cc(ido,4,k)
         ch(1,k,1) = cc(1,1,k)+tr2+tr3
         cr2 = cc(1,1,k)+TR11*tr2+TR12*tr3
         cr3 = cc(1,1,k)+TR12*tr2+TR11*tr3
         ci5 = TI11*ti5+TI12*ti4
         ci4 = TI12*ti5-TI11*ti4
         ch(1,k,2) = cr2 - ci5
         ch(1,k,3) = cr3 - ci4
         ch(1,k,4) = cr3 + ci4
         ch(1,k,5) = cr2 + ci5
      enddo
      if (ido .eq. 1) return
c
      idp2 = ido+2
      do k=1,l1
         do i=3,ido,2
            ic = idp2-i
            ti5 = cc(i,3,k) + cc(ic,2,k)
            ti2 = cc(i,3,k) - cc(ic,2,k)
            ti4 = cc(i,5,k) + cc(ic,4,k)
            ti3 = cc(i,5,k) - cc(ic,4,k)
            tr5 = cc(i-1,3,k) - cc(ic-1,2,k)
            tr2 = cc(i-1,3,k) + cc(ic-1,2,k)
            tr4 = cc(i-1,5,k) - cc(ic-1,4,k)
            tr3 = cc(i-1,5,k) + cc(ic-1,4,k)
            ch(i-1,k,1) = cc(i-1,1,k)+tr2+tr3
            ch(i,k,1) = cc(i,1,k)+ti2+ti3
            cr2 = cc(i-1,1,k)+TR11*tr2+TR12*tr3
            ci2 = cc(i,1,k)+TR11*ti2+TR12*ti3
            cr3 = cc(i-1,1,k)+TR12*tr2+TR11*tr3
            ci3 = cc(i,1,k)+TR12*ti2+TR11*ti3
            cr5 = TI11*tr5 + TI12*tr4
            ci5 = TI11*ti5 + TI12*ti4
            cr4 = TI12*tr5 - TI11*tr4
            ci4 = TI12*ti5 - TI11*ti4
            dr3 = cr3 - ci4
            dr4 = cr3 + ci4
            di3 = ci3 + cr4
            di4 = ci3 - cr4
            dr5 = cr2 + ci5
            dr2 = cr2 - ci5
            di5 = ci2 - cr5
            di2 = ci2 + cr5
            ch(i-1,k,2) = wa1(i-2)*dr2-wa1(i-1)*di2
            ch(i,k,2) = wa1(i-2)*di2+wa1(i-1)*dr2
            ch(i-1,k,3) = wa2(i-2)*dr3-wa2(i-1)*di3
            ch(i,k,3) = wa2(i-2)*di3+wa2(i-1)*dr3
            ch(i-1,k,4) = wa3(i-2)*dr4-wa3(i-1)*di4
            ch(i,k,4) = wa3(i-2)*di4+wa3(i-1)*dr4
            ch(i-1,k,5) = wa4(i-2)*dr5-wa4(i-1)*di5
            ch(i,k,5) = wa4(i-2)*di5+wa4(i-1)*dr5
         enddo
      enddo
c
      return
      end
