      subroutine dradf5 (ido,l1,cc,ch,wa1,wa2,wa3,wa4)

*     Subroutine arguments
      integer ido, l1
      double precision cc(ido,l1,5), ch(ido,5,l1), wa1(*), wa2(*),
     1  wa3(*), wa4(*)

*     Local parameters
      double precision TR11, TI11, TR12, TI12
      PARAMETER ( TR11 =  0.3090169943 74947424 d0,
     &            TI11 =  0.9510565162 95153572 d0,
     &            TR12 = -0.8090169943 74947424 d0,
     &            TI12 =  0.5877852522 92473129 d0 )

*     Local variables
      double precision ci2, ci3, ci4, ci5, cr2, cr3, cr4, cr5, di2,
     2  di3, di4, di5, dr2, dr3, dr4, dr5, ti2, ti3, ti4,
     3  ti5, tr2, tr3, tr4, tr5
      integer k, idp2, i
c
      do k=1,l1
         cr2 = cc(1,k,5) + cc(1,k,2)
         ci5 = cc(1,k,5) - cc(1,k,2)
         cr3 = cc(1,k,4) + cc(1,k,3)
         ci4 = cc(1,k,4) - cc(1,k,3)
         ch(1,1,k) = cc(1,k,1)+cr2+cr3
         ch(ido,2,k) = cc(1,k,1)+TR11*cr2+TR12*cr3
         ch(1,3,k) = TI11*ci5+TI12*ci4
         ch(ido,4,k) = cc(1,k,1)+TR12*cr2+TR11*cr3
         ch(1,5,k) = TI12*ci5-TI11*ci4
      enddo
c
      if (ido .eq. 1) return
      idp2 = ido+2
      do k=1,l1
         do i=3,ido,2
            dr2 = wa1(i-2)*cc(i-1,k,2)+wa1(i-1)*cc(i,k,2)
            di2 = wa1(i-2)*cc(i,k,2)-wa1(i-1)*cc(i-1,k,2)
            dr3 = wa2(i-2)*cc(i-1,k,3)+wa2(i-1)*cc(i,k,3)
            di3 = wa2(i-2)*cc(i,k,3)-wa2(i-1)*cc(i-1,k,3)
            dr4 = wa3(i-2)*cc(i-1,k,4)+wa3(i-1)*cc(i,k,4)
            di4 = wa3(i-2)*cc(i,k,4)-wa3(i-1)*cc(i-1,k,4)
            dr5 = wa4(i-2)*cc(i-1,k,5)+wa4(i-1)*cc(i,k,5)
            di5 = wa4(i-2)*cc(i,k,5)-wa4(i-1)*cc(i-1,k,5)
            cr2 = dr2 + dr5
            ci5 = dr5 - dr2
            cr5 = di2 - di5
            ci2 = di2 + di5
            cr3 = dr3 + dr4
            ci4 = dr4 - dr3
            cr4 = di3 - di4
            ci3 = di3 + di4
            ch(i-1,1,k) = cc(i-1,k,1)+cr2+cr3
            ch(i,1,k) = cc(i,k,1)+ci2+ci3
            tr2 = cc(i-1,k,1)+TR11*cr2+TR12*cr3
            ti2 = cc(i,k,1)+TR11*ci2+TR12*ci3
            tr3 = cc(i-1,k,1)+TR12*cr2+TR11*cr3
            ti3 = cc(i,k,1)+TR12*ci2+TR11*ci3
            tr5 = TI11*cr5 + TI12*cr4
            ti5 = TI11*ci5 + TI12*ci4
            tr4 = TI12*cr5 - TI11*cr4
            ti4 = TI12*ci5 - TI11*ci4
            ch(i-1,3,k) = tr2+tr5
            ch(idp2-i-1,2,k) = tr2-tr5
            ch(i,3,k) = ti2+ti5
            ch(idp2-i,2,k) = ti5-ti2
            ch(i-1,5,k) = tr3+tr4
            ch(idp2-i-1,4,k) = tr3-tr4
            ch(i,5,k) = ti3+ti4
            ch(idp2-i,4,k) = ti4-ti3
         enddo
      enddo
c
      return
      end
