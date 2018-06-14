      subroutine dradf4 (ido,l1,cc,ch,wa1,wa2,wa3)

*     Subroutine arguments
      integer ido, l1
      double precision cc(ido,l1,4), ch(ido,4,l1), wa1(*), wa2(*),
     1  wa3(*)

*     Local parameters
      double precision HSQT2
      PARAMETER ( HSQT2 = .7071067811 86547524 d0 )

*     Local variables
      integer k, idp2, i
      double precision ci2, ci3, ci4, cr2, cr3, cr4,
     2  ti1, ti2, ti3, ti4, tr1, tr2, tr3, tr4
c
      do k=1,l1
         tr1 = cc(1,k,2) + cc(1,k,4)
         tr2 = cc(1,k,1) + cc(1,k,3)
         ch( 1 ,1,k) = tr1 + tr2
         ch(ido,4,k) = tr2 - tr1
         ch(ido,2,k) = cc(1,k,1) - cc(1,k,3)
         ch( 1 ,3,k) = cc(1,k,4) - cc(1,k,2)
      enddo

*     Cases tested with most likely first
      if ( (ido-2).gt.0 ) then
         idp2 = ido+2
         do k=1,l1
            do i=3,ido,2
               cr2 = wa1(i-2)*cc(i-1,k,2)+wa1(i-1)*cc(i,k,2)
               ci2 = wa1(i-2)*cc(i,k,2)-wa1(i-1)*cc(i-1,k,2)
               cr3 = wa2(i-2)*cc(i-1,k,3)+wa2(i-1)*cc(i,k,3)
               ci3 = wa2(i-2)*cc(i,k,3)-wa2(i-1)*cc(i-1,k,3)
               cr4 = wa3(i-2)*cc(i-1,k,4)+wa3(i-1)*cc(i,k,4)
               ci4 = wa3(i-2)*cc(i,k,4)-wa3(i-1)*cc(i-1,k,4)
               tr1 = cr2 + cr4
               tr4 = cr4 - cr2
               ti1 = ci2 + ci4
               ti4 = ci2 - ci4
               ti2 = cc( i ,k,1) + ci3
               ti3 = cc( i ,k,1) - ci3
               tr2 = cc(i-1,k,1) + cr3
               tr3 = cc(i-1,k,1) - cr3
               ch(i-1,1,k) = tr1+tr2
               ch(idp2-i-1,4,k) = tr2-tr1
               ch(i,1,k) = ti1+ti2
               ch(idp2-i,4,k) = ti1-ti2
               ch(i-1,3,k) = ti4+tr3
               ch(idp2-i-1,2,k) = tr3-ti4
               ch(i,3,k) = tr4+ti3
               ch(idp2-i,2,k) = tr4-ti3
            enddo
         enddo
         if (mod(ido,2) .eq. 1) return
       else if ( (ido-2).lt.0 ) then
         return
      endif
c
      do k=1,l1
         ti1 = -HSQT2*( cc(ido,k,2)+cc(ido,k,4) )
         tr1 =  HSQT2*( cc(ido,k,2)-cc(ido,k,4) )
         ch(ido,1,k) = tr1 + cc(ido,k,1)
         ch(ido,3,k) = cc(ido,k,1) - tr1
         ch(1,2,k) = ti1 - cc(ido,k,3)
         ch(1,4,k) = ti1 + cc(ido,k,3)
      enddo
c
      return
      end
