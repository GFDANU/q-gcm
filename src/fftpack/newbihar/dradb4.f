      subroutine dradb4 (ido,l1,cc,ch,wa1,wa2,wa3)

*     Subroutine arguments
      integer ido, l1
      double precision cc(ido,4,l1), ch(ido,l1,4), wa1(*), wa2(*),
     1  wa3(*)

*     Local parameters
      double precision SQRT2
      PARAMETER ( SQRT2 = 1.414213562 37309505 d0 )

*     Local variables
      integer k, idp2, i, ic
      double precision ci2, ci3, ci4, cr2, cr3, cr4,
     2  ti1, ti2, ti3, ti4, tr1, tr2, tr3, tr4
c
      do k=1,l1
         tr1 = cc( 1 ,1,k) - cc(ido,4,k)
         tr2 = cc( 1 ,1,k) + cc(ido,4,k)
         tr3 = cc(ido,2,k) + cc(ido,2,k)
         tr4 = cc( 1 ,3,k) + cc( 1 ,3,k)
         ch(1,k,1) = tr2 + tr3
         ch(1,k,2) = tr1 - tr4
         ch(1,k,3) = tr2 - tr3
         ch(1,k,4) = tr1 + tr4
      enddo

*     Cases tested with most likely first
      if ( (ido-2).gt.0 ) then
         idp2 = ido+2
         do k=1,l1
            do i=3,ido,2
               ic = idp2-i
               ti1 = cc(i,1,k) + cc(ic,4,k)
               ti2 = cc(i,1,k) - cc(ic,4,k)
               ti3 = cc(i,3,k) - cc(ic,2,k)
               tr4 = cc(i,3,k) + cc(ic,2,k)
               tr1 = cc(i-1,1,k) - cc(ic-1,4,k)
               tr2 = cc(i-1,1,k) + cc(ic-1,4,k)
               ti4 = cc(i-1,3,k) - cc(ic-1,2,k)
               tr3 = cc(i-1,3,k) + cc(ic-1,2,k)
               ch(i-1,k,1) = tr2+tr3
               cr3 = tr2-tr3
               ch(i,k,1) = ti2+ti3
               ci3 = ti2 - ti3
               cr2 = tr1 - tr4
               cr4 = tr1 + tr4
               ci2 = ti1 + ti4
               ci4 = ti1 - ti4
               ch(i-1,k,2) = wa1(i-2)*cr2-wa1(i-1)*ci2
               ch(i,k,2) = wa1(i-2)*ci2+wa1(i-1)*cr2
               ch(i-1,k,3) = wa2(i-2)*cr3-wa2(i-1)*ci3
               ch(i,k,3) = wa2(i-2)*ci3+wa2(i-1)*cr3
               ch(i-1,k,4) = wa3(i-2)*cr4-wa3(i-1)*ci4
               ch(i,k,4) = wa3(i-2)*ci4+wa3(i-1)*cr4
            enddo
         enddo
         if (mod(ido,2) .eq. 1) return
       else if ( (ido-2).lt.0 ) then
         return
      endif
c
      continue
      do k=1,l1
         ti1 = cc(1,2,k)+cc(1,4,k)
         ti2 = cc(1,4,k)-cc(1,2,k)
         tr1 = cc(ido,1,k)-cc(ido,3,k)
         tr2 = cc(ido,1,k)+cc(ido,3,k)
         ch(ido,k,1) = tr2 + tr2
         ch(ido,k,2) = SQRT2*(tr1-ti1)
         ch(ido,k,3) = ti2 + ti2
         ch(ido,k,4) = -SQRT2*(tr1+ti1)
      enddo
c
      return
      end
