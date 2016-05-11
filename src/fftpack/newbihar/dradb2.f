      subroutine dradb2 (ido,l1,cc,ch,wa1)

*     Subroutine arguments
      integer ido, l1
      double precision cc(ido,2,l1), ch(ido,l1,2), wa1(*)

*     Local variables
      integer k, idp2, i, ic
      double precision ti2, tr2
c
      do k=1,l1
         ch(1,k,1) = cc(1,1,k) + cc(ido,2,k)
         ch(1,k,2) = cc(1,1,k) - cc(ido,2,k)
      enddo
c
      if ( ido.ge.2 ) then
         if ( ido.eq.2 ) goto 105
         idp2 = ido+2
         do k=1,l1
            do i=3,ido,2
               ic = idp2-i
               ch(i-1,k,1) = cc(i-1,1,k)+cc(ic-1,2,k)
               tr2 = cc(i-1,1,k)-cc(ic-1,2,k)
               ch(i,k,1) = cc(i,1,k)-cc(ic,2,k)
               ti2 = cc(i,1,k)+cc(ic,2,k)
               ch(i-1,k,2) = wa1(i-2)*tr2-wa1(i-1)*ti2
               ch(i,k,2) = wa1(i-2)*ti2+wa1(i-1)*tr2
            enddo
         enddo
c
         if (mod(ido,2) .eq. 1) return
  105    do k=1,l1
            ch(ido,k,1) = cc(ido,1,k)+cc(ido,1,k)
            ch(ido,k,2) = -(cc(1,2,k)+cc(1,2,k))
         enddo
      endif
c
      return
      end
