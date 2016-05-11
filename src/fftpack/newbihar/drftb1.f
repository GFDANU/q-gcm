      subroutine drftb1 (n,c,ch,wa,ifac)

*     Subroutine arguments
      integer n, ifac(*)
      double precision c(*), ch(*), wa(*)

*     Local variables
      integer nf, na, l1, iw, k1, ip, l2, ido, idl1, ix2, ix3, ix4, i
c
      nf = ifac(2)
      na = 0
      l1 = 1
      iw = 1
      do k1=1,nf
         ip = ifac(k1+2)
         l2 = ip*l1
         ido = n/l2
         idl1 = ido*l1
c
         if (ip .eq. 4) then
            ix2 = iw+ido
            ix3 = ix2+ido
            if (na .ne. 0) then
               call dradb4 (ido,l1,ch,c,wa(iw),wa(ix2),wa(ix3))
              else
               call dradb4 (ido,l1,c,ch,wa(iw),wa(ix2),wa(ix3))
            endif
            na = 1-na
          else if (ip .eq. 2) then
            if (na .ne. 0) then
               call dradb2 (ido,l1,ch,c,wa(iw))
              else
               call dradb2 (ido,l1,c,ch,wa(iw))
            endif
            na = 1-na
          else if (ip .eq. 3) then
            ix2 = iw+ido
            if (na .ne. 0) then
               call dradb3 (ido,l1,ch,c,wa(iw),wa(ix2))
              else
               call dradb3 (ido,l1,c,ch,wa(iw),wa(ix2))
            endif
            na = 1-na
          else if (ip .eq. 5) then
            ix2 = iw+ido
            ix3 = ix2+ido
            ix4 = ix3+ido
            if (na .ne. 0) then
               call dradb5 (ido,l1,ch,c,wa(iw),wa(ix2),wa(ix3),wa(ix4))
              else
               call dradb5 (ido,l1,c,ch,wa(iw),wa(ix2),wa(ix3),wa(ix4))
            endif
            na = 1-na
          else
            if (na .ne. 0) then
               call dradbg (ido,ip,l1,idl1,ch,ch,ch,c,c,wa(iw))
              else
               call dradbg (ido,ip,l1,idl1,c,c,c,ch,ch,wa(iw))
            endif
            if (ido .eq. 1) na = 1-na
         endif
c
         l1 = l2
         iw = iw+(ip-1)*ido
      enddo
c
      if (na .ne. 0) then
         do i=1,n
            c(i) = ch(i)
         enddo
      endif
c
      return
      end
