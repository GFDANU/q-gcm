      subroutine drftf1 (n,c,ch,wa,ifac)

*     Subroutine arguments
      integer n, ifac(*)
      double precision c(*), ch(*), wa(*)

*     Local variables
      integer nf, na, l2, iw, k1, kh, ip, l1,
     &        ido, idl1, ix2, ix3, ix4, i
c
      nf = ifac(2)
      na = 1
      l2 = n
      iw = n
      do k1=1,nf
         kh = nf-k1
         ip = ifac(kh+3)
         l1 = l2/ip
         ido = n/l2
         idl1 = ido*l1
         iw = iw-(ip-1)*ido
         na = 1-na
c
         if (ip .eq. 4) then
            ix2 = iw+ido
            ix3 = ix2+ido
            if (na .ne. 0) then
               call dradf4 (ido,l1,ch,c,wa(iw),wa(ix2),wa(ix3))
              else
               call dradf4 (ido,l1,c,ch,wa(iw),wa(ix2),wa(ix3))
            endif
          else if (ip .eq. 2) then
            if (na .ne. 0) then
               call dradf2 (ido,l1,ch,c,wa(iw))
              else
               call dradf2 (ido,l1,c,ch,wa(iw))
            endif
          else if (ip .eq. 3) then
            ix2 = iw+ido
            if (na .ne. 0) then
               call dradf3 (ido,l1,ch,c,wa(iw),wa(ix2))
              else
               call dradf3 (ido,l1,c,ch,wa(iw),wa(ix2))
            endif
          else if (ip .eq. 5) then
            ix2 = iw+ido
            ix3 = ix2+ido
            ix4 = ix3+ido
            if (na .ne. 0) then
               call dradf5 (ido,l1,ch,c,wa(iw),wa(ix2),wa(ix3),wa(ix4))
              else
               call dradf5 (ido,l1,c,ch,wa(iw),wa(ix2),wa(ix3),wa(ix4))
            endif
          else
            if (ido .eq. 1) na = 1-na
            if (na .ne. 0) then
               call dradfg (ido,ip,l1,idl1,ch,ch,ch,c,c,wa(iw))
               na = 0
              else
               call dradfg (ido,ip,l1,idl1,c,c,c,ch,ch,wa(iw))
               na = 1
            endif
         endif
c
         l2 = l1
      enddo
c
      if (na .ne. 1) then
         do i=1,n
            c(i) = ch(i)
         enddo
      endif
c
      return
      end
