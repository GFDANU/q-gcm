      subroutine dradbg (ido,ip,l1,idl1,cc,c1,c2,ch,ch2,wa)

*     Subroutine arguments
      integer ido, ip, l1, idl1
      double precision cc(ido,ip,l1), c1(ido,l1,ip), c2(idl1,ip),
     1  ch(ido,l1,ip), ch2(idl1,ip), wa(*)

*     Local parameters
      double precision TPI
      PARAMETER ( TPI = 6.2831853071 7958648 d0 )

*     Local variables
      double precision ai1, ai2, ar1, ar1h, ar2,
     2  ar2h, arg, dc2, dcp, ds2, dsp
      integer idp2, nbd, ipp2, ipph, k, i, j,
     &        jc, j2, ic, l, lc, ik, is, idij
c
      arg = TPI/dble(ip)
      dcp = dcos(arg)
      dsp = dsin(arg)
      idp2 = ido+2
      nbd = (ido-1)/2
      ipp2 = ip+2
      ipph = (ip+1)/2
      if (ido .lt. l1) go to 103
      do k=1,l1
         do i=1,ido
            ch(i,k,1) = cc(i,1,k)
         enddo
      enddo
      go to 106
c
  103 do i=1,ido
         do k=1,l1
            ch(i,k,1) = cc(i,1,k)
         enddo
      enddo
c
  106 do j=2,ipph
         jc = ipp2-j
         j2 = j+j
         do k=1,l1
            ch(1,k,j) = cc(ido,j2-2,k)+cc(ido,j2-2,k)
            ch(1,k,jc) = cc(1,j2-1,k)+cc(1,j2-1,k)
         enddo
      enddo
c
      if (ido .eq. 1) go to 116
      if (nbd .lt. l1) go to 112
      do j=2,ipph
         jc = ipp2-j
         do k=1,l1
            do i=3,ido,2
               ic = idp2-i
               ch(i-1,k,j) = cc(i-1,2*j-1,k)+cc(ic-1,2*j-2,k)
               ch(i-1,k,jc) = cc(i-1,2*j-1,k)-cc(ic-1,2*j-2,k)
               ch(i,k,j) = cc(i,2*j-1,k)-cc(ic,2*j-2,k)
               ch(i,k,jc) = cc(i,2*j-1,k)+cc(ic,2*j-2,k)
            enddo
         enddo
      enddo
      go to 116
c
  112 do j=2,ipph
         jc = ipp2-j
         do i=3,ido,2
            ic = idp2-i
            do k=1,l1
               ch(i-1,k,j) = cc(i-1,2*j-1,k)+cc(ic-1,2*j-2,k)
               ch(i-1,k,jc) = cc(i-1,2*j-1,k)-cc(ic-1,2*j-2,k)
               ch(i,k,j) = cc(i,2*j-1,k)-cc(ic,2*j-2,k)
               ch(i,k,jc) = cc(i,2*j-1,k)+cc(ic,2*j-2,k)
            enddo
         enddo
      enddo
c
  116 ar1 = 1.0d0
      ai1 = 0.0d0
      do l=2,ipph
         lc = ipp2-l
         ar1h = dcp*ar1-dsp*ai1
         ai1 = dcp*ai1+dsp*ar1
         ar1 = ar1h
         do ik=1,idl1
            c2(ik,l) = ch2(ik,1)+ar1*ch2(ik,2)
            c2(ik,lc) = ai1*ch2(ik,ip)
         enddo
         dc2 = ar1
         ds2 = ai1
         ar2 = ar1
         ai2 = ai1
         do j=3,ipph
            jc = ipp2-j
            ar2h = dc2*ar2-ds2*ai2
            ai2 = dc2*ai2+ds2*ar2
            ar2 = ar2h
            do ik=1,idl1
               c2(ik,l) = c2(ik,l)+ar2*ch2(ik,j)
               c2(ik,lc) = c2(ik,lc)+ai2*ch2(ik,jc)
            enddo
         enddo
      enddo
c
      do j=2,ipph
         do ik=1,idl1
            ch2(ik,1) = ch2(ik,1)+ch2(ik,j)
         enddo
      enddo
c
      do j=2,ipph
         jc = ipp2-j
         do k=1,l1
            ch(1,k,j) = c1(1,k,j)-c1(1,k,jc)
            ch(1,k,jc) = c1(1,k,j)+c1(1,k,jc)
         enddo
      enddo
c
      if (ido .eq. 1) go to 132
      if (nbd .lt. l1) go to 128
      do j=2,ipph
         jc = ipp2-j
         do k=1,l1
            do i=3,ido,2
               ch(i-1,k,j) = c1(i-1,k,j)-c1(i,k,jc)
               ch(i-1,k,jc) = c1(i-1,k,j)+c1(i,k,jc)
               ch(i,k,j) = c1(i,k,j)+c1(i-1,k,jc)
               ch(i,k,jc) = c1(i,k,j)-c1(i-1,k,jc)
            enddo
         enddo
      enddo
      go to 132
c
  128 do j=2,ipph
         jc = ipp2-j
         do i=3,ido,2
            do k=1,l1
               ch(i-1,k,j) = c1(i-1,k,j)-c1(i,k,jc)
               ch(i-1,k,jc) = c1(i-1,k,j)+c1(i,k,jc)
               ch(i,k,j) = c1(i,k,j)+c1(i-1,k,jc)
               ch(i,k,jc) = c1(i,k,j)-c1(i-1,k,jc)
            enddo
         enddo
      enddo
  132 continue
c
      if (ido .eq. 1) return
      do ik=1,idl1
         c2(ik,1) = ch2(ik,1)
      enddo
c
      do j=2,ip
         do k=1,l1
            c1(1,k,j) = ch(1,k,j)
         enddo
      enddo
c
      if (nbd .gt. l1) go to 139
      is = -ido
      do j=2,ip
         is = is+ido
         idij = is
         do i=3,ido,2
            idij = idij+2
            do k=1,l1
               c1(i-1,k,j) = wa(idij-1)*ch(i-1,k,j)-wa(idij)*ch(i,k,j)
               c1(i,k,j) = wa(idij-1)*ch(i,k,j)+wa(idij)*ch(i-1,k,j)
            enddo
         enddo
      enddo
      go to 143
c
  139 is = -ido
      do j=2,ip
         is = is+ido
         do k=1,l1
            idij = is
            do i=3,ido,2
               idij = idij+2
               c1(i-1,k,j) = wa(idij-1)*ch(i-1,k,j)-wa(idij)*ch(i,k,j)
               c1(i,k,j) = wa(idij-1)*ch(i,k,j)+wa(idij)*ch(i-1,k,j)
            enddo
         enddo
      enddo
c
  143 return
      end
