      subroutine dradfg (ido,ip,l1,idl1,cc,c1,c2,ch,ch2,wa)

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
      integer ipph, ipp2, idp2, nbd, ik, j, k,
     &        is, idij, i, jc, l, lc, j2, ic
c
      arg = TPI/dble(ip)
      dcp = dcos(arg)
      dsp = dsin(arg)
      ipph = (ip+1)/2
      ipp2 = ip+2
      idp2 = ido+2
      nbd = (ido-1)/2
      if (ido .eq. 1) go to 119
      do ik=1,idl1
         ch2(ik,1) = c2(ik,1)
      enddo
      do j=2,ip
         do k=1,l1
            ch(1,k,j) = c1(1,k,j)
         enddo
      enddo
c
      if (nbd .gt. l1) go to 107
      is = -ido
      do j=2,ip
         is = is+ido
         idij = is
         do i=3,ido,2
            idij = idij+2
            do k=1,l1
               ch(i-1,k,j) = wa(idij-1)*c1(i-1,k,j)+wa(idij)*c1(i,k,j)
               ch(i,k,j) = wa(idij-1)*c1(i,k,j)-wa(idij)*c1(i-1,k,j)
            enddo
         enddo
      enddo
      go to 111
c
  107 is = -ido
      do j=2,ip
         is = is+ido
         do k=1,l1
            idij = is
            do i=3,ido,2
               idij = idij+2
               ch(i-1,k,j) = wa(idij-1)*c1(i-1,k,j)+wa(idij)*c1(i,k,j)
               ch(i,k,j) = wa(idij-1)*c1(i,k,j)-wa(idij)*c1(i-1,k,j)
            enddo
         enddo
      enddo
c
  111 if (nbd .lt. l1) go to 115
      do j=2,ipph
         jc = ipp2-j
         do k=1,l1
            do i=3,ido,2
               c1(i-1,k,j) = ch(i-1,k,j)+ch(i-1,k,jc)
               c1(i-1,k,jc) = ch(i,k,j)-ch(i,k,jc)
               c1(i,k,j) = ch(i,k,j)+ch(i,k,jc)
               c1(i,k,jc) = ch(i-1,k,jc)-ch(i-1,k,j)
            enddo
         enddo
      enddo
      go to 121
c
  115 do j=2,ipph
         jc = ipp2-j
         do i=3,ido,2
            do k=1,l1
               c1(i-1,k,j) = ch(i-1,k,j)+ch(i-1,k,jc)
               c1(i-1,k,jc) = ch(i,k,j)-ch(i,k,jc)
               c1(i,k,j) = ch(i,k,j)+ch(i,k,jc)
               c1(i,k,jc) = ch(i-1,k,jc)-ch(i-1,k,j)
            enddo
         enddo
      enddo
      go to 121
c
  119 do ik=1,idl1
         c2(ik,1) = ch2(ik,1)
      enddo
c
  121 do j=2,ipph
         jc = ipp2-j
         do k=1,l1
            c1(1,k,j) = ch(1,k,j)+ch(1,k,jc)
            c1(1,k,jc) = ch(1,k,jc)-ch(1,k,j)
         enddo
      enddo
c
      ar1 = 1.d0
      ai1 = 0.d0
      do 127 l=2,ipph
         lc = ipp2-l
         ar1h = dcp*ar1-dsp*ai1
         ai1 = dcp*ai1+dsp*ar1
         ar1 = ar1h
         do ik=1,idl1
            ch2(ik,l) = c2(ik,1)+ar1*c2(ik,2)
            ch2(ik,lc) = ai1*c2(ik,ip)
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
               ch2(ik,l) = ch2(ik,l)+ar2*c2(ik,j)
               ch2(ik,lc) = ch2(ik,lc)+ai2*c2(ik,jc)
            enddo
         enddo
  127 continue
c
      do j=2,ipph
         do ik=1,idl1
            ch2(ik,1) = ch2(ik,1)+c2(ik,j)
         enddo
      enddo
c
      if (ido .lt. l1) go to 132
      do k=1,l1
         do i=1,ido
            cc(i,1,k) = ch(i,k,1)
         enddo
      enddo
      go to 135
c
  132 do i=1,ido
         do k=1,l1
            cc(i,1,k) = ch(i,k,1)
         enddo
      enddo
c
  135 do j=2,ipph
         jc = ipp2-j
         j2 = j+j
         do k=1,l1
            cc(ido,j2-2,k) = ch(1,k,j)
            cc(1,j2-1,k) = ch(1,k,jc)
         enddo
      enddo
c
      if (ido .eq. 1) return
      if (nbd .lt. l1) go to 141
      do j=2,ipph
         jc = ipp2-j
         j2 = j+j
         do k=1,l1
            do i=3,ido,2
               ic = idp2-i
               cc(i-1,j2-1,k) = ch(i-1,k,j)+ch(i-1,k,jc)
               cc(ic-1,j2-2,k) = ch(i-1,k,j)-ch(i-1,k,jc)
               cc(i,j2-1,k) = ch(i,k,j)+ch(i,k,jc)
               cc(ic,j2-2,k) = ch(i,k,jc)-ch(i,k,j)
            enddo
         enddo
      enddo
      return
c
  141 do j=2,ipph
         jc = ipp2-j
         j2 = j+j
         do i=3,ido,2
            ic = idp2-i
            do k=1,l1
               cc(i-1,j2-1,k) = ch(i-1,k,j)+ch(i-1,k,jc)
               cc(ic-1,j2-2,k) = ch(i-1,k,j)-ch(i-1,k,jc)
               cc(i,j2-1,k) = ch(i,k,j)+ch(i,k,jc)
               cc(ic,j2-2,k) = ch(i,k,jc)-ch(i,k,j)
            enddo
         enddo
      enddo
c
      return
      end
