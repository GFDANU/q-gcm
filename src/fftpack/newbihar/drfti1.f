      subroutine drfti1 (n,wa,ifac)

*     Subroutine arguments
      integer n, ifac(*)
      double precision wa(*)

*     Local parameters
      double precision TPI
      PARAMETER ( TPI = 6.2831853071 7958648 d0 )

*     Local variables
      integer ntryh(4), nl, nf, j, ntry, nq, nr, i, ib, is,
     &        nfm1, l1, k1, ip, ld, l2, ido, ipm, ii
      double precision arg, argh, argld, fi
c
      data ntryh(1), ntryh(2), ntryh(3), ntryh(4) /4, 2, 3, 5/
c
      nl = n
      nf = 0
      j = 0
c
  101 j = j+1
      if (j.le.4) ntry = ntryh(j)
      if (j.gt.4) ntry = ntry + 2
  104 nq = nl/ntry
      nr = nl-ntry*nq
      if (nr.ne.0) go to 101
c
      nf = nf+1
      ifac(nf+2) = ntry
      nl = nq
      if (ntry .ne. 2) go to 107
      if (nf .eq. 1) go to 107
      do i=2,nf
         ib = nf-i+2
         ifac(ib+2) = ifac(ib+1)
      enddo
      ifac(3) = 2
  107 if (nl .ne. 1) go to 104
      ifac(1) = n
      ifac(2) = nf
c
      argh = TPI/dble(n)
      is = 0
      nfm1 = nf-1
      l1 = 1
      if (nfm1 .eq. 0) return
      do 110 k1=1,nfm1
         ip = ifac(k1+2)
         ld = 0
         l2 = l1*ip
         ido = n/l2
         ipm = ip-1
         do j=1,ipm
            ld = ld+l1
            i = is
            argld = dble(ld)*argh
            fi = 0.d0
            do ii=3,ido,2
               i = i+2
               fi = fi+1.d0
               arg = fi*argld
               wa(i-1) = dcos(arg)
               wa(i) = dsin(arg)
            enddo
            is = is+ido
         enddo
c
         l1 = l2
  110 continue
c
      return
      end
