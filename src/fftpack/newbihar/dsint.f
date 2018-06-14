      subroutine dsint (n,x,wsave)

*     Subroutine arguments
      integer n
      double precision x(n+1), wsave(*)

*     Local parameters
      double precision SQRT3
      PARAMETER ( SQRT3 = 1.7320508075 6887729 d0 )

*     Local variables
      integer np1,ns2,k,modn,i
      double precision xkc, t1, t2, x1, xh, xim1

*     Cases tested with most likely first
      if ( (n-2).gt.0 ) then
         np1 = n+1
         ns2 = n/2
         x1 = x(1)
         x(1) = 0.0d0
         do k=1,ns2
            xkc = x(np1-k)
            t1 = x1 - xkc
            t2 = wsave(k)*( x1 + xkc )
            x1 = x(k+1)
            x(k+1) = t1+t2
            x(np1+1-k) = t2-t1
         enddo
         modn = mod(n,2)
         if (modn .ne. 0) x(ns2+2) = 4.0d0*x1
c
         call drfftf (np1,x,wsave(ns2+1))
c
         x(1) = 0.5d0*x(1)
         do i=3,n,2
            xim1 = x(i-1)
            x(i-1) = -x(i)
            x(i) = x(i-2)+xim1
         enddo
         if (modn.eq.0) x(n) = -x(n+1)
       else if ( (n-2).eq.0 ) then
         xh = SQRT3*(x(1)+x(2))
         x(2) = SQRT3*(x(1)-x(2))
         x(1) = xh
       else
         x(1) = x(1)+x(1)
      endif

      return
      end
