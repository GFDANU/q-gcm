      subroutine dsinti (n,wsave)

*     Subroutine arguments
      integer n
      double precision wsave(*)

*     Local parameters
      double precision PI
      PARAMETER ( PI = 3.141592653 58979324 d0 )

*     Local variables
      integer np1,ns2,k
      double precision dt, fk
c
      if (n .le. 1) return
      np1 = n+1
      ns2 = n/2
      dt = PI/dble(np1)
      fk = 0.d0
      do k=1,ns2
         fk = fk+1.d0
         wsave(k) = 2.d0*dsin(fk*dt)
      enddo
c
      call drffti (np1,wsave(ns2+1))
c
      return
      end
