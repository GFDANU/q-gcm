      subroutine drfftf (n,r,wsave)

*     Subroutine arguments
      integer n
      double precision r(n), wsave(*)
c
      if (n .eq. 1) return
c
      call drftf1 (n,r,wsave,wsave(n+1),wsave(2*n+1))
c
      return
      end
