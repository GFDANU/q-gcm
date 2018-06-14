      subroutine drffti (n,wsave)

*     Subroutine arguments
      integer n
      double precision wsave(*)
c
      if (n .eq. 1) return
c
      call drfti1 (n,wsave(n+1),wsave(2*n+1))
c
      return
      end
