c***********************************************************************
c     Q-GCM Version 1.5.0 : last modified 02/07/2013
c
c     11/10/2017: added parameters related to high-res. atmospheric output
c                 over ocean (Sergey Kravtsov)
c
c***********************************************************************

      MODULE highresout

      IMPLICIT NONE

      PUBLIC
      SAVE

      integer :: nxbd,nybd,ndxr2
      PARAMETER ( nxbd = 0, nybd = 0, ndxr2 = 1)

*     nxbd, nybd are the numbers of land points to be added, in each direction,
*     to the atmospheric find-grid output over ocean
*     ndxr2 is the ratio of the atmospheric output gridlength over ocean
*     to the fine-grid atmospheric gridlength. Can be 2 (for 10-km atmospheric
*     resolution), otherwise 1.

*
*     Additional derived parameters for output dimensions:
*
      integer :: nxtout,nytout,nxpout,nypout
      PARAMETER ( nxtout = nxaooc/ndxr2 + 2*nxbd,
     &            nytout = nyaooc/ndxr2 + 2*nybd,
     &            nxpout = nxtout + 1, nypout = nytout + 1 )

*     nxtout, nytout are the dimensions of the high-resolution atmospheric
*     output over ocean; nxpout, nypout - the same on the p-grid
*
      END MODULE highresout 
c
c***********************************************************************


