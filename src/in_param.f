c***********************************************************************
c     Q-GCM Version 1.5.0 : last modified 24/09/2013
c***********************************************************************
c
c     Copyright 2013 Jeff Blundell, Andy Hogg and Bill Dewar.
c     This file is part of Q-GCM.
c
c     Q-GCM is free software: you can redistribute it and/or modify
c     it under the terms of the GNU General Public License as published
c     by the Free Software Foundation, either version 3 of the License,
c     or (at your option) any later version.
c
c     Q-GCM is distributed in the hope that it will be useful,
c     but WITHOUT ANY WARRANTY; without even the implied warranty
c     of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
c     See the GNU General Public License for more details.
c
c     You should have received a copy of the GNU General Public License
c     along with Q-GCM.  If not, see <http://www.gnu.org/licenses/>.
c
c***********************************************************************
*
*     Read name of output directory
*     -----------------------------
      open (odunit, file='./outdata.dat')
      read (odunit,'(A)') outdir
      close(odunit)
*
*     Read input parameters
*     ---------------------
      open (ipunit, file='./input.params')
      call ipbget (inpbuf, ipunit)
      read (inpbuf,*) trun
      call ipbget (inpbuf, ipunit)
      read (inpbuf,*) dta
      call ipbget (inpbuf, ipunit)
      read (inpbuf,*) nstr
      call ipbget (inpbuf, ipunit)
      read (inpbuf,*) dxo
      call ipbget (inpbuf, ipunit)
      read (inpbuf,*) delek
      call ipbget (inpbuf, ipunit)
      read (inpbuf,*) cdat
      call ipbget (inpbuf, ipunit)
      read (inpbuf,*) rhoat
      call ipbget (inpbuf, ipunit)
      read (inpbuf,*) rhooc
      call ipbget (inpbuf, ipunit)
      read (inpbuf,*) cpat
      call ipbget (inpbuf, ipunit)
      read (inpbuf,*) cpoc
      call ipbget (inpbuf, ipunit)
      read (inpbuf,*) bccoat
      call ipbget (inpbuf, ipunit)
      read (inpbuf,*) bccooc
      call ipbget (inpbuf, ipunit)
      read (inpbuf,*) xcexp
      call ipbget (inpbuf, ipunit)
      read (inpbuf,*) ycexp
      call ipbget (inpbuf, ipunit)
      read (inpbuf,*) valday
      call ipbget (inpbuf, ipunit)
      read (inpbuf,*) odiday
      call ipbget (inpbuf, ipunit)
      read (inpbuf,*) adiday
      call ipbget (inpbuf, ipunit)
      read (inpbuf,*) dgnday
      call ipbget (inpbuf, ipunit)
      read (inpbuf,*) prtday
      call ipbget (inpbuf, ipunit)
      read (inpbuf,*) resday
      call ipbget (inpbuf, ipunit)
      read (inpbuf,*) nsko
      call ipbget (inpbuf, ipunit)
      read (inpbuf,*) nska
      call ipbget (inpbuf, ipunit)
      read (inpbuf,*) dtavat
      call ipbget (inpbuf, ipunit)
      read (inpbuf,*) dtavoc
      call ipbget (inpbuf, ipunit)
      read (inpbuf,*) dtcovat
      call ipbget (inpbuf, ipunit)
      read (inpbuf,*) dtcovoc
      call ipbget (inpbuf, ipunit)
      read (inpbuf,*) xlamda
      call ipbget (inpbuf, ipunit)
      read (inpbuf,*) hmoc
      call ipbget (inpbuf, ipunit)
      read (inpbuf,*) st2d
      call ipbget (inpbuf, ipunit)
      read (inpbuf,*) st4d
      call ipbget (inpbuf, ipunit)
      read (inpbuf,*) hmat
      call ipbget (inpbuf, ipunit)
      read (inpbuf,*) hmamin
      call ipbget (inpbuf, ipunit)
      read (inpbuf,*) ahmd
      call ipbget (inpbuf, ipunit)
      read (inpbuf,*) at2d
      call ipbget (inpbuf, ipunit)
      read (inpbuf,*) at4d
      call ipbget (inpbuf, ipunit)
      read (inpbuf,*) hmadmp
      call ipbget (inpbuf, ipunit)
      read (inpbuf,*) fsbar
      call ipbget (inpbuf, ipunit)
      read (inpbuf,*) fspamp
      call ipbget (inpbuf, ipunit)
      read (inpbuf,*) zm
      call ipbget (inpbuf, ipunit)
      read (inpbuf,*) zopt
      call ipbget (inpbuf, ipunit)
      read (inpbuf,*) gamma
      call ipbget (inpbuf, ipunit)
      read (inpbuf,*) ah2oc
      call ipbget (inpbuf, ipunit)
      read (inpbuf,*) ah4oc
      call ipbget (inpbuf, ipunit)
      read (inpbuf,*) tabsoc
      call ipbget (inpbuf, ipunit)
      read (inpbuf,*) hoc
      call ipbget (inpbuf, ipunit)
      read (inpbuf,*) gpoc
      call ipbget (inpbuf, ipunit)
      read (inpbuf,*) ah4at
      call ipbget (inpbuf, ipunit)
      read (inpbuf,*) tabsat
      call ipbget (inpbuf, ipunit)
      read (inpbuf,*) hat
      call ipbget (inpbuf, ipunit)
      read (inpbuf,*) gpat
      call ipbget (inpbuf, ipunit)
      read (inpbuf,'(A)') name
      call ipbget (inpbuf, ipunit)
      read (inpbuf,'(A)') topocname
      call ipbget (inpbuf, ipunit)
      read (inpbuf,'(A)') topatname
      call ipbget (inpbuf, ipunit)
      read (inpbuf,*) outfloc
      call ipbget (inpbuf, ipunit)
      read (inpbuf,*) outflat
      close(ipunit)
*
