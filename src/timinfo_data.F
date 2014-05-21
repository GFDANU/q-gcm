c***********************************************************************
c     Q-GCM Version 1.5.0 : last modified 24/09/2013
c***********************************************************************
c
c     Copyright 2013 Jeff Blundell, Andy Hogg and Bill Dewar.
c     This file is part of Q-GCM.
c
c     Q-GCM is free software: you can redistribute it and/or modify
c     it under the terms of the GNU General Public License as
c     published by the Free Software Foundation, either version 3
c     of the License, or (at your option) any later version.
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
c
      MODULE timinfo

*     Contains parameters and variables controlling run length,
*     timestepping, diagnostics, validity testing, monitoring, etc.

      IMPLICIT NONE

      PUBLIC
      SAVE

      double precision, PARAMETER :: secday = 86400.0d0
      double precision, PARAMETER :: daysyr = 365.0d0
      double precision, PARAMETER :: secsyr = secday*daysyr

*     Double precision variables defining time intervals for
*     run length, diagnostics, validity testing, monitoring, etc.
*     (apart from the actual timestep intervals, dto & dta
*      which are in occonst.mod & atconst.mod respectively)
      double precision :: trun,tini,tend,tday,tyrs,
     &                    valday,odiday,adiday,dgnday,prtday,
     &                    resday,dtavat,dtavoc,dtcovat,dtcovoc

*     Integer variables for controlling timestepping, averaging,
*     diagnostic calculations, monitoring, validity checking, etc.
      integer :: nstr,nvalid,noutoc,noutat,nocmon,nocprt,noutre,
     &           ntavat,ntavoc,nmidat,nmidoc,nsteps0,nsteps,
     &           ntsrun,ntdone,noutstepoc,noutstepat,numoutsteps

#ifdef get_covar
*     Extra variables for covariance computation
      integer :: ntcovat,ntcovoc
#endif /* get_covar */

*     Double precision variables:
*     ---------------------------
*     trun = length of the current run (years)                  read from input.params
*     tini = time at the start of the of the current run (years);
*          read from a restart file, or set to zero for a 'zero' or 'rbal' start
*     tend = time at the end of the of the current run (years)  set in main program
*     tday = current time from the start of the total run,      set in main program
*            not just the current run segment (days)
*     tyrs = current time from the start of the total run,      set in main program
*            not just the current run segment (years)
*     valday  = validity checking interval (days)               read from input.params
*     odiday  = ocean  data dump interval  (days)               read from input.params
*     adiday  = atmos. data dump interval  (days)               read from input.params
*     dgnday  = diagnostics dump interval  (days)               read from input.params
*     prtday  = print to std output interval (days)             read from input.params
*     resday  = restart data dump interval (days)               read from input.params
*     dtavat  = Atmos. averaging interval  (days) (zero => off) read from input.params
*     dtavoc  = Ocean  averaging interval  (days) (zero => off) read from input.params
*     dtcovat = Atmos. covariance interval (days) (zero => off) read from input.params
*     dtcovoc = Ocean  covariance interval (days) (zero => off) read from input.params

*     Integer variables:
*     ------------------
*     (Note that the model counts in units of atmospheric timesteps, since it is
*      assumed that these will always be shorter than oceanic timesteps, by the factor nstr)
*     nstr   = timestep ratio dto/dta  (derives dto from dta)   read from input.params
*     nvalid = validity checking interval         (timesteps)   derived from valday
*     noutoc = ocean  data dump interval          (timesteps)   derived from odiday
*     noutat = atmos. data dump interval          (timesteps)   derived from adiday
*     nocmon = diagnostic data dump interval      (timesteps)   derived from dgnday
*     nocprt = print to std output interval       (timesteps)   derived from prtday
*     noutre = restart data dump interval         (timesteps)   derived from resday
*     ntavat = ocean  averaging interval          (timesteps)   derived from dtavat
*     ntavoc = atmos. averaging interval          (timesteps)   derived from dtavoc
*     nmidat = half-width of atmos. avg. interval (timesteps)   derived from ntavat
*     nmidoc = half-width of ocean  avg. interval (timesteps)   derived from ntavoc
*     nsteps0 = no. of timesteps already completed at restart   derived from tini
*     (inferred assuming the runs up to restart used same timestep
*      length as the current run, which need not have been the case)
*     nsteps = no. of timesteps to be completed at end of run   derived from tend
*     ntsrun = no. of timesteps to be performed within the run  derived from trun
*     ntdone = no. of timesteps completed within current run    set in main program
*     noutstepoc = no. of ocean  data dumps within current run  set in main program
*     noutstepat = no. of atmos. data dumps within current run  set in main program
*     numoutsteps = no. of diagnostic dumps within current run  set in main program
*
      END MODULE timinfo
c
c***********************************************************************
