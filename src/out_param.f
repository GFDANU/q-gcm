c***********************************************************************
c     Q-GCM Version 1.5.0 : last modified 18/09/2013
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
*     Write parameters out in matlab format file called
*     input_parameters.m in the output directory

      write(10,'(a)') '%%Matlab script to read in parameters'

*     Preprocessor options set in Makefile
*     ------------------------------------
      write(10,122) 'oceanonly = ',oceanonly,';  %% Ocean only run?'
      write(10,122) 'atmosonly = ',atmosonly,';  %% Atmos only run?'
      write(10,122) 'getcovar = ',getcovar,';   %% Get covar in run?'
      write(10,122) 'cyclicoc = ',cyclicoc,';   %% Cyclic ocean?'
      write(10,122) 'hflxsb = ',hflxsb,';     %% S boundary heat flux?'
      write(10,122) 'hflxnb = ',hflxnb,';     %% N boundary heat flux?'
      write(10,122) 'tauudiff = ',tauudiff,';   %% Use oc. vel. in tau?'

*     Input parameters from parameters_data.F
*     ---------------------------------------
      write(10,122) 'nxto= ',nxto,';        %% Ocean x gridcells'
      write(10,122) 'nyto= ',nyto,';        %% Ocean y gridcells'
      write(10,122) 'nlo= ',nlo,';         %% Ocean QG layers'
      write(10,122) 'nxta= ',nxta,';        %% Atmos. x gridcells'
      write(10,122) 'nyta= ',nyta,';        %% Atmos. y gridcells'
      write(10,122) 'nla= ',nla,';         %% Atmos. QG layers'
      write(10,122) 'nxaooc= ',nxaooc,
     &              ';      %% Atmos. x gridcells over ocean'
      write(10,122) 'nyaooc= ',nyaooc,
     &              ';      %% Atmos. y gridcells over ocean'
      write(10,122) 'ndxr= ',ndxr,
     &              ';        %% Atmos./Ocean gridlength ratio'
      write(10,122) 'nx1= ',nx1,';         %% Starting index for',
     &              ' ocean in atmospheric grid'
      write(10,122) 'ny1= ',ny1,';         %% Starting index for',
     &              ' ocean in atmospheric grid'

      write(10,121) 'fnot= ',fnot,';     %% Coriolis parameter'
      write(10,121) 'beta= ',beta,';     %% Beta'

      if ( getcovar.eq.1 ) then
        write(10,122) 'nscvoc= ',nscvoc,';      %% Ocean x subsampling'
        write(10,122) 'nvcvoc= ',nvcvoc,';      %% Ocean subsamp vector'
        write(10,122) 'nmcvoc= ',nmcvoc,';      %% Ocean covar components'
        write(10,122) 'nscvat= ',nscvat,';      %% Atmos. x subsampling'
        write(10,122) 'nvcvat= ',nvcvat,';      %% Atmos. subsamp vector'
        write(10,122) 'nmcvat= ',nmcvat,';      %% Atmos. covar components'
      endif

*     Run parameters from input.params
*     --------------------------------
      write(10,121) 'tini= ',tini,';     %% Start time in years'
      write(10,121) 'trun= ',trun,';     %% Run length in years'
      write(10,121) 'tend= ',tend,';     %% Final time in years'
      write(10,121) 'dto= ',dto,';      %% Ocean timestep in seconds'
      write(10,121) 'dta= ',dta,';      %% Atmos. timestep in seconds'
      write(10,121) 'dxo= ',dxo,';      %% Ocean grid spacing in km'
      write(10,121) 'dxa= ',dxa,';      %% Atmos. grid spacing in km'

      write(10,121) 'delek= ',delek,
     &              ';    %% Ocean bottom Ekman thickness'
      write(10,121) 'cdat= ',cdat,';     %% Air-Sea',
     &              ' momentum drag coefficient'
      write(10,121) 'rhoat= ',rhoat,';    %% Atmos. density'
      write(10,121) 'rhooc= ',rhooc,';    %% Ocean  density'
      write(10,121) 'cpat= ',cpat,';     %% Atmos. spec. heat capacity'
      write(10,121) 'cpoc= ',cpoc,';     %% Ocean spec. heat capacity'
      write(10,121) 'bccoat= ',bccoat,
     &              ';   %% Mixed BC coefficient for atmos (nondim.)'
      write(10,121) 'bccooc= ',bccooc,
     &              ';   %% Mixed BC coefficient for ocean (nondim.)'
      write(10,121) 'xcexp= ',xcexp,
     &              ';    %% coupling coefficient x'
      write(10,121) 'ycexp= ',ycexp,
     &              ';    %% coupling coefficient y'

*     Data dumping parameters
*     -----------------------
      write(10,121) 'valday= ',valday,
     &              ';   %% Solution check interval (days)'
      write(10,121) 'odiday= ',odiday,
     &              ';   %% Ocean  data dump interval (days)'
      write(10,121) 'adiday= ',adiday,
     &              ';   %% Atmos. data dump interval (days)'
      write(10,121) 'dgnday= ',dgnday,
     &              ';   %% Diagnostics dump interval (days)'
      write(10,121) 'resday= ',resday,
     &              ';   %% Restart dump interval (days)'
      write(10,122) 'noutoc= ',noutoc,';     %% Output interval: ocean'
      write(10,122) 'noutat= ',noutat,';     %% Output interval: atmos.'
      write(10,122) 'nsko= ',nsko,';        %% Subsampling interval',
     &              ' for ocean output'
      write(10,122) 'nska= ',nska,';        %% Subsampling interval',
     &              ' for atmos. output'
      write(10,121) 'dtavat= ',dtavat,
     &              ';   %% Atmos. averaging int. (days)'
      write(10,121) 'dtavoc= ',dtavoc,
     &              ';   %% Ocean  averaging int. (days)'

*     Mixed layer parameters
*     ----------------------
      write(10,121) 'hmoc= ',hmoc,';     %% Fixed ocean  ml depth'
      write(10,121) 'hmat= ',hmat,';     %% Fixed atmos. ml depth'
      write(10,121) 'st2d= ',st2d,';     %% sst del-sqd diffusivity'
      write(10,121) 'st4d= ',st4d,';     %% sst del-4th diffusivity'
      write(10,121) 'ahmd= ',ahmd,';     %% hmixa lateral diffusivity'
      write(10,121) 'at2d= ',at2d,';     %% ast del-sqd diffusivity'
      write(10,121) 'at4d= ',at4d,';     %% ast del-4th diffusivity'
      write(10,121) 'tsbdy= ',tsbdy,';    %% o.m.l. S. bdy. temp (rel)'
      write(10,121) 'xlamda= ',xlamda,';   %% Sensible/latent transfer'
      write(10,121) 'hmadmp= ',hmadmp,';   %% At. mixed layer h damping'

*     Radiation parameters
*     --------------------
      write(10,121) 'fsbar= ',fsbar,
     &              ';    %% Mean radiation forcing'
      write(10,121) 'fspamp= ',fspamp,';   %% Radiation perturbation'
      write(10,121) 'zm= ',zm,';       %% Optical depth in a.m.l.'
      write(10,121) 'zopt= ',zopt(1),
     &              ';     %% Optical depth in layer 1'
      do k=2,nla
        write(10,121) 'zopt= [zopt ',zopt(k),'];   %% Layers 2,n'
      enddo
      write(10,121) 'gamma= ',gamma,';    %% Adiabatic lapse rate'

*     Oceanic QG layer parameters
*     ---------------------------
*     Reduced gravities, temperatures, friction coeffs and thicknesses
      write(10,121) 'gpoc= ',gpoc(1),
     &              ';     %% Reduced gravity for ocean 1/2 interface'
      do k=2,nlo-1
        write(10,121) 'gpoc= [gpoc ',gpoc(k),'];   %% Interfaces 2,n-1'
      enddo
      write(10,121) 'ah2oc= ',ah2oc(1),';    %% Del-sqd coefft ocean'
      do k=2,nlo
        write(10,121) 'ah2oc= [ah2oc ',ah2oc(k),']; %% Layers 2,n'
      enddo
      write(10,121) 'ah4oc= ',ah4oc(1),';    %% Del-4th coefft ocean'
      do k=2,nlo
        write(10,121) 'ah4oc= [ah4oc ',ah4oc(k),']; %% Layers 2,n'
      enddo
      write(10,121) 'tabsoc= ',tabsoc(1),
     &              ';   %% Abs. temperature for ocean layer 1'
      do k=2,nlo
        write(10,121) 'tabsoc= [tabsoc ',tabsoc(k),'];    %% Layers 2,n'
      enddo
      write(10,121) 'tocc= ',toc(1),
     &              ';     %% Rel. temperature for ocean layer 1'
      do k=2,nlo
        write(10,121) 'tocc= [tocc ',toc(k),'];   %% Layers 2,n'
      enddo
      write(10,121) 'hoc= ',hoc(1),
     &              ';      %% Thickness of ocean layer 1'
      do k=2,nlo
        write(10,121) 'hoc= [hoc ',hoc(k),'];     %% Layers 2,n'
      enddo

*     Atmospheric QG layer parameters
*     -------------------------------
*     Reduced gravities, temperatures, friction coeffs and thicknesses

      write(10,121) 'gpat= ',gpat(1),
     &              ';     %% Reduced gravity for atmos 1/2 interface'
      do k=2,nla-1
        write(10,121) 'gpat= [gpat ',gpat(k),'];   %% Interfaces 2,n-1'
      enddo
      write(10,121) 'ah4at= ',ah4at(1),';    %% Del-4th coefft atmos'
      do k=2,nla
        write(10,121) 'ah4at= [ah4at ',ah4at(k),'];   %% Layers 2,n'
      enddo
      write(10,121) 'tabsat= ',tabsat(1),
     &              ';   %% Abs. temperature for atmos layer 1'
      do k=2,nla
        write(10,121) 'tabsat= [tabsat ',tabsat(k),'];    %% Layers 2,n'
      enddo
      write(10,121) 'tat= ',tat(1),
     &              ';      %% Rel. temperature for atmos layer 1'
      do k=2,nla
        write(10,121) 'tat= [tat ',tat(k),'];     %% Layers 2,n'
      enddo
      write(10,121) 'hat= ',hat(1),
     &              ';      %% Thickness of atmos layer 1'
      do k=2,nla
        write(10,121) 'hat= [hat ',hat(k),'];     %% Layers 2,n'
      enddo

      write(10,'(3a)') 'name= ''',name(1:lename),
     &              ''';            %% Initial condition file'

      write(10,125) 'outfloc= [', outfloc(1), outfloc(2), outfloc(3),
     &              outfloc(4), outfloc(5), outfloc(6), outfloc(7),
     &              '];    %% output flag vector for ocean'
      write(10,125) 'outflat= [', outflat(1), outflat(2), outflat(3),
     &              outflat(4), outflat(5), outflat(6), outflat(7),
     &              '];    %% output flag vector for atmos.'

*     Derived parameters
*     ==================
      write(10,'(a)') '%%Derived parameters'
*     Mixed layers
*     ------------
      write(10,121) 'tmbara= ',tmbara,';   %% Actually T_{mlao}'
      write(10,121) 'tmbaro= ',tmbaro,';   %% Actually T_{mloo}'
*     Ocean
*     -----
      write(10,121) 'cphsoc= ',cphsoc(2),
     &              ';   %% Baroclinic wavespeed for ocean mode 1'
      do k=3,nlo
        write(10,121) 'cphsoc= [cphsoc ',cphsoc(k),
     &                '];   %% Higher modes'
      enddo
      write(10,121) 'rdefoc= ',rdefoc(2),
     &              ';   %% Deformation radius for ocean mode 1'
      do k=3,nlo
        write(10,121) 'rdefoc= [rdefoc ',rdefoc(k),
     &                '];   %% Higher modes'
      enddo
      write(10,121) 'tsbdy= ',tsbdy,
     &              ';    %% Rel. o.m.l. S. bndry temp. (K)'
      write(10,121) 'tnbdy= ',tnbdy,
     &              ';    %% Rel. o.m.l. N. bndry temp. (K)'
*     Atmosphere
*     ----------
      write(10,121) 'cphsat= ',cphsat(2),
     &              ';   %% Baroclinic wavespeed for atmos mode 1'
      do k=3,nla
        write(10,121) 'cphsat= [cphsat ',cphsat(k),
     &                '];   %% Higher modes'
      enddo
      write(10,121) 'rdefat= ',rdefat(2),
     &              ';   %% Deformation radius for atmos mode 1'
      do k=3,nla
        write(10,121) 'rdefat= [rdefat ',rdefat(k),
     &                '];   %% Higher modes'
      enddo
      write(10,121) 'aface= ',aface(1),
     &              ';    %% eta    coefficient aface(1)'
      do k=2,nla-1
        write(10,121) 'aface= [aface ',aface(k),
     &                '];   %% Other interfaces'
      enddo
      write(10,121) 'bface= ',bface,
     &              ';    %% etam   coefficient bface'
      write(10,121) 'cface= ',cface,
     &              ';    %% topog. coefficient cface'
      write(10,121) 'dface= ',dface,
     &              ';    %% aTm    coefficient dface'

  121 format(A,1P,E13.5,A,A)
  122 format(A,I10,A,A)
  125 format(A,7I2,A)
*
