#!/usr/bin/env python

import numpy as np
import matplotlib.pyplot as plt
import scipy.io.netcdf as netcdf
plt.ion()


# x,y subsampling (1: no sub sampling, 2: every two points, etc)
sp_xy = 1
# select layer (0: upper layer) 
ilayer = 0

dir0  = '../examples/double_gyre_ocean_only/outdata/'
file0 = 'ocpo.nc'
file1 = 'ocpo.dat'
file2 = 'mean.dat'

fid1 = open(dir0 + file1,'wb')

# get sizes
f = netcdf.netcdf_file(dir0 + file0,'r')
psi = f.variables['p'][0,0,:,:].copy()
xp  = f.variables['xp'][:].copy()
yp  = f.variables['yp'][:].copy()
zi  = f.variables['zi'][:].copy()

time = f.variables['time'][:].copy()
si_t, = time.shape
newt = f.variables['p'][:,0,1,1].copy()
newt = newt[newt<1e5]
si_t2 = newt.size

si_y,si_x = psi.shape

# compute mean
psi_me  = 0.0*np.zeros((si_y,si_x))
n_me = 0
for nt in range(0,si_t2):
  psi = f.variables['p'][nt,ilayer,:,:].copy()
  psi_me = psi_me + psi
  n_me   = n_me + 1
  psi_sav = psi[0:-1:sp_xy,0:-1:sp_xy].reshape(-1)
  np.savetxt(fid1, psi_sav[None],fmt='%10.5e')

# save mean
psi_me = psi_me/n_me
np.savetxt(dir0 + file2, psi_me[0:-1:sp_xy,0:-1:sp_xy])

# close files
f.close()
fid1.close()

# plot
plt.figure()
plt.contour(psi_me)
