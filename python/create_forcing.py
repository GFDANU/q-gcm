import matplotlib.pyplot as plt
import numpy as np
from scipy.io import netcdf

dirOUT = './'

# physical constants
dx = 10e3

# hogg 10km res
si_x = 384+1;
si_y = 480+1;

taux = np.zeros((si_y,si_x))
tauy = np.zeros((si_y,si_x))
fnet = np.zeros((si_y-1,si_x-1))

tau0 = 2e-5

# Store 
fileout = dirOUT + 'avges.nc'
f = netcdf.netcdf_file(fileout,'w')

f.createDimension('ypo',si_y)
f.createDimension('xpo',si_x)
f.createDimension('yto',si_y-1)
f.createDimension('xto',si_x-1)

ypo = f.createVariable('ypo', 'd', ('ypo',))
xpo = f.createVariable('xpo', 'd', ('xpo',))
yto = f.createVariable('yto', 'd', ('yto',))
xto = f.createVariable('xto', 'd', ('xto',))

tauxo  = f.createVariable('tauxo' , 'd', ('ypo','xpo',))
tauyo  = f.createVariable('tauyo' , 'd', ('ypo','xpo',))
fnetoc = f.createVariable('fnetoc', 'd', ('yto','xto',))

ypo[:] = np.arange(si_y)
xpo[:] = np.arange(si_x)
yto[:] = np.arange(si_y-1)
xto[:] = np.arange(si_x-1)

for ny in range(0,si_y):
  taux[ny,:] = tau0*(-np.cos((ny+0.5)/si_y*2*np.pi))


tauxo [:,:] = taux
tauyo [:,:] = np.zeros((si_y,si_x))
fnetoc[:,:] = np.zeros((si_y-1,si_x-1))

f.close()
