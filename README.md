The Quasi-Geostrophic Coupled Model (Q-GCM)
# Compile and run the double gyre ocean only case

## Create the forcing file: 

in python/
```
python create_forcing.py
\cp -f avges.nc ../examples/double_gyre_ocean_only/
```

## Link the config files and compile

in src
```
ln -s ../config/make.macro.x86_64.gfortran431.netcdf3 make.macro
rm make.config parameter_data.F
ln -s ../examples/double_gyre_ocean_only/make.config .
ln -s ../examples/double_gyre_ocean_only/parameters_data.F .
make
```

in examples/double_gyre_ocean_only/
```
mkdir outdata
cp ../../src/q-gcm .
./q-gcm
```