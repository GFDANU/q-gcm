The Quasi-Geostrophic Coupled Model (Q-GCM)
# Compilation
in src/
```
ln -s ../confing/make.macro .
ln -s ../confing/make.config .
ln -s ../confing/parameter_data.F .
make
```


# Run the double gyre ocean only example
create the forcing file: 
in python/
```
python create_forcing.py
cp avges.nc ../examples/double_gyre_ocean_only/
```

in src
```
rm make.config parameter_data.F
ln -s ../examples/double_gyre_ocean_only/make.config .
ln -s ../examples/double_gyre_ocean_only/parameter_data.F .
make`
```

in examples/double_gyre_ocean_only/
```
mkdir outdata
./q-gcm
```