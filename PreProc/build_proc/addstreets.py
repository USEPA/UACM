#This file still needs adjustment so it creates new variables instead of updating old ones
from netCDF4 import Dataset
import numpy as np

#Open the files
bldg_input = Dataset('street_2dir3.nc', 'r', format='NETCDF4')
Dir1 = bldg_input.variables['Dir1'][:,:]
B1 = bldg_input.variables['B1'][:,:]
Dir2 = bldg_input.variables['Dir2'][:,:]
B2 = bldg_input.variables['B2'][:,:]
bldg_input.close()

print("new data imported")

wrf_input = Dataset('wrfinput_d04x', 'a', format='NETCDF4')

if 'STDIR1' not in wrf_input.variables:
    stdir1 = wrf_input.createVariable('STDIR1','f4',('Time','south_north','west_east'))
else:
    stdir1 = wrf_input.variables['STDIR1']
if 'BLKW1' not in wrf_input.variables:
    blkw1 = wrf_input.createVariable('BLKW1','f4',('Time','south_north','west_east'))
else:
    blkw1 = wrf_input.variables['BLKW1']
if 'STDIR2' not in wrf_input.variables:
    stdir2 = wrf_input.createVariable('STDIR2','f4',('Time','south_north','west_east'))
else:
    stdir2 = wrf_input.variables['STDIR2']
if 'BLKW2' not in wrf_input.variables:
    blkw2 = wrf_input.createVariable('BLKW2','f4',('Time','south_north','west_east'))
else:
    blkw2 = wrf_input.variables['BLKW2']


stdir1.FieldType = int(104)
blkw1.FieldType = int(104)
stdir2.FieldType = int(104)
blkw2.FieldType = int(104)
stdir1.MemoryOrder = 'XY'
blkw1.MemoryOrder = 'XY'
stdir2.MemoryOrder = 'XY'
blkw2.MemoryOrder = 'XY'
stdir1.description = 'Most dominant street direction'
blkw1.description = 'Average block width dir stdir1'
stdir2.description = 'Second most dominant street direction'
blkw2.description = 'Average block width dir stdir2'
stdir1.units = 'deg'
blkw1.units = 'm'
stdir2.units = 'deg'
blkw2.units = 'm'
stdir1.stagger = ""
blkw1.stagger = ""
stdir2.stagger = ""
blkw2.stagger = ""
stdir1.coordinates = "XLONG XLAT XTIME"
blkw1.coordinates = "XLONG XLAT XTIME"
stdir2.coordinates = "XLONG XLAT XTIME"
blkw2.coordinates = "XLONG XLAT XTIME"


stdir1[0,:,:] = Dir1
blkw1[0,:,:] = B1
stdir2[0,:,:] = Dir2
blkw2[0,:,:] = B2
wrf_input.close()


print("Script complete")
