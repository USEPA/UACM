#This file still needs adjustment so it creates new variables instead of updating old ones
from netCDF4 import Dataset
import numpy as np

#Open the files
bldg_input = Dataset('PLUTO_data_FAI.nc', 'r', format='NETCDF4')
bldgHeight = bldg_input.variables['BuildingHeight'][:,:]
bldgFAI = bldg_input.variables['FrontalAreaIndex'][:,:]
bldgPAFraction = bldg_input.variables['PlanAreaDensity'][:,:] 
bldgCount = bldg_input.variables['COUNT'][:,:] #How many buildings are in each grid cell. May be useful, is not necessary.
bldg_input.close()

print("new data imported")

wrf_input = Dataset('wrfinput_d01', 'a', format='NETCDF4')
BUILD_HEIGHT = wrf_input.createVariable('BUILD_HEIGHT','f4',('Time','south_north','west_east'))
LAMF = wrf_input.createVariable('LAMF','f4',('Time','south_north','west_east'))
LAMP = wrf_input.createVariable('BUILD_AREA_FRACTION','f4',('Time','south_north','west_east'))

BUILD_HEIGHT.FieldType = int(104)
LAMF.FieldType = int(104)
LAMP.FieldType = int(104)
BUILD_HEIGHT.MemoryOrder = 'XY'
LAMF.MemoryOrder = 'XY'
LAMP.MemoryOrder = 'XY'
BUILD_HEIGHT.description = 'Average building height'
LAMF.description = 'Frontal area index'
LAMP.description = 'BUILDING PLAN AREA DENSITY'
BUILD_HEIGHT.units = 'm'
LAMF.units = 'dimensionless'
LAMP.units = 'dimensionless'
BUILD_HEIGHT.stagger = ""
LAMF.stagger = ""
LAMP.stagger = ""
BUILD_HEIGHT.coordinates = "XLONG XLAT XTIME"
LAMF.coordinates = "XLONG XLAT XTIME"
LAMP.coordinates = "XLONG XLAT XTIME"


BUILD_HEIGHT[0,:,:] = bldgHeight
LAMF[0,:,:] = bldgFAI
LAMP[0,:,:] = bldgPAFraction
wrf_input.close()

#Old code block for updating WRF variables
#wrfBldgHeight = wrf_input.variables['BUILD_HEIGHT'][0,:,:]
#wrfFAI = wrf_input.variables['LAMF'][0,:,:] #is this what the variable should be called?
#wrfPAFraction = wrf_input.variables['BUILD_AREA_FRACTION'][0,:,:]
#print("wrf data imported")
#for i in range(0,162):
#    for j in range(0,162):
#        if (bldgHeight[i,j] > 0):
#            wrfBldgHeight[i,j] = bldgHeight[i,j]
#            wrfFAI[i,j] = bldgFAI[i,j]
#            wrfPAFraction[i,j] = bldgPAFraction[i,j]
#print("Loop complete")
#write the new variables in to the existing file, be sure to do so for all times (which for wrfinput is 1)
#times = np.size(landmask, axis=0)
#print("Applying changes to " + str(times) + " timesteps")
#for i in range(0,times): 
#    wrf_input['BUILD_HEIGHT'][i,:,:] = wrfBldgHeight[:,:]
#    wrf_input['LAMF'][i,:,:] = wrfSARatio[:,:]
#    wrf_input['BUILD_AREA_FRACTION'][i,:,:] = wrfPAFraction[:,:]
#wrf_input.close()


print("Script complete")
