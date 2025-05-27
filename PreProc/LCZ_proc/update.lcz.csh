#! /bin/csh -f

# Setup Pointers
##############################################################################################################
# All processing relative to NLCD40-LCZ proc directory
setenv NLCDLCZ_DIR /work/MOD3DEV/pleim/WPS/NLCD40_LCZ

# Geogrid settings. Executable. Rob's GEODATA dir, index files for 40 and 51 NLCD40
setenv GEOEXE      ${NLCDLCZ_DIR}/geogrid/geogrid.exe
setenv GEODATA_DIR ${NLCDLCZ_DIR}/GEODATA/nlcd2011_ll_9s
setenv INDEX40     ${GEODATA_DIR}/index.40
setenv INDEX51     ${GEODATA_DIR}/index.51
##############################################################################################################


set DO61    = F
set DO51    = F 
set DO40    = F
set DOBLEND = T
setenv DO d04

### Steps to generate 51 class NLCD40 geogrid file for a domain

##############################################################################################################
# 1. Modeler needs to create the namelists defined below (40 & 61 class) for their domains.
#    The only differences in these files are the LCZ LU spec for the 61 class file and 40 class
#    NLCD40 for the 40 and 51 class namelists. Domain settings need to be identical. 
setenv WPSNL_NLCD40 ${NLCDLCZ_DIR}/namelist.wps.LCZF40
setenv WPSNL_LCZ61  ${NLCDLCZ_DIR}/namelist.wps.LCZF61
##############################################################################################################

##############################################################################################################
# 2. Run geogrid with a 61 class MODIS LCZ specification
if($DO61 == T) then
 cd ${NLCDLCZ_DIR}/RUNDIR
 cp ${WPSNL_LCZ61} namelist.wps
 cp ${NLCDLCZ_DIR}/geogrid/GEOGRID.TBL.ARW_LCZ ${NLCDLCZ_DIR}/geogrid/GEOGRID.TBL.ARW
 ${GEOEXE}
 mv geo_em.d0* ../GEO61/.
endif
##############################################################################################################

##############################################################################################################
# 3. Run geogrid with a 40 class NLCD40 specification using default 40 LU class index file 
##############################################################################################################
if($DO40 == T) then
 cd ${NLCDLCZ_DIR}/RUNDIR
 cp ${WPSNL_NLCD40} namelist.wps
 cp ${INDEX40} ${GEODATA_DIR}/index
 cp ${NLCDLCZ_DIR}/geogrid/GEOGRID.TBL.ARW_DEFAULT ${NLCDLCZ_DIR}/geogrid/GEOGRID.TBL.ARW
 ${GEOEXE}
 mv geo_em.d0* ../GEO40/.
endif


##############################################################################################################
# 4. Run geogrid with a 51 class NLCD40 specification using index file in raw NLCD40 data directory
##############################################################################################################
if($DO51 == T) then
 cd ${NLCDLCZ_DIR}/RUNDIR
 cp ${WPSNL_NLCD40} namelist.wps
 cp ${INDEX51} ${GEODATA_DIR}/index
 cp ${NLCDLCZ_DIR}/geogrid/GEOGRID.TBL.ARW_DEFAULT ${NLCDLCZ_DIR}/geogrid/GEOGRID.TBL.ARW
 ${GEOEXE}
 mv geo_em.d0* ../GEO51/.
endif


##############################################################################################################
# 5. Take LCZ parts (51-61) class MODIS LCZ geogrid output and put in the 51 class NLCD40 geogrid 41-51
##############################################################################################################
if($DOBLEND == T) then

 setenv GEO_MODIS61 ${NLCDLCZ_DIR}/GEO61/geo_em.${DO}.nc
 setenv GEO_NLCD40  ${NLCDLCZ_DIR}/GEO40/geo_em.${DO}.nc
 setenv GEO_NLCD51  ${NLCDLCZ_DIR}/GEO51/geo_em.${DO}.nc

 R BATCH --no-save < ${NLCDLCZ_DIR}/nlcd_lcz.R

endif
##############################################################################################################
##############################################################################################################
#    
#  After geogrid step, proceed with WPS like any other NLCD40-based input file preperation.
#  Keep all number of landuse cat in namelist.input as 40, not 51. This will keep PX LSM running
#  normal, but allow access to extra LCZ for any urban considerations.
#
##############################################################################################################








