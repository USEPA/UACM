require(ncdf4)

 geof_modis61    <- Sys.getenv("GEO_MODIS61")
 geof_nlcd40     <- Sys.getenv("GEO_NLCD40")
 geof_nlcf51     <- Sys.getenv("GEO_NLCD51")

writeLines(paste("Reading input MODIS 61 class landuse file:",geof_modis61))
f1  <-nc_open(geof_modis61)
 luf61 <-ncvar_get(f1, varid="LANDUSEF")
nc_close(f1)
lcz <-luf61[,,51:61]


writeLines(paste("Reading input NLCD40 class landuse file:",geof_nlcd40))
f1 <-nc_open(geof_nlcd40)
 luf40 <-ncvar_get(f1, varid="LANDUSEF")
 lm    <-ncvar_get(f1, varid="LANDMASK")
 lui   <-ncvar_get(f1, varid="LU_INDEX")
nc_close(f1)

writeLines(paste("Updating WRF Geogrid landuse, impervious sfc and canopy fraction:",geof_nlcf51))
f1  <-nc_open(geof_nlcf51, write=T)
luf51 <-ncvar_get(f1, varid="LANDUSEF")
luf51[,,1:40]<-luf40
luf51[,,41:51]<-lcz

ncvar_put(f1, varid="LU_INDEX", lui)
ncvar_put(f1, varid="LANDMASK", lm)
ncvar_put(f1, varid="LANDUSEF", luf51)
nc_close(f1)


#setenv GEO_MODIS61 /work/MOD3DEV/grc/WPS4.0/domains/listos12.4.1.33/GEO_LCZ/61/geo_em.d01.nc
#setenv GEO_NLCD40 /work/MOD3DEV/grc/WPS4.0/domains/listos12.4.1.33/GEO_LCZ/40/geo_em.d01.nc
#setenv GEO_NLCD51 /work/MOD3DEV/grc/WPS4.0/domains/listos12.4.1.33/GEO_LCZ/51/geo_em.d01.nc


