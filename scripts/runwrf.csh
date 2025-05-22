#!/bin/csh 
#SBATCH -J test_WRF2023
#SBATCH -p compute
#SBATCH -t 168:00:00
#SBATCH -n 240          
#SBATCH --constraint=cascadelake
#SBATCH -A mod3dev
# -----------------------------------------------------------------------
source /home/jpleim/.cshrc
# NEW Files

set WRF_DIR     = /work/MOD3DEV/pleim/UACM

cd /work/MOD3DEV/pleim/UACM47bench

      if ( -f ETAMPNEW_DATA ) rm -f ETAMPNEW_DATA
      if ( -f GENPARM.TBL   ) rm -f GENPARM.TBL
      if ( -f landFilenames ) rm -f landFilenames
      if ( -f LANDUSE.TBL   ) rm -f LANDUSE.TBL
      if ( -f RRTM_DATA     ) rm -f RRTM_DATA
      if ( -f SOILPARM.TBL  ) rm -f SOILPARM.TBL
      if ( -f tr49t67       ) rm -f tr49t67
      if ( -f tr49t85       ) rm -f tr49t85
      if ( -f tr67t85       ) rm -f tr67t85
      if ( -f VEGPARM.TBL   ) rm -f VEGPARM.TBL

#      ln -s $EXECDIR/wrf.exe               wrf.exe

      ln -s $WRF_DIR/test/em_real/ETAMPNEW_DATA ETAMPNEW_DATA
      ln -s $WRF_DIR/test/em_real/GENPARM.TBL   GENPARM.TBL
#     ln -s $WRF_DIR/test/em_real/landFilenames landFilenames
      ln -s $WRF_DIR/test/em_real/LANDUSE.TBL   LANDUSE.TBL
      ln -s $WRF_DIR/test/em_real/RRTM_DATA     RRTM_DATA
      ln -s $WRF_DIR/test/em_real/RRTMG_SW_DATA RRTMG_SW_DATA
      ln -s $WRF_DIR/test/em_real/RRTMG_LW_DATA RRTMG_LW_DATA
      ln -s $WRF_DIR/test/em_real/SOILPARM.TBL  SOILPARM.TBL
      ln -s $WRF_DIR/test/em_real/tr49t67       tr49t67
      ln -s $WRF_DIR/test/em_real/tr49t85       tr49t85
      ln -s $WRF_DIR/test/em_real/tr67t85       tr67t85
      ln -s $WRF_DIR/test/em_real/VEGPARM.TBL   VEGPARM.TBL
      ln -s $WRF_DIR/test/em_real/ozone_plev.formatted  ozone_plev.formatted
      ln -s $WRF_DIR/test/em_real/ozone_lat.formatted   ozone_lat.formatted
      ln -s $WRF_DIR/test/em_real/ozone.formatted       ozone.formatted
      ln -s $WRF_DIR/test/em_real/CAMtr_volume_mixing_ratio CAMtr_volume_mixing_ratio

# ln -sf /work/MOD3DEV/pleim/WPS/met_em.d0* .
ln -sf wrfinputs_May1-5/wrfbdy_d01 wrfbdy_d01
ln -sf wrfinputs_May1-5/wrffdda_d01 wrffdda_d01
ln -sf wrfinputs_May1-5/wrffdda_d02 wrffdda_d02
ln -sf wrfinputs_May1-5/wrffdda_d03 wrffdda_d03
ln -sf wrfinputs_May1-5/wrffdda_d04 wrffdda_d04
ln -sf wrfinputs_May1-5/wrfinput_d01 wrfinput_d01
ln -sf wrfinputs_May1-5/wrfinput_d02 wrfinput_d02
ln -sf wrfinputs_May1-5/wrfinput_d03 wrfinput_d03
ln -sf wrfinputs_May1-5/wrfinput_d04 wrfinput_d04
ln -sf wrfinputs_May1-5/wrflowinp_d01 wrflowinp_d01
ln -sf wrfinputs_May1-5/wrflowinp_d02 wrflowinp_d02
ln -sf wrfinputs_May1-5/wrflowinp_d03 wrflowinp_d03
ln -sf wrfinputs_May1-5/wrflowinp_d04 wrflowinp_d04
ln -sf wrfinputs_May1-5/wrfsfdda_d01 wrfsfdda_d01
ln -sf wrfinputs_May1-5/wrfsfdda_d02 wrfsfdda_d02
ln -sf wrfinputs_May1-5/wrfsfdda_d03 wrfsfdda_d03


time mpirun wrf.exe
#time mpirun -r ssh -np 128 ./real.exe

exit(0)




