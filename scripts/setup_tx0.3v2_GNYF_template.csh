#!/bin/bash

# Make sure to set:
# the correct advective scheme ('lw_lim' in user_nl_pop2);
# diffusive settings ('del2', ah&hmix_del2t_nml = 0.0e7 in user_nl_pop2);
# the local salinity fix for run-off (in forcing_coupled.F90); 
# and the correct value of ksno (in ice_constants_colpkg.F90)
#
# The decomposition is appropriate for 1920 OCN pes, and 480 ICE pes
# 
export CaseName=t32_GNYF_Griz
echo "casename = $CaseName"
export CaseDir=/turquoise/usr/projects/climate/jiaxu/tx0.3v2_IAF/
echo "casedir = $CaseDir/$CaseName"

./create_newcase -case $CaseDir/$CaseName -res T62_t32 -compset G -mach grizzly -proj climatehilat
cd $CaseDir/$CaseName

./xmlchange -file env_build.xml -id POP_AUTO_DECOMP -val false
./xmlchange -file env_build.xml -id POP_BLCKX -val 25
./xmlchange -file env_build.xml -id POP_BLCKY -val 20
./xmlchange -file env_build.xml -id POP_NX_BLOCKS -val 48
./xmlchange -file env_build.xml -id POP_NY_BLOCKS -val 40
./xmlchange -file env_build.xml -id POP_MXBLCKS -val 4
./xmlchange -file env_build.xml -id POP_DECOMPTYPE -val cartesian
./xmlchange -file env_build.xml -id CICE_MODE -val prognostic
./xmlchange -file env_build.xml -id CICE_AUTO_DECOMP -val false
./xmlchange -file env_build.xml -id CICE_BLCKX -val 5
./xmlchange -file env_build.xml -id CICE_BLCKY -val 400
./xmlchange -file env_build.xml -id CICE_MXBLCKS -val 1
./xmlchange -file env_build.xml -id CICE_DECOMPTYPE -val cartesian
./xmlchange -file env_build.xml -id CICE_DECOMPSETTING -val slenderX2

./xmlchange -file env_run.xml -id ATM_NCPL -val 24
./xmlchange -file env_run.xml -id LND_NCPL -val 24
./xmlchange -file env_run.xml -id ICE_NCPL -val 24
./xmlchange -file env_run.xml -id OCN_NCPL -val 24
./xmlchange -file env_run.xml -id GLC_NCPL -val 24
./xmlchange -file env_run.xml -id ROF_NCPL -val 24
./xmlchange -file env_run.xml -id WAV_NCPL -val 24
./xmlchange -file env_run.xml -id CPL_SEQ_OPTION -val RASM_OPTION1

./xmlchange -file env_run.xml -id PIO_TYPENAME -val pnetcdf
./xmlchange -file env_run.xml -id OCN_PIO_STRIDE -val -1
./xmlchange -file env_run.xml -id OCN_PIO_ROOT -val 0
./xmlchange -file env_run.xml -id OCN_PIO_NUMTASKS -val -1
./xmlchange -file env_run.xml -id OCN_PIO_TYPENAME -val pnetcdf
./xmlchange -file env_run.xml -id LND_PIO_STRIDE -val -1
./xmlchange -file env_run.xml -id LND_PIO_ROOT -val -1
./xmlchange -file env_run.xml -id LND_PIO_NUMTASKS -val -1
./xmlchange -file env_run.xml -id LND_PIO_TYPENAME -val pnetcdf
./xmlchange -file env_run.xml -id ROF_PIO_STRIDE -val -1
./xmlchange -file env_run.xml -id ROF_PIO_ROOT -val -1
./xmlchange -file env_run.xml -id ROF_PIO_NUMTASKS -val -1
./xmlchange -file env_run.xml -id ROF_PIO_TYPENAME -val pnetcdf
./xmlchange -file env_run.xml -id ICE_PIO_STRIDE -val -1
./xmlchange -file env_run.xml -id ICE_PIO_ROOT -val -1
./xmlchange -file env_run.xml -id ICE_PIO_NUMTASKS -val -1
./xmlchange -file env_run.xml -id ICE_PIO_TYPENAME -val pnetcdf
./xmlchange -file env_run.xml -id ATM_PIO_STRIDE -val -1
./xmlchange -file env_run.xml -id ATM_PIO_ROOT -val -1
./xmlchange -file env_run.xml -id ATM_PIO_NUMTASKS -val -1
./xmlchange -file env_run.xml -id ATM_PIO_TYPENAME -val pnetcdf
./xmlchange -file env_run.xml -id CPL_PIO_STRIDE -val -1
./xmlchange -file env_run.xml -id CPL_PIO_ROOT -val -1
./xmlchange -file env_run.xml -id CPL_PIO_NUMTASKS -val -1
./xmlchange -file env_run.xml -id CPL_PIO_TYPENAME -val pnetcdf
./xmlchange -file env_run.xml -id GLC_PIO_STRIDE -val -1
./xmlchange -file env_run.xml -id GLC_PIO_ROOT -val -1
./xmlchange -file env_run.xml -id GLC_PIO_NUMTASKS -val -1
./xmlchange -file env_run.xml -id GLC_PIO_TYPENAME -val pnetcdf
./xmlchange -file env_run.xml -id WAV_PIO_STRIDE -val -1
./xmlchange -file env_run.xml -id WAV_PIO_ROOT -val -1
./xmlchange -file env_run.xml -id WAV_PIO_NUMTASKS -val -1
./xmlchange -file env_run.xml -id WAV_PIO_TYPENAME -val pnetcdf
