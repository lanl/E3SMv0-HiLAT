#!/bin/bash

# Make sure to set:
# the correct advective scheme ('lw_lim' in user_nl_pop2);
# diffusive settings ('del2', ah&hmix_del2t_nml = 0.0e7 in user_nl_pop2);
# the local salinity fix for run-off (in forcing_coupled.F90); 
# and the correct value of ksno (in ice_constants_colpkg.F90)
#
# Add  init_ts_option = 'file' to user_nl_pop2
#
# The decomposition is appropriate for 1920 OCN pes, and 480 ICE pes
# 
export CaseName=JRA55_t32_Test2_grizzly
echo "casename = $CaseName"
export CaseDir=/turquoise/usr/projects/climate/wilbert/HiLAT-JRA55/
echo "casedir = $CaseDir/$CaseName"

./create_newcase -case $CaseDir/$CaseName -res TL319_t32 -compset GJRA55 -mach grizzly -proj climatehilat
cd $CaseDir/$CaseName

./xmlchange -file env_build.xml -id POP_AUTO_DECOMP -val false
./xmlchange -file env_build.xml -id POP_BLCKX -val 20
./xmlchange -file env_build.xml -id POP_BLCKY -val 20
./xmlchange -file env_build.xml -id POP_NX_BLOCKS -val 60
./xmlchange -file env_build.xml -id POP_NY_BLOCKS -val 40
./xmlchange -file env_build.xml -id POP_MXBLCKS -val 4
./xmlchange -file env_build.xml -id POP_DECOMPTYPE -val cartesian
./xmlchange -file env_build.xml -id CICE_MODE -val prognostic
./xmlchange -file env_build.xml -id CICE_AUTO_DECOMP -val false
./xmlchange -file env_build.xml -id CICE_BLCKX -val 10
./xmlchange -file env_build.xml -id CICE_BLCKY -val 400
./xmlchange -file env_build.xml -id CICE_MXBLCKS -val 1
./xmlchange -file env_build.xml -id CICE_DECOMPTYPE -val cartesian
./xmlchange -file env_build.xml -id CICE_DECOMPSETTING -val slenderX2

./xmlchange -file env_run.xml -id ATM_NCPL -val 32
./xmlchange -file env_run.xml -id LND_NCPL -val 32
./xmlchange -file env_run.xml -id ICE_NCPL -val 32
./xmlchange -file env_run.xml -id OCN_NCPL -val 32
./xmlchange -file env_run.xml -id GLC_NCPL -val 32
./xmlchange -file env_run.xml -id ROF_NCPL -val 32
./xmlchange -file env_run.xml -id WAV_NCPL -val 32
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

./xmlchange -file env_mach_pes.xml -id NTASKS_ATM -val 240
./xmlchange -file env_mach_pes.xml -id ROOTPE_ATM -val 2400
./xmlchange -file env_mach_pes.xml -id NTASKS_LND -val 16
./xmlchange -file env_mach_pes.xml -id ROOTPE_LND -val 2640
./xmlchange -file env_mach_pes.xml -id NTASKS_ICE -val 240
./xmlchange -file env_mach_pes.xml -id ROOTPE_ICE -val 2400
./xmlchange -file env_mach_pes.xml -id NTASKS_OCN -val 2400
./xmlchange -file env_mach_pes.xml -id ROOTPE_OCN -val 0
./xmlchange -file env_mach_pes.xml -id NTASKS_CPL -val 240
./xmlchange -file env_mach_pes.xml -id ROOTPE_CPL -val 2400
./xmlchange -file env_mach_pes.xml -id NTASKS_GLC -val 16
./xmlchange -file env_mach_pes.xml -id ROOTPE_GLC -val 2640
./xmlchange -file env_mach_pes.xml -id NTASKS_ROF -val 240
./xmlchange -file env_mach_pes.xml -id ROOTPE_ROF -val 2400
./xmlchange -file env_mach_pes.xml -id NTASKS_WAV -val 16
./xmlchange -file env_mach_pes.xml -id ROOTPE_WAV -val 2640
