#!/bin/csh -f

#cat >! $CASEBUILD/pop2conf/tracegas.tavg.nml << EOF
#tavg_freq_opt             = 'nday'   'nyear'
#tavg_freq                 =  1       1
#tavg_file_freq_opt        = 'nmonth' 'nyear'
#tavg_file_freq            =  1       1
#tavg_start_opt            = 'nstep'  'nstep'
#tavg_start                =  0       0
#tavg_fmt_in               = 'nc'     'nc'
#tavg_fmt_out              = 'nc'     'nc'
#ltavg_has_offset_date     = .false.  .false.
#tavg_offset_years         =  1       1
#tavg_offset_months        =  1       1
#tavg_offset_days          =  2       2
#ltavg_one_time_header     = .false.  .false.
#tavg_stream_filestrings   = 'tracegas.nday1' 'tracegas.nyear1'
#EOF

#------------------------------------------------------------------------------------
# For now, set streams manually. You must only set as many streams as are declared
#  in the tavg_nml section. For example, if there are three streams:
#  @ s1 = $my_stream
#  @ s2 = $s1 + 1
#  @ s3 = $s2 + 1
#  (SW) for two group of Phaeocystis
#------------------------------------------------------------------------------------

@ my_stream = $1
if ($my_stream < 1) then
   echo invalid my_stream number  ($my_stream)
   exit 5
endif

@ s1 = 1             # use the base-model stream 1
@ s2 = $my_stream    # use an tracegastem-defined stream
@ s3 = $s2 + 1       # use an tracegastem-defined stream

cat >! $CASEROOT/Buildconf/pop2conf/tracegas_tavg_contents << EOF
1  DMS
1  DMSP
1  DMS_S_TOTAL
1  DMS_S_DMSP
1  DMS_R_TOTAL
1  DMS_R_B
1  DMS_R_PHOT
1  DMS_R_BKGND
1  DMSP_S_TOTAL
1  DMSP_S_PHAEO
1  DMSP_S_NONPHAEO
1  DMSP_S_ZOO
1  DMSP_R_TOTAL
1  DMSP_R_B
1  DMSP_R_BKGND
1  DMS_WS
1  DMS_SURF
1  DMS_SAT
1  DMS_SCHMIDT
1  STF_DMS
1  Cyano_frac
1  Cocco_frac
1  Eukar_frac
1  diatS
1  diatN
1  phytoN
1  coccoS
1  cyanoS
1  eukarS
1  diazS
1  phaeoS
1  phaeonS
EOF
