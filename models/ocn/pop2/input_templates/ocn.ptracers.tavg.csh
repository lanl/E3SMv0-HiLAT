#!/bin/csh -f

#------------------------------------------------------------------------------------
# For now, set streams manually. You must only set as many streams as are declared
#  in the tavg_nml section. For example, if there are three streams:
#  @ s1 = $my_stream
#  @ s2 = $s1 + 1
#  @ s3 = $s2 + 1
#------------------------------------------------------------------------------------

@ my_stream = $1
if ($my_stream < 1) then
   echo invalid my_stream number $my_stream
   exit 5
endif

@ s1 = 1   # use base-model stream 1

cat >! $CASEROOT/Buildconf/pop2conf/ptracers_tavg_contents << EOF
$s1  pTEMP
$s1  pSALT
$s1  STF_pTEMP
$s1  STF_pSALT
$s1  VDC_pT
$s1  VDC_pS
$s1  KPP_SRC_pTEMP
$s1  KPP_SRC_pSALT
$s1  pHMXL
$s1  pXMXL
$s1  pTMXL
$s1  pHBLT
$s1  pXBLT
$s1  pTBLT
EOF
