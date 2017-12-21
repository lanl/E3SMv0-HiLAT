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

cat >! $CASEROOT/Buildconf/pop2conf/pptracers_tavg_contents << EOF
$s1  ppTEMP
$s1  ppSALT
$s1  STF_ppTEMP
$s1  STF_ppSALT
$s1  KPP_SRC_ppTEMP
$s1  KPP_SRC_ppSALT
EOF
