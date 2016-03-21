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

cat >! $CASEROOT/Buildconf/pop2conf/surfdyes_tavg_contents << EOF
$s1  DYE_PREC
$s1  DYE_MELT
$s1  DYE_ROFF
$s1  DYE_IOFF
$s1  J_DYE_MELT
$s1  J_DYE_ROFF
$s1  J_DYE_PREC
$s1  STF_DYE_MELT
$s1  STF_DYE_ROFF
$s1  STF_DYE_PREC
$s1  FvPER_DYE_MELT
$s1  FvPER_DYE_ROFF
$s1  FvPER_DYE_PREC
$s1  DYE_MELT_zint
$s1  DYE_ROFF_zint
$s1  DYE_PREC_zint
$s1  tend_zint_DYE_MELT
$s1  tend_zint_DYE_ROFF
$s1  tend_zint_DYE_PREC
EOF
