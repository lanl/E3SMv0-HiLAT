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

cat >! $CASEROOT/Buildconf/pop2conf/sectdyes_tavg_contents << EOF
$s1  DYE01
$s1  DYE02
$s1  DYE03
$s1  DYE04
$s1  DYE05
$s1  DYE06
$s1  STF_DYE01
$s1  STF_DYE02
$s1  STF_DYE03
$s1  STF_DYE04
$s1  STF_DYE05
$s1  STF_DYE06
$s1  FvPER_DYE01
$s1  FvPER_DYE02
$s1  FvPER_DYE03
$s1  FvPER_DYE04
$s1  FvPER_DYE05
$s1  FvPER_DYE06
EOF
