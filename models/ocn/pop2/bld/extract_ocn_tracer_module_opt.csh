#!/bin/csh -f

#
# extract values from key-value pairs in OCN_TRACER_MODULES_OPT
#
# extracted value is written to last line of stdout, with no trailing newline
# usage is to pipe stdout through 'tail -n 1' to get extracted value
# this allows debugging write statements to be added to script without disrupting its usage
#
# known keys, with default value in parentheses, are:
#   IRF_MODE (offline_transport)
#   IRF_NT (depends on IRF_MODE and OCN_GRID)
#   ECOSYS_NT (27)
#   ZOOPLANKTON_CNT (1)
#   AUTOTROPH_CNT (3)
#   GRAZER_PREY_CNT (3)
#

set key = $1

#
# set default value
#

set retval = unknown

# IRF module option defaults

if ($key == IRF_MODE) set retval = offline_transport
if ($key == IRF_NT) then
  set retval = 125 # default for grids with no overflows
  set IRF_MODE = `$0 IRF_MODE | tail -n 1`
  if ($OCN_GRID == gx3v7) set retval = 156
  if ($OCN_GRID == gx1v6) set retval = 178
  if ($OCN_GRID == tx0.3v2) set retval = 15
endif

#
# parse OCN_TRACER_MODULES_OPT for specified value
#

foreach module_opt ( `echo $OCN_TRACER_MODULES_OPT` )
  if ( `echo $module_opt | cut -f 1 -d =` == $key ) then
    set retval = `echo $module_opt | cut -f 2 -d =`
  endif
end

echo -n $retval

exit 0
