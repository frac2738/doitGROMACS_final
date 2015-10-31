#!/bin/bash

#-------
GGplot() {
  # function that calls the R script doitRGROMACS.R to plot in ggplot: 
  # rmsd - gyration radius - rmsf - simulation conditions - structure analysis
  if [[ -f $optionRprog ]]; then
    $RscriptEXE $optionRprog -o=$nameprod -d=$nameprod"_rmsd.xvg"              \
      -g=$nameprod"_rgyration.xvg" -ss=$nameprod"_ss_count.xvg"                \
      -x=$nameprod"_density.xvg" -t=$nameprod"_temperature.xvg"                \
      -u=$nameprod"_potential.xvg" -p=$nameprod"_pressure.xvg"                 \
      -fb=$nameprod"_rmsf_bb.xvg" -fsc=$nameprod"_rmsf_sc.xvg"                 \
      -hm="pes_profile.txt"
    $RscriptEXE $DIR/scriptino.R $DIR/schematics_bindingSites                  \
      $DIR/schematics_motileRegions $nameprod"_rmsf_bb.xvg"                    \
      $nameprod"_rmsf_schem"
  else 
    error_exit " the function "GGplot" requires a R script located in $DIR."
  fi 
} &> >(tee doitgromacs_ggplot.log) >&2


