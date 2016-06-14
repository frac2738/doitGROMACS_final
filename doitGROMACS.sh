#!/bin/bash 

#------------------------------------------------------------------------------
#
#   File:       doitGROMACS.sh          
#   Version:    V2.1                                                    
#   Update:     13.Jun.16                                                  
#
#   Copyright:  (c) Francesco Carbone, UCL, 2013-2016
#   Author:     Francesco Carbone, UCL                                    
#   Function:   Script to execute a bunch of stuff with gromacs           
#   Address:    Institute of Structural and Molecular Biology
#               Division of Biosciences
#               University College London
#               Gower Street
#               WC1E 6BT
#   EMail:      f.carbone.12@ucl.ac.uk
#
#------------------------------------------------------------------------------
#
#   Decription:
#   This script was design to automatise some parts of md simualtions with 
#   Gromacs (v4.6 and v5). This script is capable to perform energy minimisation
#   equilibration steps (NVT and NPT) and run some standards analyses on the 
#   trajectories. 
#
#------------------------------------------------------------------------------
#
#   Usage:
#   ./doitGROMACS.sh [-u xxx] [-b xxx] [-n xxx] [-t xxx] [-k xxx] [-s xxx] 
#                    [-f xxx] [-c xxx] [-e xxx]
#   -g              - if present gromacs 5 syntax will be used
#   -u              - Analyse an unrest trajectory
#   -b              - Set binary location (gromacs and R)
#   -n              - Set the name
#   -t              - Set the simulation length
#   -s              - Set the tpr file
#   -f              - Set the trajectory file
#   -c              - Set the pdb file
#   -e              - Set the energy file
#
#------------------------------------------------------------------------------

#---------------------------- The program begins here --------------------------

# find from where the script is lunch.
DIR=$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )
FUNCTIONS_BIN="$DIR/functions"
optionRprog="$FUNCTIONS_BIN/doitRGROMACS.R"

# source all the functions ifiles
source $FUNCTIONS_BIN/doitGROMACS_errorsHandling.sh
source $FUNCTIONS_BIN/doitGROMACS_routines.sh
source $FUNCTIONS_BIN/doitGROMACS_messages.sh
source $FUNCTIONS_BIN/doitGROMACS_equilib.sh
source $FUNCTIONS_BIN/doitGROMACS_analyses.sh
source $FUNCTIONS_BIN/doitGROMACS_Rcalls.sh
source $FUNCTIONS_BIN/doitGROMACS_extra.sh

while getopts "hgzb:n:t:s:f:c:e:" opt; do
 case $opt in
    h) helpMessage; exit    ;;
    z) unres=$OPTARG        ;;
    b) cpu=$OPTARG          ;;
    n) name1=$OPTARG        ;;
    t) timens=$OPTARG       ;;
    s) tpr=$OPTARG          ;;
    f) trj=$OPTARG          ;;
    c) pdb1=$OPTARG         ;;
    e) energy=$OPTARG       ;;
    g) gromacs_ver='5'      ;;
    \?) helpMessage;  exit  ;;
 esac
done

checkFlags
# check existance of the parameters CONFIG_FILE and source it or create a new one
export CONFIG_FILE="$(find . -maxdepth 1 -name doitGROMACS.config)"
if [[ -f $CONFIG_FILE ]]; then
  source $CONFIG_FILE 
else
  echo "Configuration file not found, a new file will be created"
  cp $DIR/doitGROMACS_default.config ./doitGROMACS.config
  source ./doitGROMACS.config
fi

# check existance of the binary BINARY_FILE and source it or create a new one
export BINARY_FILE="$(find . -maxdepth 1 -name doitGROMACS_binaries.config)" 
if [[ -f $BINARY_FILE ]]; then
  source $BINARY_FILE
else
  echo "Binary file not found, a new file will be created"
  cp $DIR/doitGROMACS_binaries.config . 
  source ./doitGROMACS_binaries.config
fi
groPATH4="groPATH4_${cpu}"; eval groPATH4=\$$groPATH4
groPATH5="groPATH5_${cpu}"; groPATH5=${!groPATH5}
REXE="REXE_${cpu}"; REXE=${!REXE}
RscriptEXE="RscriptEXE_${cpu}"; RscriptEXE=${!RscriptEXE}

# set the gromacs syntax (version 4 or 5)
setGROMACSbinaries
# list the options 
doitOptions
# check the existance of the selected option 
read -e -p "execute option  " choice
case $choice in
  all|emin|nvt|npt|h20|cond|rmsdfg|cluster|pca|sas|sas-sites|dssp|contact|hb|hb-sites|ggplot|ggplot-bis|indexCreator|modvim+|mean|mean_multi|omega)
    if [ -z ${timens} ]; then
      timens="X"
      nameprod=${name1}_${timens}
    else
      nameprod=${name1}_${timens} 
      fi ;;
  catomain)
    checkFlags_pdb; CAtoMAIN $pdb1  ;;
  split_states)
    checkFlags_pdb; split_states $pdb1  ;;
  PDBtoXTC)
    checkFlags_pdb; PDBtoXTC  ;;
  dssp_perc)
    gromDSSPpercentage  ;;
  dssp_perc)
    gromDSSPpercentage  ;;
  *)
  doitOptions
  error_exit " line $LINENO, An error has occurred. Execution halted! choice '$choice' not recognised."  ;;
esac

case $choice in
  all )
    inputs && energy_minimization && nvt && npt  ;;
  emin)
    energy_minimization && nvt && npt   ;;
  nvt)
    nvt && npt  ;;
  npt)
    npt   ;;
  h20)
    clean_trj   ;;
  cond) 
    sim_conditions ;;
  rmsdfg)
    rmsdf ;;
  cluster)
    gromCLUSTER ;;
  pca)
    gromPCA   ;;
  sas)
    gromSAS  ;;
  sas-sites)
    gromSAS-sites  ;;
  dssp)
    gromDSSP ;;
  contact)
    gromCONTACT ;;
  hb)
    gromHB ;;         
  hb-sites)
    gromHB-sites ;;         
  ggplot)
      GGplot ;;
  ggplot-bis) # hidden function
    sim_conditions && rmsdf && gromDSSP && GGplot ;;
  indexCreator)
      indexCreator ;;
  modvim+)
      modVim_plus  ;;
  mean)
      mean_single  ;;
  mean_multi)
      mean_multi  ;;
  omega)
      gromOMEGA   ;;
esac
#----------------------------- The program ends here ---------------------------
