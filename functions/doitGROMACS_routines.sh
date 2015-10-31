#!/bin/bash 

# set gromacs binaries depending on the version called (-g flag)
setGROMACSbinaries() {                                                          
  if [ -n "${gromacs_ver}" ]; then                                              
    groPATH=$groPATH5; gromacs_ver_mex                                          
    pdb2gmx='gmx pdb2gmx'; editconf='gmx editconf';                             
    genbox='gmx solvate'; grompp='gmx grompp'; genion='gmx genion'              
    g_energy='gmx energy'; trjconv='gmx trjconv'; tpbconv='gmx convert-tpr'     
    g_rms='gmx rms'; g_rmsf='gmx rmsf'; g_gyrate='gmx gyrate'                   
    xpm2ps='gmx xpm2ps'; g_cluster='gmx cluster'; g_covar='gmx covar'           
    g_anaeig='gmx anaeig'; g_analyze='gmx analyze'; g_sham='gmx sham'           
    g_sas='gmx sasa'; g_hbond='gmx hbond'; do_dssp='gmx do_dssp'                
  else                                                                          
    groPATH=$groPATH4                                                           
    gromacs_ver='4.6'; gromacs_ver_mex                                          
    pdb2gmx='pdb2gmx'; editconf='editconf';                                     
    genbox='genbox'; grompp='grompp'; genion='genion'                           
    g_energy='g_energy'; trjconv='trjconv'; tpbconv='tpbconv'                   
    g_rms='g_rms'; g_rmsf='g_rmsf'; g_gyrate='g_gyrate'                         
    xpm2ps='xpm2ps'; g_cluster='g_cluster'; g_covar='g_covar'                   
    g_anaeig='g_anaeig'; g_analyze='g_analyze'; g_sham='g_sham'                 
    g_sas='g_sas'; g_hbond='g_hbond'; do_dssp='do_dssp'                         
  fi                                                                            
}                                                                               

# description : this function takes a ".stat" file from an UNRES simulation and 
#               passes it to an R script for analyses.  
unresAnalyses() {                                                               
   $RSexePATH $rScripts/doitUNRES.R $unres                                      
}                                                                               
                                                                                
# description: Replace all the "@" with "#" using sed. Usefull to avoid conflicts
#              between R and grace.                                             
modVim() {                                                                      
  sed -i 's/@/#/g' "$1" # -i modify and save                                    
}                                                                               
                                                                                
modVim_plus() {                                                                 
  read -e -p "che file vuoi? " file1                                            
  modVim $file1                                                                 
}                                                                                     
                
checkFlags() {                                                                  
  # check the definition of -b and -n                                           
  if [[ -z "${cpu+x}" || -z "${name1+x}" ]]; then                               
    helpMessage                                                                 
    error_exit "Execution halted! the script was called without any flags."     
  fi                                                                            
}                                                                               
                                                                                
checkFlags_t() {                                                                
  if [[ -z "${tpr+x}" || -z "${trj+x}" ]]; then                                 
    helpMessage; error_exit " execution halted: tpr and trj files are required (flags -s -f)"
  fi                                                                            
}                                                                               
                                                                                
checkFlags_pdb() {                                                              
  if [ -z "${pdb1}" ]; then                                                     
    helpMessage; error_exit " execution halted: a pdb file is required (-c)"    
  fi                                                                            
}                                                                               
                                                                                
checkFlags_energy() {                                                           
  if [ -z "${energy}" ]; then                                                   
    helpMessage; error_exit " execution halted: an energy file is required (-e)"
  fi                                                                            
}  

