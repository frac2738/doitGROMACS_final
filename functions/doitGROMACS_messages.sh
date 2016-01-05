#!/bin/bash 

helpMessage() {                                                                 
  cat <<EOF                                                                     
                                                                                
                        doitGROMACS.sh -  v 2.0                               
                                                                                
    Copyright (c) 2013-2014, Francesco Carbone, University College London (UCL) 
                                                                                
    This script is designed to automatise the first step of a molecular dynamics
    experiment ( solvation and equilibration) and to run some basic analyses on 
    the trajectory files (.xtc), using GROMACS tools. This script allow also a  
    basic analyses of UNRES (coarse grained ff) trajectories (-u).              
    This script  is written with GROMACS 4.6 in mind and although it should be  
    compatible with any previous versions, it is advise  to check the  validity 
    of each commands before the use with an older version.                      
                                                                                
    At every run the scripts checks the existace of a config file located in the
    working directory (doitGROMACS.config). If a config file is not found, the  
    script uses doitGROMACS_default.config as template for a new configuration  
    file specific for the machine in which the script is run from by adding the 
    correct path to both gromacs and R. If the machine is not specify with the  
    "-b" flag, the config filw will contain standard paths.                     
    NOTE 1: The script recognise four configurations: acrm/emerald/lappy/default.
            If gromacs or R is installed in a non standard location, the user   
            have to manually edit the config file to match its system, otherwise
            doitGROMACS will not find the binaries.                             
                                                                                
    USAGE:  ./doitGROMACS.sh -h                     -->   HELP                  
            ./doitGROMACS.sh -u                     -->   UNRES analysis        
            ./doitGROMACS.sh -b arg1 -n arg2 -...   -->   GROMACS analyses      
                                                                                
                                                                                
    Option   Type     Value       Description                                   
    --------------------------------------------------------------------------- 
    -[no]h   bool     yes         Print help info                               
    -u       string   txt file    Analyses of an unres trajectory               
    -g       bool     bool        Set gromacs 5 syntax                          
                                                                                
                                          [ALWAYS REQUIRED in absence of -u]    
    --------------------------------------------------------------------------- 
    -b       string   acrm        Set the location of binaries                  
                                  acrm      -> Darwin building computer         
                                  emerald   -> Emerald cluster                  
                                  lappy     -> Personal laptop                  
    -n       int      wt          Set the name                                  
                                                                                
                                          [OPTIONAL (function dependant)]       
    --------------------------------------------------------------------------- 
    -t       int      200         Set the simulation length                     
    -s       string   .tpr        .tpr file                                     
    -f       string   .xtc        Trajectory file                               
    -c       string   .pdb        Pdb file to use to start a simulation         
    -e       string   .edr        Energy file                                   
                                                                                
    NOTE 2:  In my simualtions all the output are printed in this format:       
                              NAMErX_TIME                                       
            where  NAME is the name of the mutation ( 306r ), rX is the replica 
            number (r1,r2,...) and TIME is the simulation time.As a consequence 
            this script takes and process output names in this form.            
                                                                                
EOF
}

doitOptions() {                                                                 
   cat <<EOF                                                                    
                                 -----------                                    
    ---------------------------- doitOPTIONS ----------------------------       
                                 -----------                                    
                                                                                
        all           - Starting from scratch **                                
        emin          - Starting from E-minimization **                         
        nvt           - Starting from NVT **                                    
        npt           - Starting from NPT **                                    
        h20           - Remove water from a trajectory file                     
        cond          - Check the simulation conditions (U-T-P-density)         
        rmsdfg        - Calculate RMSD, GYRATION RADIUS and RMSF [backbone & sidechains]
        dssp          - DSSP analysis    
        contact       - Compute contact maps                                       
        cluster       - Cluster analysis                                        
        pca           - PCA analysis                                            
        sas           - SAS analysis                                            
        sas-sites     - SAS analysis on only the binding sites                  
        hb            - Hydrogen bonds analysis [not yet implemented]           
        hb-sites      - Hydrogen bonds analysis on binding sites
        ggplot        - Plot with ggplot (R)                                    
        indexCreator  - Create binding sites index for the mutant       
          
        mean          - calculate the mean of values [2nd column]                    
        mean_multi    - calculate the mean of values [several columns]           
        modvim+       - replace "@" with "#" in a file                          
        catomain      - rebuild a full atoms structure from CA structure        
        split_states  - given an unres trj extract all the frames and convert   
                        into all-atom structures             
        PDBtoXTC      - convert a pdb file (several structures) into a xtc                    
                                                                                
    --------------------------------------------------------------------        
                               -----------                                      
                                                                                
    **  Options that require a parameter file (.mdp) that MUST be placed in the 
        same directory; the functions only accept these files:                  
                                                                                
                          TEMP_min.mdp                                          
                          TEMP_nvt.mdp                                          
                          TEMP_npt.mdp                                          
                          TEMP_md.mdp                                           
                                                                                
        with TEMP = temperature in Kelvin (310, 400, ...)                       
        The doitGROMACS tarball contains some examples of mdp file that can be  
        used after editing.                                                     
                                                                                
EOF
}

gromacs_ver_mex() {                                                             
  cat <<EOF                                                                     
                        -----------------------------                           
    ------------------- GROMACS $gromacs_ver syntax will be used -------------------
                        -----------------------------                           
EOF
}

