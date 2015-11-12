#!/bin/bash
#---------------
# description : Given an energy file, it extracts energy, temerature, pressure
#               and density profiles.
# requirements: energy.edr
sim_conditions() {
  # if -e is not set prompt and exit
  checkFlags_energy
  # export potential energy, temperature, pressure and density profiles 
  echo Potential | $groPATH/$g_energy -f $energy -o $nameprod"_potential.xvg"
  modVim $nameprod"_potential.xvg"
  echo Temperature | $groPATH/$g_energy -f $energy -o $nameprod"_temperature.xvg"
  modVim $nameprod"_temperature.xvg"
  echo Pressure | $groPATH/$g_energy -f $energy -o $nameprod"_pressure.xvg"
  modVim $nameprod"_pressure.xvg"
  echo Density | $groPATH/$g_energy -f $energy -o $nameprod"_density.xvg"
  modVim $nameprod"_density.xvg"
  # to plot using ggplot check option the function "GGplot"
} &> >(tee doitgromacs_simcond.log) >&2


#----------
# description : It removes the water molecules from a trajectory (.xtc) and 
#               removes the pbc effects (pbc = periodic buonday conditions).
# requirements: .tpr + .trj 
clean_trj() {
   # if -s and -f are not set prompt and exit
   checkFlags_t
   # removing water molecules from the trajectory file
   echo Protein | $groPATH/$trjconv -s $tpr -f $trj                             \
      -o $nameprod"_only.xtc" || checkExitCode
   # removing water molecules from the tpr file 
   echo Protein | $groPATH/$tpbconv -s $tpr -o $nameprod"_only.tpr" || checkExitCode
   # creating a protein only .gro file 
   echo Protein | $groPATH/$trjconv -s $nameprod"_only.tpr"                     \
      -f $nameprod"_only.xtc" -o $nameprod"_only.gro" -dump $optiondump || checkExitCode
   # account for the periodicity (nojump)
   echo Protein | $groPATH/$trjconv -s $nameprod"_only.tpr"                     \
      -f $nameprod"_only.xtc" -o $nameprod"_nojump.xtc" -pbc nojump || checkExitCode
   # account for the periodicity (fitting) 
   (echo "Backbone"; echo "Protein") | $groPATH/$trjconv -s $nameprod"_only.tpr"\
      -f $nameprod"_nojump.xtc" -o $nameprod"_fit.xtc" -fit rot+trans || checkExitCode
   # remove intermediate files and rename them in a less complicated way
   rm $nameprod"_nojump.xtc" $nameprod".xtc" $nameprod"_only.xtc" || checkExitCode
   mv $nameprod"_fit.xtc" $nameprod".xtc" || checkExitCode
   mv $nameprod"_only.tpr" $nameprod".tpr" || checkExitCode
   mv $nameprod"_only.gro" $nameprod".gro" || checkExitCode
} &> >(tee doitgromacs_cleantrj.log) >&2

#------
# description : It calculates rmsd (backbone), gyration radius and rmsf (both (echo "r 213 243 234 234"; echo "name 10 prova"; echo "q") | 
#               for the backbone and the sidechains) 
# requirements: .tpr + .xtc 
rmsdf() {
   checkFlags_t
   # calculating the RMSD 
   (echo "$optionRMSD"; echo "$optionRMSD") | $groPATH/$g_rms -s $tpr -f $trj   \
      -o $nameprod"_rmsd.xvg" -tu ns
   modVim $nameprod"_rmsd.xvg"
   # calculating the radius of Gyration 
   echo $optionGYRATION | $groPATH/$g_gyrate -s $tpr -f $trj                    \
      -o $nameprod"_rgyration.xvg"
   modVim $nameprod"_rgyration.xvg"
   echo $optionRMSFb | $groPATH/$g_rmsf -s $tpr -f $trj                         \
      -o $nameprod"_rmsf_bb.xvg" -oq $nameprod"_rmsf_bb.pdb" -res
   modVim $nameprod"_rmsf_bb.xvg" 
      # res=averages for each residues
   echo $optionRMSFsc | $groPATH/$g_rmsf -s $tpr -f $trj                        \
      -o $nameprod"_rmsf_sc.xvg" -oq $nameprod"_rmsf_sc.pdb" -res 
   modVim $nameprod"_rmsf_sc.xvg" 
} &> >(tee doitgromacs_rmsf.log) >&2

#----------------
# description : read the name of the function
# requirements: .tpr + .xtc
gromCLUSTER() {
   # if -s and -f are not set prompt and exit
   checkFlags_t
   if [ ! -d ./clusters_$name1 ] ; then
      mkdir clusters_$name1
   fi
   cd clusters_$name1
   # check if the rmsd matrix exists
   while [ ! -f ./rmsd-matrix.xpm ]
   do
      # create the rmsd-matrix (using the backbone for the calculation)
      (echo "$optionCLUSTER"; echo "$optionCLUSTER") | $groPATH/$g_rms -s ../$tpr -f ../$trj   \
         -f2 ../$trj -m rmsd-matrix.xpm -dist rmsd-distribution.xvg            \
         -skip $optionSKIPcluster -b $optionSTARTime
      # improve rmsd matrix plot
      $groPATH/$xpm2ps -f rmsd-matrix.xpm -o rmsd-matrix.eps
   done
   #check the distribution file and decide the cutoff
   xmgrace rmsd-distribution.xvg
   read -e -p "Which cutoff do you want to use? " cutoff
   # cluster analysis on the rmsd matrix (Using the backbone for the calculation)
   (echo "$optionCLUSTER"; echo "$optionCLUSTER") | $groPATH/$g_cluster -s ../$tpr -f ../$trj  \
      -dm rmsd-matrix.xpm -o clusters.xpm -sz clusters-size.xvg                \
      -clid clusters-ovt.xvg -cl clusters.pdb -cutoff $cutoff                  \
      -method $optionCLUSTERMETHOD -skip $optionSKIPcluster -b $optionSTARTime
   # to visualize in pymol use
   # split_states clusters
   # delete clusters
   # dss
   # show cartoon
   cd ..
} &> >(tee doitgromacs_cluster.log) >&2

#----------------------
# description : The name says everything
# requirements: see "clusterAnalysis"
repeatgromCLUSTER() {
   mv clusters_$name1 clusters"$i"_$name1
   clusterAnalysis
   read -e -p "Do you want to rerun the analysis with a different method? [yes/no] " ramen
}

#----
# description : The name says everything
# requirements: .tpr + .xtc
gromPCA() {
  # if -s and -f are not set prompt and exit
  checkFlags_t
  if [ ! -d ./PCA_$name1 ] ; then
     mkdir PCA_$name1
  fi
  cd PCA_$name1
  # check if the covariance matrix exists
  while [ ! -f ./covariance.xpm ]
  do
     # calculating the covariance matrix (C-alpha), eigenvalues and eigenvectors
     (echo "$optionPCA"; echo "$optionPCA") | $groPATH/$g_covar -s ../$tpr -f ../$trj   \
        -o eigenvalues.xvg -v eigenvectors.trr -xpma covariance.xpm -tu ps    \
        -dt $optionDTpca -b $optionSTARTime
     # improve covariance matrix ploting
     $groPATH/$xpm2ps -f covariance.xpm -o covariance.eps
  done
  xmgrace eigenvalues.xvg
  # determinare su quali autovalori lavorare
  read -e -p "how many eigenvalues do you want to use for the analysis? " range
  (echo "$optionPCA"; echo "$optionPCA") | $groPATH/$g_anaeig -v eigenvectors.trr       \
     -s ../$tpr -f ../$trj -first 1 -last $range -proj "projection-1"$range".xvg" \
     -tu ps -b $optionSTARTime
  for ((i = 1; i <= range; i=i+1))
  do
     (echo "$optionPCA"; echo "$optionPCA") | $groPATH/$g_anaeig -v eigenvectors.trr \
        -s ../$tpr -f ../$trj -first $i -last $i -nframes 100                 \
        -extr "ev"$i".pdb" -filt "ev"$i".xtc" -proj "ev"$i".xvg" -tu ps       \
        -b $optionSTARTime
     $groPATH/$g_analyze -f "ev"$i".xvg" -cc "ev"$i"-cc.xvg" -ac "ev"$i"-ac.xvg"
  done
  # calculate 2d projections and FES
  (echo "$optionPCA"; echo "$optionPCA") | $groPATH/$g_anaeig -v eigenvectors.trr   \
     -s ../$tpr -f ../$trj -first 1 -last 2 -2d 2d-12.xvg
  $groPATH/$g_sham -f 2d-12.xvg -ls gibbs-12.xpm -notime
  $groPATH/$xpm2ps -f gibbs-12.xpm -o gibbs-12.eps -rainbow red
  perl $DIR/doitGROMACS_xpm2txt.pl gibbs-12.xpm
  GGplot; mv $nameprod"_pes.png" ..
  cd ..
} &> >(tee doitgromacs_pca.log) >&2

#------------
# description : The name says everything
# requirements: .tpr + .xtc
# NOTE: USES A BIGGER PROBE
gromSAS() {   
  # if -s and -f are not set prompt and exit
  checkFlags_t
  # create the directory
  if [ ! -d ./sas_$name1 ] ; then
     mkdir sas_$name1
  fi
  cd sas_$name1
  if [ -n "${gromacs_ver}" ]; then    # if gromacs 5
    $groPATH/$g_sas -s ../$tpr -f ../$trj -n ../$name1.ndx                    \
    -o $name1"_area.xvg" -or $name1"_resarea.xvg" -dt $optionDTsas            \
    -b $optionSTARTime  -probe $optionPROBE -surface Protein -output Protein
  modVim $name1"_area.xvg" ; modVim $name1"_resarea.xvg"
  else
   (echo "Protein"; echo "$optionSAS") | $groPATH/$g_sas -s ../$tpr -f ../$trj  \
      -o $name1"_area.xvg" -or $name1"_resarea.xvg"    \
      -dt $optionDTsas -b $optionSTARTime -probe $optionPROBE || checkExitCode
   # probe 0.7 nm perchè 16A -> 1.6 nm -> raggio 0.7 nm
  modVim $name1"_area.xvg" ; modVim $name1"_resarea.xvg"
  fi ; cd ..
} &> >(tee doitgromacs_sas.log) >&2

gromSAS-sites() {
  # if -s and -f are not set prompt and exit
  checkFlags_t;
  # create the directory
  if [ ! -d ./sas_$name1 ] ; then
     mkdir sas_$name1
  fi
  cd sas_$name1
  if [ -n "${gromacs_ver}" ]; then    # if gromacs 5
    (echo "G6P"; echo "Coenzyme", echo "strNADPplus" ) | $groPATH/$g_sas      \
      -s ../$tpr -f ../$trj -n ../$name1.ndx -o $name1"_area_sites.xvg"       \
      -or $name1"_resarea_sites.xvg" -dt $optionDTsas -b $optionSTARTime      \
      -surface Protein -output 
  modVim $name1"_area_sites.xvg" ; modVim $name1"_resarea_sites.xvg"
  else
    (echo "Protein"; echo "G6P") | $groPATH/$g_sas -s ../$tpr -f ../$trj -n ../$name1.ndx \
      -o $name1"_g6p_area.xvg" -or $name1"_g6p_resarea.xvg" -dt $optionDTsas              \
      -b $optionSTARTime  || checkExitCode 
    (echo "Protein"; echo "Coenzyme") | $groPATH/$g_sas -s ../$tpr -f ../$trj  \
      -n ../$name1.ndx -o $name1"_Coenzyme_area.xvg" -b $optionSTARTime        \
      -or $name1"_Coenzyme_resarea.xvg" -dt $optionDTsas || checkExitCode
    (echo "Protein"; echo "strNADPplus") | $groPATH/$g_sas -s ../$tpr -f ../$trj  \
      -n ../$name1.ndx -o $name1"_strNADPplus_area.xvg"  -b $optionSTARTime       \
      -or $name1"_strNADPplus_resarea.xvg" -dt $optionDTsas || checkExitCode
  fi ; cd ..
} &> >(tee doitgromacs_sasites.log) >&2

#-------
# description :
# requirements: .tpr + .xtc
gromHB() {
  checkFlags_t
  # create the directory
  if [ ! -d ./hydrogenBonds_$name1 ] ; then
     mkdir hydrogenBonds_$name1
  fi
  cd hydrogenBonds_$name1
  (echo "$optionHB"; echo "$optionHB") | $groPATH/$g_hbond -s ../$tpr -f ../$trj\
    -num $name1"_hb_count.xvg" -dist $name1"_hb_dist.xvg"                     \
    -hbm $name1"_hb_matrix" -tu ps -dt $optionDThb -b $optionSTARTime 
  modVim $name1"_hb_count.xvg" ; modVim $name1"_hb_dist.xvg"
  cd ..
} &> >(tee doitgromacs_hb.log) >&2

#-------
# description :
# requirements: .tpr + .xtc
gromHB-sites() {
  checkFlags_t;
  # create the directory
  if [ ! -d ./hydrogenBonds_$name1 ] ; then
    mkdir hydrogenBonds_$name1
  fi
  cd hydrogenBonds_$name1
  (echo "G6P"; echo "G6P") | $groPATH/$g_hbond -s ../$tpr                     \
    -f ../$trj -n ../$name1.ndx -tu ps -dt $optionDThb -b $optionSTARTime     \
    -num $name1"_G6P_count.xvg" -dist $name1"_G6P_dist.xvg"                   \
    -hbm $name1"_G6P_matrix" 
  modVim $name1"_G6P_count.xvg" ; modVim $name1"_G6P_dist.xvg"
  (echo "$Coenzyme"; echo "$Coenzyme") | $groPATH/$g_hbond -s ../$tpr         \
    -f ../$trj -n ../$name1.ndx -tu ps -dt $optionDThb -b $optionSTARTime     \
    -num $name1"_Coenzyme_count.xvg" -dist $name1"_Coenzyme_dist.xvg"         \
    -hbm $name1"_Coenzyme_matrix" 
  modVim $name1"_Coenzyme_count.xvg" ; modVim $name1"_Coenzyme_dist.xvg"
  (echo "strNADPplus"; echo "strNADPplus") | $groPATH/$g_hbond -s ../$tpr     \
   -f ../$trj -n ../$name1.ndx -tu ps -dt $optionDThb -b $optionSTARTime      \
  -num $name1"_strNADPplus_count.xvg" -dist $name1"_strNADPplus_dist.xvg"     \
  -hbm $name1"_strNADPplus_matrix"   
  modVim $name1"_strNADPplus_count.xvg" ; modVim $name1"_strNADPplus_dist.xvg"
  cd ..
} &> >(tee doitgromacs_hbsites.log) >&2

#---------
# description :
# requirements: .tpr + .xtc
gromDSSP() {
  # if -s and -f are not set prompt and exit
  checkFlags_t
  echo $optionDSSP | $groPATH/$do_dssp -ver 1 -f $trj -s $tpr -o $nameprod"_ss.xpm" \
    -sc $nameprod"_ss_count.xvg" -tu ps -dt $optionDTdssp -b $optionSTARTime  
  modVim $nameprod"_ss_count.xvg"
} &> >(tee doitgromacs_dssp.log) >&2

