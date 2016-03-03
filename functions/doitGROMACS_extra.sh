#!/bin/bash

# DESCRIPTION: Extra functions, 

mean_single() {
  ls && echo "select some files please: " &&  read -a FILES_ARRAY
  for i in "${FILES_ARRAY[@]}" ; do
    modVim $i
    # grep: avoid commented lines
    # awk: select 4th column and calculate the mean
    grep -v '^#' $i | awk '{total += $2} END {print total/NR}'
  done
}

mean_multi() {
  ls && echo "select some files please: " &&  read -a FILES_ARRAY
  for i in "${FILES_ARRAY[@]}" ; do
    modVim $i
    # grep: avoid commented lines
    # awk: select 2nd column and calculate the mean
    grep -v '^#' $i | awk '{total += $2} END {print total/NR}'
    grep -v '^#' $i | awk '{total += $3} END {print total/NR}'
    grep -v '^#' $i | awk '{total += $4} END {print total/NR}'
    grep -v '^#' $i | awk '{total += $5} END {print total/NR}'
  done
}

PDBtoXTC() {
  # input is passed using the -c flag
  filenameFULL=$(basename $pdb1)
  filename="${filenameFULL%.*}"
  $groPATH/$trjconv -f $pdb1 -o $filename".xtc"
  $groPATH/gmxcheck -f $filename".xtc" 
  echo "done, check if the number of frame matches."

  echo Protein | $groPATH/$trjconv -s $optionUNRESfilename -f $filename".xtc" -o $filename".gro"       \
    -dump $optiondump || checkExitCode
   # account for the periodicity (nojump)
  echo Protein | $groPATH/$trjconv -s $filename".gro" -f $filename".xtc"        \
    -o $filename"_nojump.xtc" -pbc nojump || checkExitCode
  (echo "Backbone"; echo "Protein") | $groPATH/$trjconv -s $filename".gro" -f $filename"_nojump.xtc" \
    -o $filename"_fit.xtc" -fit rot+trans || checkExitCode
}&> >(tee doitgromacs_PDBtoXTC.log) >&2

CAtoMAIN() {
  filenameFULL=$(basename $1)
  filename="${filenameFULL%.*}"
  grep -B1000000 -m1 -e "TER" $1 > start.pdb
  ~martin/bin/catomain start.pdb $filename"_start.pdb"
  grep -A1000000 -m2 -e "TER" $1 > end.pdb
  sed -e '/^TER/d' end.pdb
  ~martin/bin/catomain end.pdb $filename"_end.pdb"
  sed -i 's/ A / B /g' $filename"_end.pdb" 

  cat $filename"_start.pdb" $filename"_end.pdb" > $filename"_all.pdb"
  rm start.pdb; rm end.pdb
  rm $filename"_start.pdb"; rm $filename"_end.pdb"
}&> >(tee doitgromacs_catomain.log) >&2

split_states() {
  if [ ! -d ./$1"_states" ]; then
    mkdir $1_states
  fi
  cp $1 $1"_states"
  cd $1"_states"
  touch catomain_states.pdb # file containing all the pdbs together

  # grep the times and save them in a list
  #times=($(grep '^REMARK' $1 | grep -Eo '[0-9]+.[0-9]+\s'))
  #lines=($(grep -n 'REMARK' $1 | grep -Eo '^[0-9]+'))
  #length_times=${#times[@]}

  IFS="
"
  structures=($(grep '^REMARK' $1))   # save lines in array
  tot_structures=${#structures[@]}    # return total elements of array
  nbr_structures=($(grep '^REMARK' $1 | grep -Eo ' ([0-9].) '))
  
  for (( i=0; i<$((tot_structures-1)); i++ )); do
    echo ${structures[i]} >> catomain_states.pdb
    
    if [[ $i  == $(($tot_structures - 1)) ]]; then
      grep -A1000000 -m1 -e "${structures[i]}" $1 > "$i.pdb"
      CAtoMAIN "$i.pdb"
      cat $i"_all.pdb" >> catomain_states.pdb
      echo "ENDMDL" >> catomain_states.pdb 
    else
      grep -A100000 "${structures[i]}" $1 | grep -B1000000 "${structures[i+1]}" > "$i.pdb"
      sed -i '$ d' "$i.pdb"
      CAtoMAIN "$i.pdb"
      cat $i"_all.pdb" >> catomain_states.pdb
      echo "ENDMDL" >> catomain_states.pdb
    fi
  done
  rm $1
  cd ..
}&> >(tee doitgromacs_split_states.log) >&2

# version without the final concatenation (catomain_states.pdb)
split_statesORIGINAL() {
  if [ ! -d ./$1"_states" ]; then
    mkdir $1_states
    cp $1 $1"_states"
    cd $1"_states"
  fi

  # check how many structure you have
  #grep '^REMARK' $1

  # grep the times and save them in a list
  times=($(grep '^REMARK' $1 | grep -Eo '[0-9]+.[0-9]+\s'))
  #lines=($(grep -n 'REMARK' $1 | grep -Eo '^[0-9]+'))
  length_times=${#times[@]}

  for (( i=0; i<$length_times; i++ )); do
    if [[ $i == $((length_times - 1)) ]]; then
      grep -A1000000 -m1 -e "REMARK.time\s*${times[$i]}" $1 > "${times[i]}.pdb"
      CAtoMAIN "${times[i]}.pdb"
    else
      grep -A100000 "REMARK.time\s*${times[$i]}" $1 | grep -B1000000 "REMARK.time\s*${times[$i+1]}" > "${times[i]}.pdb"
      sed -i '$ d' "${times[$i]}.pdb"
      CAtoMAIN "${times[i]}.pdb"
    fi  
  done

  rm $1
  cd ..
}&> >(tee doitgromacs_split_states.log) >&2

