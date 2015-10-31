#!/bin/bash

# DESCRIPTION: Extra functions, 

mean_sas() {
  ls && echo "select some files please: " &&  read -a FILES_ARRAY
  for i in "${FILES_ARRAY[@]}" ; do
    modVim $i
    # grep: avoid commented lines
    # awk: select 4th column and calculate the mean
    grep -v '^#' $i | awk '{total += $4} END {print total/NR}'
  done
}

mean_hb() {
  ls && echo "select some files please: " &&  read -a FILES_ARRAY
  for i in "${FILES_ARRAY[@]}" ; do
    modVim $i
    # grep: avoid commented lines
    # awk: select 2nd column and calculate the mean
    grep -v '^#' $i | awk '{total += $2} END {print total/NR}'
  done
}

CAtoMAIN () {
  filenameFULL=$(basename $1)
  filename="${filenameFULL%.*}"
  echo $filename
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

