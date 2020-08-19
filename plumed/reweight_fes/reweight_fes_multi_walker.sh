#!/bin/bash

rw_dir='reweighting_discrete'
filenames1='fes.b1.iter-'
filenames2='.data'

first=10000
last=500000
step=10000

walkers=25

if [ "$#" -ne 2 ]; then
  echo "Wrong number of arguments.
        Usage:
        reweight_fes.sh $scriptfile $directory"
fi


plumedfile=$(realpath $1)


# first create the basic structure
rootdir=$2
cd $rootdir
echo 'Currently working in '`pwd`
mkdir "$rw_dir"
cd "$rw_dir"
for time in $(seq ${first} ${step} ${last})
do
  # create dir and concatenate COLVAR files until current time
  mkdir "$time"
  cd "$time"
  echo "Processing time $time"
  last_line=`expr $time \/ 10 + 2`
  head -n "$last_line" ../../COLVAR.0 > colvar.data
  for i in $(seq 0 $(expr $walkers - 1)) 
  do
    sed -n "2,${last_line}p" ../../COLVAR.$i >> colvar.data
  done
  cp $plumedfile plumed.dat
  # run plumed and move output to parent directory
  plumed --no-mpi driver --noatoms > /dev/null
  cp 'fes' ../${filenames1}${time}$filenames2
  cd ../
  rm -r $time
  echo "FES at time $time created.\n"
done
echo ''
