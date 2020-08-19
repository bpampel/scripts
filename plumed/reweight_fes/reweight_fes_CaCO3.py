#!/bin/bash

rw_dir='reweighting_long'
filenames1='fes.b1.iter-'
filenames2='.data'

iter_per_colvar=10 # note that iterations are more frequent than colvar printouts here!
first_line=2002  # first line of colvar files to consider

# times (iterations) to create
first_iter=30000
last_iter=300000
step_iter=1000

walkers=25

if [ "$#" -ne 2 ]; then
  echo 'Wrong number of arguments.
        Usage:
        reweight_fes.sh $scriptfile $directory'
    exit -1
fi


plumedfile=$(realpath $1)


# first create the basic structure
rootdir=$2
cd $rootdir
echo 'Currently working in '`pwd`
mkdir "$rw_dir" || exit -1
cd "$rw_dir"
for time in $(seq ${first_iter} ${step_iter} ${last_iter})
do
  # create dir and concatenate COLVAR files until current time
  mkdir "$time" || exit -1
  cd "$time"
  echo "Processing time $time"
  last_line=$(expr ${time} \/ ${iter_per_colvar} + 2)
  head -n 1 ../../COLVAR.0 > colvar.data
  for i in $(seq 0 $(expr $walkers - 1)) 
  do
    sed -n "${first_line},${last_line}p" ../../COLVAR.$i >> colvar.data
  done
  cp $plumedfile plumed.dat
  # run plumed and move output to parent directory
  plumed --no-mpi driver --noatoms > /dev/null
  cp 'fes' ../${filenames1}${time}$filenames2
  cd ../
  rm -r $time
  echo "FES at time $time created."
done
echo ''
