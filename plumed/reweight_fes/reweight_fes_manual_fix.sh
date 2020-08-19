#!/bin/bash

# has to be set correctly to match the FES output stride
timestep=500
filenames1='fes.b1.iter-'
filenames2='.data'

first=62500
last=97500
step=2500


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
#filenum=$(ls fes* | wc -l) # just to monitor progress
#counter=0
#mkdir 'reweighting'
cd 'reweighting'
for time in $(seq ${first} ${step} ${last})
do
  # create dir and move files
  mkdir $time
  cd $time
  last_line=`expr $time \* $timestep + 3`
  sed "$last_line"',$d' ../../colvar.data > colvar.data
  cp $plumedfile plumed.dat
  # run plumed and move output to parent directory
  plumed --no-mpi driver --noatoms > /dev/null
  cp 'fes' ../${filenames1}${time}$filenames2
  cd ../
  #rm -rf $time
  echo "fes at time $time finished"
done
echo ''
