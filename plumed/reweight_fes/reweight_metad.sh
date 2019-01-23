#!/bin/bash

plumedfile='/home/theorie/pampel/scripts/plumed/reweight_fes/plumed.dat'

refdir="/usr/data/pampel/averaging/chebyshev_ord20/1"
cd $refdir
feslist=`ls fes*`

# first create the basic structure
rootdir=$1
cd $rootdir
echo 'Currently working in '`pwd`
mkdir 'reweighting'
cd 'reweighting'
for fes in $( ls ../fes* )
do
  time=$(echo $fes | sed -s 's/\([0-9]*\).*-\([0-9]*\).*/\2/')
  # create dir and move files
  mkdir $time
  cd $time
  last_line=`expr $time \* 5 + 3`
  sed "$last_line"',$d' ../../COLVAR > COLVAR
  cp $plumedfile .
  # run plumed and move output to parent directory
  plumed --no-mpi driver --noatoms > /dev/null
  cp 'fes' $fes # variable already contains ../
  cd ../
  rm -rf $time
  # show progress
  ((counter+=1))
  echo -ne $((counter * 100 / filenum)) '% done\r'
done
echo ''
