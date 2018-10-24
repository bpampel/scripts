#!/bin/bash

plumedfile='/home/theorie/pampel/scripts/plumed/reweight_fes/plumed.dat'

# first create the basic structure
rootdir=$1
cd $rootdir
echo 'working in '`pwd`
mkdir 'reweighting'
cd 'reweighting'
for time in {0..20000..100}
do
  # to show progress because it takes a bit
  progress=`expr $time \* 100 \/ 20000`
  echo -ne $progress' % done\r'
  # create dir and move files
  mkdir $time
  cd $time
  last_line=`expr $time \* 5 + 3`
  sed "$last_line"',$d' ../../COLVAR > COLVAR
  cp $plumedfile .
  # run plumed and move output to parent directory
  plumed --no-mpi driver --noatoms > /dev/null
  filename='fes.b1.iter-'$time'.data'
  cp 'fes' '../'$filename
  cd ../
  rm -rf $time
done
echo ''
