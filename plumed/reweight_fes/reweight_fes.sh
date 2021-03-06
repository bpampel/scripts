#!/bin/bash

# has to be set correctly to match the FES output stride
timestep=5


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
filenum=$(ls fes* | wc -l) # just to monitor progress
counter=0
mkdir 'reweighting'
cd 'reweighting'
for fes in $( ls ../fes* )
do
  time=$(echo $fes | sed -s 's/\([0-9]*\).*-\([0-9]*\).*/\2/')
  # create dir and move files
  mkdir $time
  cd $time
  last_line=`expr $time \* $timestep + 3`
  sed "$last_line"',$d' ../../colvar.data > colvar.data
  cp $plumedfile plumed.dat
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
