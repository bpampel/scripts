#!/bin/bash

filenames1='fes.b1.iter-'
filenames2='.data'

iter_per_colvar=10 # note that iterations are more frequent than colvar printouts here!
first_line=2002  # first line of colvar files to consider

# times (iterations) to create
first_iter=30000
last_iter=500000
step_iter=10000

walkers=25

if [ "$#" -ne 3 ]; then
  echo 'Wrong number of arguments.
        Usage:
        reweight_fes.sh $scriptfile $directory $num_traj'
    exit -1
fi


plumedfile=$(realpath $1)
num_traj=$3

rw_dir="reweighting_individual_${num_traj}"

# first create the basic structure
rootdir=$2
cd $rootdir
echo 'Currently working in '`pwd`
mkdir "$rw_dir" || exit -1
cd "$rw_dir"
for i in $(seq 0 ${num_traj} 24)
do
    final_walker=$(expr $i + $num_traj - 1)
    if [ $final_walker -ge 24 ]; then break; fi
    mkdir "walker${i}" || exit -1
    cd "walker${i}"
    for time in $(seq ${first_iter} ${step_iter} ${last_iter})
    do
      # create dir and concatenate COLVAR files until current time
      mkdir "$time" || exit -1
      cd "$time"
      echo "Processing time $time"
      last_line=$(expr ${time} \/ ${iter_per_colvar} + 2)
      head -n 1 ../../../COLVAR.${i} > colvar.data
      for j in $(seq $i $final_walker)
      do
          sed -n "${first_line},${last_line}p" ../../../COLVAR.${j} >> colvar.data
      done
      cp $plumedfile plumed.dat
      # run plumed and move output to parent directory
      plumed --no-mpi driver --noatoms > /dev/null
      cp 'fes' ../${filenames1}${time}$filenames2
      cd ../
      rm -r $time
      echo "FES at time $time created."
    done
    /home/theorie/pampel/scripts/python/delta_F.py . -d -kT 0.025852 -d -m /usr/data/pampel/CaCO3/mask* -c 2
    echo "delta_F file for walker $i created"
    cd ../
done
echo ''
