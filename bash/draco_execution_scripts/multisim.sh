#!/bin/bash -l
# Standard output and error:
#SBATCH -o ./tjob.out.%j
#SBATCH -e ./tjob.err.%j
# Initial working directory:
#SBATCH -D .
# Job Name:
#SBATCH -J ves_puremd_1
# Queue (Partition):
#SBATCH --partition=short
# Number of nodes and tasks per node:
#SBATCH --ntasks=16
#SBATCH --ntasks-per-core=1
#SBATCH --mem=1000
#
#SBATCH --mail-type=end
#SBATCH --mail-user=pampel@mpip-mainz.mpg.de
#
# Wall clock limit:
#SBATCH --time=04:00:00


# create dir from command line argument
dir=$1
iter=$2
echo 'Starting iteration '$iter' with plumed.dat from folder '$dir
mkdir -p /ptmp/bpampel/$dir/$iter
cd /ptmp/bpampel/$dir/$iter

# move needed files here and remember them
rsync -az /u/bpampel/input/ ./
startfiles=`ls`
cp '/u/bpampel/'$dir'/plumed.dat' ./
echo 'Copied files to '`pwd`

# change seed using python random
seed=`python -c "import random; print random.randint(10000000,99999999)"`
sed -i 's/seed\ equal.*/seed\ equal\ '"$seed"'/' 'start.lmp'

# run job
srun ./lmp_mpi < start.lmp > out.lmp
echo 'Job finished, now syncing files'

# copy output files to userspace and clean up
rm $startfiles
cd ../
rsync -az $iter /u/bpampel/$dir/
rm -r $iter

# invoke next simulation or stop if iterations are reached
iter=`expr $iter - 1`
if (( iter==0 )); then
  cd ../
  rm -r $dir
  echo 'All iterations finished!'
else
  cd /u/bpampel/
  echo 'Executing next iteration'
  sed -i 's/-J\ ves_.*/-J\ ves_'"$dir"'_'"$iter"'/' 'multisim.sh'
  sbatch multisim.sh $dir $iter
fi
