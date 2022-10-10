#!/bin/bash -l
# Standard output and error:
#SBATCH -o ./tjob.out.%j
#SBATCH -e ./tjob.err.%j
# Initial working directory:
#SBATCH -D .
# Job Name:
#SBATCH -J ves_splines_N22_tdwt_20
# Queue (Partition):
#SBATCH --partition=short
# Number of nodes and tasks per node:
#SBATCH --ntasks=32
#SBATCH --ntasks-per-core=1
#SBATCH --mem=1000
#
#SBATCH --mail-type=end
#SBATCH --mail-user=pampel@mpip-mainz.mpg.de
#
# Wall clock limit:
#SBATCH --time=04:00:00


# create dir from command line argument
DIR=$1
ITER=$2
STARTDIR=/u/bpampel/${DIR}
TMPDIR=/ptmp/bpampel/${DIR}/${ITER}
echo "Starting simulation #${ITER} with plumed.dat from folder ${DIR}"
mkdir -p $TMPDIR
cd $TMPDIR

# move needed files here and remember them
rsync -az ${STARTDIR}/input/ ./
STARTFILES=`ls`
cp ${STARTDIR}/plumed.dat ./
echo 'Copied files to '`pwd`

# change seed using python random
SEED=`python -c "import random; print random.randint(10000000,99999999)"`
sed -i 's/seed\ equal.*/seed\ equal\ '"$SEED"'/' 'start.lmp'

# run job
srun ./lmp_mpi < start.lmp > out.lmp
echo 'Job finished, now syncing files'

# copy output files to userspace and clean up
rm $STARTFILES
cd ../
rsync -az $ITER /u/bpampel/$DIR/
rm -r $ITER
