#!/bin/bash
#$ -pe PE_8 8
#$ -m e
#$ -M pampel@mpip-mainz.mpg.de
#$ -S /bin/bash

############################################################################
# Definition of variables
############################################################################
LAMMPS=/people/thnfs/homes/pampel/opt/bin/lmp_mpi_THINC
MPI=/usr/lib64/mpi/gcc/openmpi/bin/mpirun
############################################################################

export LD_LIBRARY_PATH=/usr/lib64/mpi/gcc/openmpi/lib64/libmpi_cxx.so:$LD_LIBRARY_PATH

printf "Starting...\n"
date

printf "This job runs on $HOSTNAME\n"

WRITEDIR=/usr/scratch/pampel

cd /usr/scratch/

if [ -d pampel ]
then
    :
else
    mkdir pampel
fi

cd pampel

# copy needed files here (use rsync to avoid copying twice)
rsync -az /people/thnfs/homes/pampel/bin/Example2/ .

${LAMMPS} < start.lmp > out.lmp

rsync -az ./* /people/thnfs/homes/pampel/bin/Example2/

# uncomment to automatically remove all files in the scratch folder
cd ../
rm -rf pampel/

printf "Job done!\n"
date
