# Install LAMMPS on Draco
# Last changed by Benjamin Pampel on Apr 08 2019
#
# Notes:
#  - this links lammps with a static plumed library which has to be generated first
#

LAMMPS_DIR="${PWD}/lammps"
PLUMED2_INSTALL_DIR="${HOME}/code/plumed2/br-ves-localized-bfs/install" # could also be set by some enviromental variable
PLUMED2_LINK_MODE="static"
CORES=4

module load vmd

# if shared check if plumed is installed at given dir
#if [ ! -e "${PLUMED2_INSTALL_DIR}/lib/libplumed.a" ]; \
  #then echo "Shared plumed library not found. Was plumed compiled and installed?"; \
  #exit 1;
#fi

cd "${LAMMPS_DIR}/src"
make lib-plumed args="-p ${PLUMED2_INSTALL_DIR} -m ${PLUMED2_LINK_MODE}"
make yes-class2 yes-kspace yes-manybody yes-molecule yes-rigid yes-user-misc yes-user-plumed yes-user-molfile
make -j $CORES mpi


echo "
#
# Added in LAMMPS installation on $(date)
# ------------------------------------------------------------
export PATH=${LAMMPS_DIR}/src:\${PATH}
# ------------------------------------------------------------
#
" > ~/source_lammps.sh


