# Install LAMMPS locally on pckr (very basic)

STARTDIR=$(pwd)
LAMMPS_DIR="/usr/data/pampel/lammps/lammps-24Mar2022/"
PLUMED2_INSTALL_DIR="/usr/data/pampel/plumed-code/v2.8/install" # could also be set by some enviromental variable
PLUMED2_LINK_MODE="shared"
CORES=6

module load vmd

# if shared check if plumed is installed at given dir
#if [ ! -e "${PLUMED2_INSTALL_DIR}/lib/libplumed.a" ]; \
  #then echo "Shared plumed library not found. Was plumed compiled and installed?"; \
  #exit 1;
#fi

cd "${LAMMPS_DIR}/src"
make lib-plumed args="-p ${PLUMED2_INSTALL_DIR} -m ${PLUMED2_LINK_MODE}"
make yes-class2 yes-kspace yes-manybody yes-molecule yes-rigid yes-misc yes-plumed yes-molfile yes-extra-fix
make -j $CORES mpi


echo "
#
# Added in LAMMPS installation on $(date)
# ------------------------------------------------------------
export PATH=${LAMMPS_DIR}/src:\${PATH}
# ------------------------------------------------------------
#
" > $STARTDIR/source_lammps.sh
