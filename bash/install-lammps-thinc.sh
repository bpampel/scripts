# Install LAMMPS on Draco
# Last changed by Benjamin Pampel on Apr 08 2019
#

LAMMPS_DIR=`readlink -f $1`
PLUMED2_INSTALL_DIR="${HOME}/code/plumed2/source_thinc/install" # could also be set by some enviromental variable
PLUMED2_LINK_MODE="static"
CORES=4

VERSION=XE19u3
INITPATH=parallel_studio_xe_2019
INITSCRIPT=psxevars.sh
. /sw/linux/intel/${VERSION}/${INITPATH}/${INITSCRIPT} intel64 -python 2

# if shared check if plumed is installed at given dir
if [ ! -e "${PLUMED2_INSTALL_DIR}/lib/libplumed.a" ]; \
  then echo "Shared plumed library not found. Was plumed compiled and installed?"; \
  exit 1;
fi

cd "${LAMMPS_DIR}/src"
make lib-plumed args="-p ${PLUMED2_INSTALL_DIR} -m ${PLUMED2_LINK_MODE}"
make yes-class2 yes-kspace yes-manybody yes-molecule yes-partition yes-rigid yes-user-misc yes-user-plumed

cp ${LAMMPS_DIR}/Makefile.local MAKE
make -j $CORES mpi


echo "
#
# Added in LAMMPS installation on $(date)
# ------------------------------------------------------------
export PATH=${LAMMPS_DIR}/src:\${PATH}
# ------------------------------------------------------------
#
" > ~/source_lammps_thinc.sh
