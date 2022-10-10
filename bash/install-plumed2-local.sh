# Install PLUMED 2 locally on pckr
# Last changed by Benjamin Pampel on Mar 06 2019
#
# Notes:
#  - This uses intel MPI and gcc with internal blas/LAPACK
#

CORES=4


PLUMED_DIR=`readlink -f $1`
INSTALL_DIR=${PLUMED_DIR}/install

# check if dir exists
if [ ! -e "${PLUMED_DIR}" ]; \
  then echo "Please specify the root directory of the plumed code as first argument"; \
  exit 1;
fi

cd $PLUMED_DIR
mkdir $INSTALL_DIR 2> /dev/null

export CXX=mpicxx
export CC=mpicc
export FC=mpifort

./configure --enable-modules=all --prefix=${PLUMED_DIR}/install

make -j $CORES
cd src/lib
make install


echo "
#
# Added in PLUMED 2 installation on $(date)
# ------------------------------------------------------------
export PLUMED2_INSTALL=${INSTALL_DIR}
export PATH=\${PLUMED2_INSTALL}/bin:\${PATH}
export INCLUDE=\${PLUMED2_INSTALL}/include:\${INCLUDE}
export CPATH=\${PLUMED2_INSTALL}/include:\${CPATH}
export LIBRARY_PATH=\${PLUMED2_INSTALL}/lib:\${LIBRARY_PATH}
export LD_LIBRARY_PATH=\${PLUMED2_INSTALL}/lib:\${LD_LIBRARY_PATH}
export PKG_CONFIG_PATH=\${PLUMED2_INSTALL}/lib/pkgconfig:\${PKG_CONFIG_PATH}
export PLUMED_KERNEL=\${PLUMED2_INSTALL}/lib/libplumedKernel.so
# ------------------------------------------------------------
#
_plumed() { eval "$(plumed --no-mpi completion 2>/dev/null)";}
complete -F _plumed -o default plumed
" > $INSTALL_DIR/source_plumed.sh
