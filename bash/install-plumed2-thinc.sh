#!/usr/bin/env bash
# Authors: Bin Song

# Install PLUMED2 on MPIP local cluster
# intel compiler setup. Use -python 3 for python3

CORES=6

PLUMED_DIR=`readlink -f $1`
INSTALL_DIR=${PLUMED_DIR}/install

VERSION=XE19u3
INITPATH=parallel_studio_xe_2019
INITSCRIPT=psxevars.sh
. /sw/linux/intel/${VERSION}/${INITPATH}/${INITSCRIPT} intel64 -python 2


# check if dir exists
if [ ! -e "${PLUMED_DIR}" ]; \
  then echo "Please specify the root directory of the plumed code as first argument"; \
  exit 1;
fi

cd $PLUMED_DIR
mkdir $INSTALL_DIR 2> /dev/null

export CXX=mpiicpc
export CC=mpiicc
export FC=mpiifort

## -fp-model precise -fma -align -finline-functions are optional optimization flags that increases your time of compilation. You are strongly advised to use -no-prec-div with intel compilers.
## Intel compilers support single excutable for multiple code paths. Hence, multiple SIMD sets: avx, avx2 and avx512 can be used.
## Parallel16.q and Debian20.q nodes support avx2 and below; Parallel20-3.q nodes also support avx512 and below.

./configure CXXFLAGS="-fp-model precise -fma -align -finline-functions -no-prec-div -xcore-avx2 -axcore-avx512 -qopt-zmm-usage=high" CPPFLAGS="-I${MKLROOT}/include -I${MKLROOT}/include/fftw" LDFLAGS="-L${MKLROOT}/lib/intel64" LIBS=-mkl --enable-fftw=yes --enable-modules=ves --enable-boost_serialization --prefix=${INSTALL_DIR}

make -j4
cd src/lib
make install
cd -

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
" >> ~/source_plumed_thinc.sh
