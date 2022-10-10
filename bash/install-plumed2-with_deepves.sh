#!/usr/bin/env bash

# Install PLUMED2 with libtorch

# set the correct directories of plumed and libtorch
libtorch_root="/usr/data/pampel/plumed_ves_ml/libtorch"
plumed_root="/usr/data/pampel/plumed_ves_ml/plumed2"
# from plumed-nest 19/060
neural_network_ves_file="/usr/data/pampel/plumed_ves_ml/deep_ves-data/0_code/NeuralNetworkVes.cpp"

INSTALL_DIR=${PWD}/install
mkdir install

cd $plumed_root || exit 1
cp $neural_network_ves_file src/bias/ || exit 1

# all flags are to include libtorch
./configure CPPFLAGS="-I${libtorch_root}/include -I${libtorch_root}/include/torch/csrc/api/include" CXXFLAGS="-D_GLIBCXX_USE_CXX11_ABI=0 -std=c++14" LDFLAGS="-L${libtorch_root}/lib" --enable-modules=ves --prefix=${INSTALL_DIR}


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
" > source_plumed.sh
