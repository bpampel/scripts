#/bin/bash

# this script is used to sparsify colvar files of sets of simulations (for smaller storage size)
# the folders are numbered 1-20, and only every 100th point of the colvar files is kept

for i in {1..20}; do
  cd $i
  awk '(!((NR-2)%100)) || (NR<2)' colvar.data > colvar_sparse.data
  mv colvar_sparse.data colvar.data
  cd ../
done
