#!/bin/bash

print_usage() {
  echo 'Usage: batch_simulations -d dir [-s first_iteration] -e last_iteration'
}

# set variables from flags
START_ITER=1
while getopts ':d:s:e:' flag; do
  case "${flag}" in
    d) DIR="${OPTARG}" ;;
    s) START_ITER="${OPTARG}" ;;
    e) END_ITER="${OPTARG}" ;;
    *) print_usage
       exit 1 ;;
  esac
done

# check if arguments were given
if [ -z "$DIR" ]; then echo "Please specify directory"; print_usage; exit 1; fi
if [ -z "$END_ITER" ]; then echo "Please specify iterations"; print_usage; exit 1; fi

# actual execution
for ITER in $(eval echo "{${START_ITER}..${END_ITER}}"); do
  sed -i 's/-J\ ves_.*/-J\ ves_'"$DIR"'_'"$ITER"'/' 'run_simulation.sh'
  sbatch run_simulation.sh $DIR $ITER
done
