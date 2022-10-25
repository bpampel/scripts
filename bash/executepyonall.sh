#!/bin/sh
script=$1

for i in $(find . -name 'order*');
do
  dir="${i}/"
  cp 'Python_scripts/'$script $dir
  cd $dir
  './'$script
  rm $script
  cd ../../
done
