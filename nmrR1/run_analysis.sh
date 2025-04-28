#!/bin/bash

# get the number of frames from the number of lines in the
# file with box dimensions
line="$(wc -l ../box2.txt)"
N="$(echo $line | cut -d ' ' -f 1)"

# go through each frame and calculate the directions of the CH
# vectors defined in ch_vectors.txt
for (( k=0; k<$N; k+=1 ))
do
  echo $k > cf.txt
  tr -d '\n' < cf.txt > curr_frame.txt  
  /usr/local/bin/vmd -dispdev text -e get_CH_vectors.tcl ../struct.psf
done


