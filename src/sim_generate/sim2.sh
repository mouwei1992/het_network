#! /bin/bash

for i in `seq 1 10`;
  do
    ./sim2_no_proc.out
    ./sim2_per10.out
    ./sim2_per100.out
    ./sim2_per1000.out
  done
