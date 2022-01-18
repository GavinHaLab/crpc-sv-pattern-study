#!/bin/bash
# For covariate data, if there is any gap in between each bin, e.g. repliation
# timing data lifted over from hg19, it is necessary to close those gaps.
# This scripts detects gaps between every two bins and use the mean of two
# neighbors for the gap.
if [ $# -lt 2 ]
then
  echo "Usage: $0 <input bedGraph> <output bedGraph>"
  exit
fi
cat $1 \
  | awk 'BEGIN{OFS="\t"; pc=""; pe=0; pv=0;}
      {
        if ($1==pc) {
          if (pe < $2) {
            print pc, pe, $2, 0.5*(pv+$4);
            print;
          } else if (pe > $2) {
            print $1, pe, $3, $4;
          } else {
            print;
          }
        } else {
          #print $1, 0, $2, $4;
          print;
        }
        pc=$1; pe=$3; pv=$4;
      }' \
  > $2
