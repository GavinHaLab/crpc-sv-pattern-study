#!/bin/bash
TELOMERE_REF=telomere_q5m.bed
CHROMOTHRIPSIS_REF=chromothripsis.bed

set -euo pipefail
echo "[$(date) | preparing input files]"
cat input/header_sv \
  <(tail -n +2 input/sv_data.bedpe \
    | sort-alt -k12,12 -k1,1N -k2,2n \
    | awk -F "\t" 'BEGIN{OFS="\t"}{if($24 >= 1e4 || $24 == -1 \
          || $47 == "Inversion-FoldBack") \
        {print $12":"$1, $2, $3, NR, 0, $9;
         print $12":"$4, $5, $6, NR, 1, $10}}' \
    | intersectBed -a stdin -b $CHROMOTHRIPSIS_REF -v \
    | sort -k4,4n -k5,5n \
    | awk -F "\t" 'BEGIN{OFS="\t"; pid=0; pchr=""; pstart=0; pstrand="";} \
        {
          if(pid==$4) {
            split($1, f, ":");
            split(pchr, g, ":");
            print g[1], $4, (g[2]=="chrX")?23:g[2], pstart, (pstrand=="+")?0:1,
              (f[2]=="chrX")?23:f[2], $2, ($6=="+")?0:1;
          }
          pid=$4;
          pchr=$1;
          pstart=$2;
          pstrand=$6;
        }' \
    | sed 's/chr//g' \
    | awk -F "\t" 'BEGIN{OFS="\t"}{print $3, $4-1, $4, $0}' \
    | intersectBed -a stdin -b $TELOMERE_REF -v \
    | cut -f 4- \
    | awk 'BEGIN{OFS="\t"}{print $6, $7-1, $7, $0}' \
    | intersectBed -a stdin -b $TELOMERE_REF -v \
    | cut -f 4-) \
  | ./filter_blacklist.py blacklist.txt \
  > sv_merged.txt

Rscript prepare_seg_from_TITAN.r -n -i cn_data.seg -o cn_data.tmp
cat cn_data.tmp \
  | sed 's/chr\([^\t]\+\)/\1/g' \
  | awk 'BEGIN{OFS="\t"}{if($2=="X")$2=23; print}' \
  | awk 'BEGIN{OFS="\t"}{if(NR>1){$2="chr"$2}print}' \
  | ./filter_blacklist.py blacklist.txt \
  | sed 's/chr\([^\t]\)/\1/' \
  > cn_merged.seg

echo "[$(date) | starting ChainFinder]"
cd bin
# run_ChainFinder.sh is included in ChainFinder
./run_ChainFinder.sh $PATH_TO_MATLAB
