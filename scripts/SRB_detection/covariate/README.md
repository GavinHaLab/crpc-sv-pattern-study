1. Nucleotide composition
File name: `nt_composition.hg38.rds`. This file uses tiling non-overlapping
10kb bins. See `cov_nt_composition.r`

2. Replication timing
Two files, `reptiming_NHEK.hg19.bed` and `reptiming_LNCaP.hg19.bed` were used.
NHEK file is downloaded from the following link, under accession Int92817591:
https://www2.replicationdomain.com/database.php
LNCaP file is downloaded from the following link, under accession ENCFF995YGM:
https://www.encodeproject.org/replication-timing-series/ENCSR888VII/

Commands for lifting over and filtering:
```
IN=reptiming_LNCaP.hg19.bed
# some bins become too wide after liftover. Those might be problematic and
# thus we are throwing out those longer than 1.5x of the original bin size.
# For LNCaP, bin size is 5kb, so cutoff was set to 7.5kb.
# For NHEK, it uses various bin size, per-bin cutoff is needed but not
# implemented.
CUTOFF=7500
cat $IN \
  | awk 'BEGIN{OFS="\t"}{print $1, $2, $3, $1":"$2"-"$3":"$4}' > ${IN%.*}.data
liftOver ${IN%.*}.data hg19ToHg38.over.chain.gz \
  ${IN%.hg19*}.hg38.unfiltered.data ${IN%.hg19*}.hg38.unconverted.data
cat ${IN%.hg19*}.hg38.unfiltered.data \
  | awk -v c=$CUTOFF 'BEGIN{OFS="\t"}$3-$2<c{split($4, f, ":"); \
    print $1, $2, $3, f[3]}' \
  | sort -k1,1 -k2,2n -k3,3n \
  | ./filter_chr.sh \
  > ${IN%.hg19*}.hg38.filtered.bed
./smoothing.sh ${IN%.hg19*}.hg38.filtered.bed ${IN%.hg19*}.hg38.smoothed.bed
```

3. DNase hyper sensitive sites
DHS peaks are downloaded from ENCODE LNCaP, genome version hg38.
Sample accession numbers: ENCFF434GSJ (rep1) and ENCFF535OEM (rep2).
https://www.encodeproject.org/experiments/ENCSR000EPF/
Files are saved as `dhs_LNCaP_rep1.hg38.bed` and `dhs_LNCaP_rep2.hg38.bed`.

4. RepeatMasker
Five types are used: SINE, LINE, LTR, simple repeat, and DNA transposon.
Annotation is extracted from the UCSC genome browser version of RepeatMasker
for hg38.
```
# coordinates already 0-based
zcat rmsk.hg38.txt.gz \
  | awk 'BEGIN{OFS="\t"}$6 !~/_/ && $6 ~/chr[0-9X]/ && $12 == "LINE"\
    {print $6, $7, $8}' \
  | sort -k1,1 -k2,2n \
  | mergeBed \
  | ./filter_chr.sh \
  > retrotransposon.line.hg38.bed
zcat rmsk.hg38.txt.gz \
  | awk 'BEGIN{OFS="\t"}$6 !~/_/ && $6 ~/chr[0-9X]/ && $12 == "SINE"\
    {print $6, $7, $8}' \
  | sort -k1,1 -k2,2n \
  | mergeBed \
  | ./filter_chr.sh \
  > retrotransposon.sine.hg38.bed
zcat rmsk.hg38.txt.gz \
  | awk 'BEGIN{OFS="\t"}$6 !~/_/ && $6 ~/chr[0-9X]/ && $12 == "LTR"\
    {print $6, $7, $8}' \
  | sort -k1,1 -k2,2n \
  | mergeBed \
  | ./filter_chr.sh \
  > retrotransposon.ltr.hg38.bed
zcat rmsk.hg38.txt.gz \
  | awk 'BEGIN{OFS="\t"}$6 !~/_/ && $6 ~/chr[0-9X]/ && $12 == "Simple_repeat"\
    {print $6, $7, $8}' \
  | sort -k1,1 -k2,2n \
  | mergeBed \
  | ./filter_chr.sh \
  > retrotransposon.simple.hg38.bed
zcat rmsk.hg38.txt.gz \
  | awk 'BEGIN{OFS="\t"}$6 !~/_/ && $6 ~/chr[0-9X]/ && $12 == "DNA"\
    {print $6, $7, $8}' \
  | sort -k1,1 -k2,2n \
  | mergeBed \
  | ./filter_chr.sh \
  > retrotransposon.dna.hg38.bed
```

5. Chromatin states
Model is downloade from the Roadmap epigenetics project. Since H3K27ac is
available for LNCaP, the 18-state model trained from 6 chromatin marks was
used.
https://egg2.wustl.edu/roadmap/data/byFileType/chromhmmSegmentations/ChmmModels/core_K27ac/jointModel/final/model_18_core_K27ac.txt
```
# Binarize BAM files
java -mx4000M -jar ChromHMM.jar BinarizeBam hg38.chrom.sizes input cell_mark.txt binbam
java -mx4000M -jar ChromHMM.jar MakeSegmentation model_18_core_K27ac.txt binbam .

# Segment genome
java -mx4000M -jar ChromHMM.jar MakeSegmentation model_18_core_K27ac.txt binbam .

# Rename states
# trim output files for out-of-boundary bins and remove unwanted chromosomes
bedtools slop -i LNCaP_18_segments.bed \
    -g hg38.chrom.sizes -b 0 \
  | awk '$1!~/_/ && $1~/chr[0-9X]/' \
  > LNCaP_18_segments.clean.bed
# substitute state IDs
join -1 4 -2 1 -t '	' <(sort -k4,4 LNCaP_18_segments.clean.bed) \
    <(sort -k1,1 table_18states.txt) \
    -o 1.1,1.2,1.3,2.2 \
  | sort -k1,1 -k2,2n \
  > LNCaP_chromHMM_18states.clean.bed
```


6. Common fragile sites
Fragile sites are downloaded from [HGNC biomart](https://biomart.genenames.org/martform/#!/default/HGNC?datasets=hgnc_gene_mart).
Choose "fragile site" in "Locus type", and choose "Approved symbol",
"Chromosome location" to get the table `fragile_sites_table.txt`. All fragile
sites are annotated as cytoband.
*Note*:
1. There are two sites labeled as `withdrawn` by HGNC, and they are not
   included in the rest of analysis (FRA4E and FRA8F).
2. FRA10A is annotated as "10q23.3 or 10q24.2". In this analysis both are
   included.

Downloaded from the link below, saved as `common_fragile_sites.hg19.txt`
https://data.broadinstitute.org/snowman/ofer/1934ICGC_SV_RetrievedfromOferDesktop/gBreak/data-raw/common_fragile_sites.txt
This table includes common fragile sites collected from multiple studies, and
some of them were refined to more specific locations.

To combine these, two processes are taken:
1. Lift over fragile sites that have a citation. Those are refined.
2. Extract sites that only have a cytoband annotation. Those are directly taken
from hg38.
```
tail -n +2 common_fragile_sites.hg19.txt \
  | awk -F "\t" '$5!~/[0-9X][pq]/' \
  | cut -f 1-4 > cfs_refined.hg19
liftOver cfs_refined.hg19 \
  hg19ToHg38.over.chain.gz \
  cfs_refined.hg38 /dev/null
grep -v -f \
    <(cut -f 4 cfs_refined.hg19 | sed 's/[a-z]\+//g' | uniq) \
    fragile_sites_table.txt \
  > band_list.txt
./merge_bands.py hg38.cytoband.bed band_list.txt > cfs_band.hg38
cat cfs_refined.hg38 cfs_band.hg38 \
  | sed 's/chr//' \
  | sort -k1,1 -k2,2n \
  | cut -f 1-4 \
  > common_fragile_sites.hg38.txt
rm cfs_refined.hg19 cfs_refined.hg38 band_list.txt cfs_band.hg38
```
