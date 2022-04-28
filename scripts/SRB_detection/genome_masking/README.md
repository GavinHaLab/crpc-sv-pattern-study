# Universal mask for native hg38
## Genome FASTA file
File name: `GCA_000001405.15_GRCh38_no_alt_analysis_set.fna.gz`
It comes from the following blog:
https://lh3.github.io/2017/11/13/which-human-reference-genome-to-use

Download link:
ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/001/405/GCA_000001405.15_GRCh38/seqs_for_alignment_pipelines.ucsc_ids/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna.gz

For convenience, file was renamed to `hs38DH.fa`.

## RepeatMasker file
Downloaded from UCSC table browser hg38. File was saved as `rmsk.hg38.txt.gz`

## Tools
1. mdust
https://github.com/lh3/mdust
2. seqtk v1.3
https://github.com/lh3/seqtk/tree/v1.3
3. sort-alt
https://github.com/lh3/foreign/tree/master/sort

## Command line
Adopted from:
https://gist.github.com/lh3/9d6dcfc3436a735ef197
Also see the original paper, Suppl. Info. Section 4:
https://static-content.springer.com/esm/art%3A10.1038%2Fnature18964/MediaObjects/41586_2016_BFnature18964_MOESM204_ESM.pdf

### Compositional mask
1. low complexity by mdust
```
mdust hs38DH.fa -c -w7 -v28 \
  | cut -f1,3,4 \
  | gawk -vOFS="\t" '{--$2;print}' \
  | bgzip > hs38DH.mdust-w7-v28.bed.gz
```

2. long homopolymer by seqtk
```
seqtk hrun hs38DH.fa \
  | bgzip > hs38DH.hrun.bed.gz
```

3. satellite sequence by rmsk
```
zcat rmsk.txt.gz \
  | grep Satellite \
  | cut -f6,7,8 \
  | sort-alt -k1,1N -k2,2n \
  | bgzip > hs38DH.satellite.bed.gz
```

4. low complexity by rmsk
```
zcat rmsk.txt.gz \
  | grep Low_complexity \
  | cut -f6,7,8 \
  | sort-alt -k1,1N -k2,2n \
  | bgzip \
  > hs38DH.rmsk-lc.bed.gz
```

5. Rename chromosomes to make them consistent
ALT and FIX contigs in rmsk.txt are excluded.
```
for i in hs38DH.*bed.gz
do
  echo $i
  zcat $i \
    | sed 's/^chr//; s/v\([0-9]\)/.\1/; s/.\+_\(.\+\)_random/\1/' \
    | egrep -v "_alt|_fix" \
    > ${i%.*}
done
```

6. Merge these 4 files
```
cat hs38DH.hrun.bed hs38DH.mdust-w7-v28.bed hs38DH.rmsk-lc.bed hs38DH.satellite.bed \
  | sort-alt -k1,1N -k2,2n \
  | gawk -v OFS="\t" '{$2=$2<10?0:$2-10;$3+=10;print}' \
  | cut -f 1-3 \
  | mergeBed \
  | bgzip \
  > compositional.bed.gz
```

### Mappability mask
See http://lh3lh3.users.sourceforge.net/snpable.shtml

1. Extract all overlapping 75-mers
```
splitfa ../hs38DH.fa 75 \
  | split -l 20000000
```

2. Align all 75-mers to the genome
```
for i in x??
do
  bwa aln -R 1000000 -O 3 -E 3 hs38DH.fa $i \
    > ${i}.sai
  bwa samse $REF ${i}.sai $i \
    | bgzip \
    > ${i}.sam.gz
done
```

3. Generate raw mask and the final mask
```
gzip -dc x??.sam.gz \
  | ./gen_raw_mask.pl \
  > rawMask_75.fa
./gen_mask -l 75 -r 0.5 rawMask_75.fa \
  > mask_75_50.fa
```

4. Convert the final mask to BED files
```
./mask_to_bed.py mask_75_50.fa \
  > um_75_50_levels.bed
cat um_75_50_levels.bed \
  | awk '$5<2' \
  | sed 's/^chr//' \
  | cut -f 1-3 \
  | bgzip \
  > mappability.bed.gz
```

According to Heng's blog, levels have the following meanings:
- Level 3: the majortiy of overlapping k-mers are mapped uniquely and without
  1-mismatch (or 1-difference, depending on the BWA command line) hits.
- Level 2: the majority of overlapping k-mers are unique and level is not 3.
- Level 1: the majority of overlapping k-mers are non-unique.
- Level 0: all the k-mers overlapping x cannot be mapped due to excessive
  ambiguous bases.
In the SGDP paper (Mallick, et al. 2016), Level 2 and 3 are used to define the
high mappability region.

## Final mask
Merge the mask files to produce a final universal mask for the hs38DH genome.
```
cat <(zcat compositional.bed.gz) \
  <(zcat mappability.bed.gz) \
  | sort -k1,1 -k2,2n \
  | mergeBed \
  > um75-hs38DH.bed
complementBed -i \
  <(cat um75-hs38DH.bed \
    | awk 'BEGIN{OFS="\t"}$1!~/_/ && $1~/^[0-9XYM]/{print "chr"$0}') \
  -g hg38.chrom.sizes \
  | sed 's/chr//' \
  | awk '$1!~/_/ && $1~/^[0-9X]/' \
  > um75-hs38DH.complement.bed
```

Use the intersection of the final mask file with callable regions identified
by GATK `CallableLoci` (file not included) to get the eligible region file
`eligible_intersect_callable.bed` used for SRB analysis.
