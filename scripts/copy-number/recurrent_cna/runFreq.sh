# This will produce Figure S1 A and B
# run runTitanMatrix.sh before run this script.
# output RData (from runTitanMatrix.sh) is required as 1st argument.

Rscript getGeneFrequency_filterLength_groupTest_hg38.R CRPC10X_WCDT_100kb_geneMats.RData 0 ../../../metadata/100kb.bins.bed 0 chrPosn 1e7 overall fdr ../../../metadata/all_lesions.conf_99.txt ../../../metadata/combined_gene_list_159.txt CRPC10X_WCDT_100kb
