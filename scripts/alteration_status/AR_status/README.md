# Produce CRPC10X42_WCDT_AR-enh_compare_AR-vs-Enh.txt
Rscript compareCN_2regions_v2.R pathToTitanOptimalClusterSolution/ ../../../metadata/AR_coord.txt ../../../metadata/sampleList_allCases.txt CRPC10X42_WCDT_AR-enh_compare

# Fig 3-B
Rscript AR_expression_with_copy_number_status.R
