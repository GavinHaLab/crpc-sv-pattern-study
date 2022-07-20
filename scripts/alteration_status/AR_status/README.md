### 1. Produce "CRPC10X42_WCDT_AR-enh_compare_AR-vs-Enh.txt"

Rscript compareCN_2regions_v2.R {pathToTitanOptimalClusterSolution} ../../../metadata/AR_coord.txt ../../../metadata/sampleList_allCases.txt CRPC10X42_WCDT_AR-enh_compare

### 2. Create Fig 3-B 
Rscript AR_expression_with_copy_number_status.R

* This code takes "CRPC10X42_WCDT_AR-enh_compare_AR-vs-Enh.txt" from compareCN_2regions_v2.R.

