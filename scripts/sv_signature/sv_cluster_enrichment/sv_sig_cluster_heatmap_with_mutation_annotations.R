#SV Signature Cluster Heatmap for WCDT

#Given a signature exposure file output from signature.tools.lib, a number of
#clusters, and annotation files, writes a file of the cluster each sample is
#assigned to with hierarchical clustering and plots a heatmap of the clusters
#with selected annotations.


###imports
library(tidyverse)
library(data.table)
library(knitr)
library(pheatmap)
library(dendextend)
library(tseries)
library(qdap)
library(grid)

options(stringsAsFactors = FALSE)


###settings
sig_file <- "../../../metadata/SigFit_withBootstrap_Exposures_mKLD_bfmCosSim_alpha-1_tr5_p0.05.tsv"
output_folder <- "./outputs/"

tdp_file <- "../../../metadata/all_TDP_status.csv"
tmb_file <- "../../../metadata/tmb_all_combined.csv"

chromothripsis_file <- "../../../metadata/all_Chromothripsis.csv"
focal_chromothripsis_file <- "../../../metadata/all_curated_chromothripsis_focal.csv"
bfb_chromothripsis_file <- "../../../metadata/all_curated_chromothripsis_bfb.csv"
Chromoplexy_file <- "../../../metadata/all_Chromoplexy.csv"

therapy_file <- "../../../metadata/all_therapy.csv"
treatment_file <- "../../../metadata/all_treatment.csv"
mutation_file <- "../../../metadata/all_mutation_data_combined.csv"

cna_file <- "../../../metadata/all_comut_CN_cf0.8.csv"

sv_file <- "../../../metadata/all_comut_SV_combined_single.csv"

biopsy_site_file <- "../../../metadata/all_biopsy_sites.csv"
cancer_type_file <- "../../../metadata/all_cancer_type.csv"

ploidy_file <- "../../../metadata/all_ploidy.csv"
purity_file <- "../../../metadata/all_purity.csv"

ets_curation_file <- "../../../metadata/combined_ets_all.csv"
ar_category_file <- "../../../metadata/ar_combined_comut_file.csv"

###functions
##change cluster name based on SV cluster order
cluster_swap <- function(vec, from, to) {
  tmp <- to[ match(vec, from) ]
  tmp[is.na(tmp)] <- vec[is.na(tmp)]
  return(tmp)
}

save_heatmap_pdf <- function(x, filename, width = width, height = height) {
   stopifnot(!missing(x))
   stopifnot(!missing(filename))
   pdf(filename, width=width, height=height)
   grid::grid.newpage()
   grid::grid.draw(x$gtable)
   dev.off()
}


#read in signature exposure file and make into matrix
df <- fread(sig_file)
names(df)[1] <- "sig"

#remove 10X samples
df <- df %>%
          select(sig, contains("DTB"))

#remove Ref.Sig.R3 and Ref.Sig.R8, which are 0 for all samples
df <- df[-c(3, 7), ]
matrix <- as.matrix(df, rownames = "sig")
#normalize signature exposures within a sample to sum to 1
normalized_matrix <- prop(matrix, digits = 6)
saveRDS(normalized_matrix, "normalized_matrix.rds")

#mean-center exposures across all samples within a specific signature
centered_matrix <- sweep(normalized_matrix, 1, apply(normalized_matrix, 1, mean))

#cluster samples in mean-centered matrix based on Euclidean distance metric
distance_matrix <- dist(t(centered_matrix), method = "euclidean")

cls.method = "ward.D2" #"ward.D2", "complete"
clustering <- hclust(distance_matrix, method = cls.method)


num_clusters <- 9


#create dataframe of samples and cluster assignments, where k is number of clusters.
#this helps verify that clusters used in heatmap (created from non-mean-
#centered data) are the same as those using mean-centering
sample_cluster_list <- cutree(tree = as.dendrogram(clustering), k = num_clusters)
sample_cluster_df <- as.data.frame(sample_cluster_list) %>%
                         rownames_to_column()

sample_cluster_df$cluster = cluster_swap(sample_cluster_df$sample_cluster_list, c(4,5,1,7,6,2,8,9,3), c(1:9) )
names(sample_cluster_df)[1]  = "sample"

sample_cluster_df = sample_cluster_df %>%
        select(sample,cluster)

fwrite(sample_cluster_df, paste0(output_folder, "sample_by_cluster.txt"), sep = "\t")

#get count of SVs assigned to sigs for each sample
sv_count_df <- df %>%
    pivot_longer(cols = -sig, names_to = "sample", values_to = "num_svs") %>%
    group_by(sample) %>%
    summarize(sv_count = sum(num_svs))


#grab dataframes of annotations like TDP status and chromothripsis status
tdp_df <- fread(tdp_file) %>%
              select(-category) %>%
              mutate(tdp_status = ifelse(value == "TDP", 1, 0)) %>%
              select(-value)

tmb_df <- fread(tmb_file)
              
chromothripsis_df <- fread(chromothripsis_file) %>%
                    select(-value,-category) %>%
                    mutate(chromothripsis = 1)

Chromoplexy_df <- fread(Chromoplexy_file)%>%
                    select(-value,-category) %>%
                    mutate(chromoplexy = 1)

focal_chromothripsis_df <- fread(focal_chromothripsis_file) %>%
                    mutate(focal_chromothripsis = 1)

bfb_chromothripsis_df <- fread(bfb_chromothripsis_file)  %>%
                    mutate(bfb_chromothripsis = 1)

treatment_df <- fread(treatment_file, header = TRUE) %>%
                    select(-V1, -category)
names(treatment_df)[2] <- "treatment"

therapy_df <- fread(therapy_file, header = TRUE) %>%
                    select(-V1, -category)
names(therapy_df)[2] <- "therapy"

biopsy_site_df <- fread(biopsy_site_file) %>%
                      select(-category)
names(biopsy_site_df)[2] <- "biopsy_site"

cancer_type_df <- fread(cancer_type_file) %>%
                      select(-category)
names(cancer_type_df)[2] <- "cancer_type"

erg_dna_rna_df <- fread(ets_curation_file, header = TRUE) %>%
                   filter(str_detect(category, "ERG")) %>%
                   filter(str_detect(value, "DNA\\+RNA")) %>%
                   select(-category,-value) %>%
                   mutate(erg_dna_rna = 1)

erg_all_df <- fread(ets_curation_file, header = TRUE) %>%
                   filter(str_detect(category, "ERG")) %>%
                   filter(!str_detect(value, "None")) %>%
                   select(-category,-value) %>%
                   mutate(erg_all = 1)

etv1_dna_rna_df <- fread(ets_curation_file, header = TRUE) %>%
                   filter(str_detect(category, "ETV1")) %>%
                   filter(str_detect(value, "DNA\\+RNA")) %>%
                   select(-category,-value) %>%
                   mutate(etv1_dna_rna = 1)

etv1_all_df <- fread(ets_curation_file, header = TRUE) %>%
                   filter(str_detect(category, "ETV1")) %>%
                   filter(!str_detect(value, "None")) %>%
                   select(-category,-value) %>%
                   mutate(etv1_all = 1)

etv4_dna_rna_df <- fread(ets_curation_file, header = TRUE) %>%
                   filter(str_detect(category, "ETV4")) %>%
                   filter(str_detect(value, "DNA\\+RNA")) %>%
                   select(-category,-value) %>%
                   mutate(etv4_dna_rna = 1)

etv4_all_df <- fread(ets_curation_file, header = TRUE) %>%
                   filter(str_detect(category, "ETV4")) %>%
                   filter(!str_detect(value, "None")) %>%
                   select(-category,-value) %>%
                   mutate(etv4_all = 1)

etv5_dna_rna_df <- fread(ets_curation_file, header = TRUE) %>%
                   filter(str_detect(category, "ETV5")) %>%
                   filter(str_detect(value, "DNA\\+RNA")) %>%
                   select(-category,-value) %>%
                   mutate(etv5_dna_rna = 1)

etv5_all_df <- fread(ets_curation_file, header = TRUE) %>%
                   filter(str_detect(category, "ETV5")) %>%
                   filter(!str_detect(value, "None")) %>%
                   select(-category,-value) %>%
                   mutate(etv5_all = 1)

ets_dna_rna_df <- fread(ets_curation_file, header = TRUE) %>%
                   filter(!str_detect(category, "ELK4")) %>%
                   filter(str_detect(value, "DNA\\+RNA")) %>%
                   select(-category,-value) %>%
                   mutate(ets_dna_rna = 1) %>%
                   unique()

ets_all_df <- fread(ets_curation_file, header = TRUE) %>%
                   filter(!str_detect(category, "ELK4")) %>%
                   filter(!str_detect(value, "None")) %>%
                   select(-category,-value) %>%
                   mutate(ets_all = 1) %>%
                   unique()

ar_coamp_df <- fread(ar_category_file, header = TRUE) %>%
                   filter(str_detect(category, "Enhancer/AR status")) %>%
                   filter(str_detect(value, "Coamplification")) %>%
                   select(-category,-value) %>%
                   mutate(ar_coamp = 1)

ar_selective_enhancer_df <- fread(ar_category_file, header = TRUE) %>%
                   filter(str_detect(category, "Enhancer/AR status")) %>%
                   filter(str_detect(value, "Selective Enhancer")) %>%
                   select(-category,-value) %>%
                   mutate(ar_selective_enhancer = 1)

ar_selective_ar_df <- fread(ar_category_file, header = TRUE) %>%
                   filter(str_detect(category, "Enhancer/AR status")) %>%
                   filter(str_detect(value, "Selective AR")) %>%
                   select(-category,-value) %>%
                   mutate(ar_selective_ar = 1)

ar_enhancer_td_df <- fread(ar_category_file, header = TRUE) %>%
                   filter(str_detect(category, "Enhancer TD")) %>%
                   filter(str_detect(value, "yes")) %>%
                   select(-category,-value) %>%
                   mutate(ar_enhancer_td = 1)       

ar_bfb_df <- fread(ar_category_file, header = TRUE) %>%
                   filter(str_detect(category, "BFB")) %>%
                   filter(str_detect(value, "yes")) %>%
                   select(-category,-value) %>%
                   mutate(ar_bfb = 1)   

ar_chromoplexy_df <- fread(ar_category_file, header = TRUE) %>%
                   filter(str_detect(category, "Chromoplexy")) %>%
                   filter(str_detect(value, "yes")) %>%
                   select(-category,-value) %>%
                   mutate(ar_chromoplexy = 1)   

ar_chromothripsis_df <- fread(ar_category_file, header = TRUE) %>%
                   filter(str_detect(category, "chromothripsis")) %>%
                   filter(str_detect(value, "yes")) %>%
                   select(-category,-value) %>%
                   mutate(ar_chromothripsis = 1)   


ar_ecDNA_df <- fread(ar_category_file, header = TRUE) %>%
                   filter(str_detect(category, "ecDNA")) %>%
                   filter(str_detect(value, "yes")) %>%
                   select(-category,-value) %>%
                   mutate(ar_ecDNA = 1)  


ar_intragenic_df <- fread(ar_category_file, header = TRUE) %>%
                   filter(str_detect(category, "Intragenic AR")) %>%
                   filter(str_detect(value, "yes")) %>%
                   select(-category,-value) %>%
                   mutate(ar_intragenic = 1)  

ploidy_df <- fread(ploidy_file) %>%
                 select(-V1, -category)
names(ploidy_df)[2] <- "ploidy"

purity_df <- fread(purity_file) %>%
                 select(-V1, -category)
names(purity_df)[2] <- "purity"


mutation_df <- fread(mutation_file, header = TRUE) %>%
                   select(-V1, -value) %>%
                   filter(str_detect(sample, "DTB")) %>%
                   mutate(mutation = 1) %>%
                   unique()
names(mutation_df)[2] <- "gene"
#select only mutations that are found in at least 3 samples
common_mutation_genes <- mutation_df %>%
                             count(gene) %>%
                             filter(n >= 3) %>%
                             select(gene) %>%
                             deframe()

subset_mutation_df <- mutation_df %>%
                          filter(gene %in% common_mutation_genes)

subset_mutation_df <- as.data.frame(subset_mutation_df)

#pivot mutation data wider so there is a column for each gene
long_mutation_df <- subset_mutation_df %>%
                        complete(expand(subset_mutation_df, sample, gene)) %>%
                        pivot_wider(names_from = gene, values_from = mutation)
names(long_mutation_df)[2:length(names(long_mutation_df))] <- paste0(names(long_mutation_df)[2:length(names(long_mutation_df))], "_mutation")


###CNA 
cna_gain_df <- fread(cna_file, header = TRUE) %>%
              select(-V1) %>%
              filter(str_detect(sample, "DTB")) %>%
              #deep events only
              filter(grepl("CN_Amplification", value)) %>%
              mutate(cna_gain = 1) %>%
              select(sample, category, cna_gain) %>%
              unique()

cna_del_df <- fread(cna_file, header = TRUE) %>%
              select(-V1) %>%
              filter(str_detect(sample, "DTB")) %>%
              #deep events only
              filter(grepl("CN_Homozygous_Del",value)) %>%
              mutate(cna_del = 1) %>%
              select(sample, category, cna_del) %>%
              unique()

names(cna_gain_df) <- c("sample", "gene", "cna_gain")
names(cna_del_df) <- c("sample", "gene", "cna_del")

#select only genes which have CNAs in at least 3 samples
common_cna_genes <- cna_gain_df %>%
                        count(gene) %>%
                        filter(n >= 3) %>%
                        select(gene) %>%
                        deframe()
subset_cna_gain_df  <- cna_gain_df %>%
                     filter(gene %in% common_cna_genes)

#pivot cna data wider so there is a column for each gene
long_cna_gain_df <- subset_cna_gain_df %>%
                   complete(expand(subset_cna_gain_df, sample, gene)) %>%
                   pivot_wider(names_from = gene, values_from = cna_gain)
names(long_cna_gain_df)[2:length(names(long_cna_gain_df))] <- paste0(names(long_cna_gain_df)[2:length(names(long_cna_gain_df))], "_cna_gain")


#select only genes which have CNAs in at least 3 samples
common_cna_genes <- cna_del_df %>%
                        count(gene) %>%
                        filter(n >= 3) %>%
                        select(gene) %>%
                        deframe()
subset_cna_del_df  <- cna_del_df %>%
                     filter(gene %in% common_cna_genes)

#pivot cna data wider so there is a column for each gene
long_cna_del_df <- subset_cna_del_df %>%
                   complete(expand(subset_cna_del_df, sample, gene)) %>%
                   pivot_wider(names_from = gene, values_from = cna_del)
names(long_cna_del_df)[2:length(names(long_cna_del_df))] <- paste0(names(long_cna_del_df)[2:length(names(long_cna_del_df))], "_cna_del")



#Nov.17.2021
#combine sv gene + sv flanking
sv_df <- fread(sv_file, header = TRUE) %>%
              select(-V1) %>%
              filter(str_detect(sample, "DTB")) %>%
              filter(grepl("gene",source) | grepl("svLong_flank",source) | grepl("svShort_flank",source)) %>%
              mutate(sv = 1) %>%
              select(sample,category,sv) %>%
              unique()

names(sv_df) <- c("sample", "gene", "sv")

#select only genes which have CNAs in at least 3 samples
common_sv_genes <- sv_df %>%
                        count(gene) %>%
                        filter(n >= 3) %>%
                        select(gene) %>%
                        deframe()
subset_sv_df <- sv_df %>%
                     filter(gene %in% common_sv_genes)

#pivot cna data wider so there is a column for each gene
long_sv_df <- subset_sv_df %>%
                   complete(expand(subset_sv_df, sample, gene)) %>%
                   pivot_wider(names_from = gene, values_from = sv)
names(long_sv_df)[2:length(names(long_sv_df))] <- paste0(names(long_sv_df)[2:length(names(long_sv_df))], "_sv")

#merge annotation dataframes together
annotation_df <- sample_cluster_df %>%
                     left_join(ploidy_df, by.x = "sample", by.y = "sample") %>%
                     left_join(purity_df, by.x = "sample", by.y = "sample") %>%
                     left_join(tdp_df, by.x = "sample", by.y = "sample") %>%
                     left_join(tmb_df, by.x = "sample", by.y = "sample") %>%
                     left_join(therapy_df, by.x = "sample", by.y = "sample") %>%
                     left_join(treatment_df, by.x = "sample", by.y = "sample") %>%
                     left_join(biopsy_site_df, by.x = "sample", by.y = "sample") %>%
                     left_join(cancer_type_df, by.x = "sample", by.y = "sample") %>%
                     left_join(chromothripsis_df, by.x = "sample", by.y = "sample") %>%
                     left_join(focal_chromothripsis_df, by.x = "sample", by.y = "sample") %>%
                     left_join(bfb_chromothripsis_df, by.x = "sample", by.y = "sample") %>%
                     left_join(Chromoplexy_df, by.x = "sample", by.y = "sample") %>%
                     left_join(erg_dna_rna_df, by.x = "sample", by.y = "sample") %>%
                     left_join(erg_all_df, by.x = "sample", by.y = "sample") %>%
                     left_join(etv1_dna_rna_df, by.x = "sample", by.y = "sample") %>%
                     left_join(etv1_all_df, by.x = "sample", by.y = "sample") %>%
                     left_join(etv5_dna_rna_df, by.x = "sample", by.y = "sample") %>%
                     left_join(etv5_all_df, by.x = "sample", by.y = "sample") %>%
                     left_join(etv4_dna_rna_df, by.x = "sample", by.y = "sample") %>%
                     left_join(etv4_all_df, by.x = "sample", by.y = "sample") %>%
                     left_join(ets_dna_rna_df, by.x = "sample", by.y = "sample") %>%
                     left_join(ets_all_df, by.x = "sample", by.y = "sample") %>%
                     left_join(ar_coamp_df, by.x = "sample", by.y = "sample") %>%
                     left_join(ar_selective_enhancer_df, by.x = "sample", by.y = "sample") %>%
                     left_join(ar_selective_ar_df, by.x = "sample", by.y = "sample") %>%
                     left_join(ar_enhancer_td_df, by.x = "sample", by.y = "sample") %>%
                     left_join(ar_bfb_df, by.x = "sample", by.y = "sample") %>%
                     left_join(ar_chromoplexy_df, by.x = "sample", by.y = "sample") %>%
                     left_join(ar_chromothripsis_df, by.x = "sample", by.y = "sample") %>%
                     left_join(ar_ecDNA_df, by.x = "sample", by.y = "sample") %>%
                     left_join(sv_count_df, by.x = "sample", by.y = "sample") %>%
                     left_join(ar_intragenic_df, by.x = "sample", by.y = "sample") %>%
                     left_join(long_mutation_df, by.x = "sample", by.y = "sample") %>%
                     left_join(long_sv_df, by.x = "sample", by.y = "sample") %>%
                     left_join(long_cna_gain_df, by.x = "sample", by.y = "sample") %>%
                     left_join(long_cna_del_df, by.x = "sample", by.y = "sample")

                     
annotation_df %>%
    fwrite(paste0(output_folder, "compiled_annotation_data_by_sample.txt"), sep = "\t")

#create heatmap of normalized sig exposures along with mean-centered exposure clustering
#and annotations like gene mutations and TDP status.
#you can create more heatmaps by simply selecting different columns of the annotation
#dataframe to plot. you can also set annotation_legend to TRUE, as well as control
#the color scheme of the annotations.


heatmap_output_file <- paste0(output_folder, sprintf("SVcluster_heatmap_%s.pdf",cls.method))

heatmap_annotation_df <- annotation_df %>%
                              select(sample, cluster, sv_count,
                                CCND1_cna_gain,CDK12_mutation,
                                TP53_mutation,  
                                SPOP_mutation,
                                BRCA2_mutation, 
                                chromoplexy,
                                chromothripsis,
                                TMB, ets_all
                                )

my_colour = list(
      cluster = c('1' = "#E7298A", '2' = "#66A61E",'3' = "#1B9E77", '4' = "#A6761D",
                  '5' = "#E6AB02", '6'= "#D95F02", '7'="#666666", '8'="red", '9'="#7570B3"),   
      sv_count =  brewer.pal(8, "Greens"),
      CCND1_cna_gain = c('1'='#90AECE'), CDK12_mutation = c('1'='#EEA393'), 
      TP53_mutation = c('1'='#EEA393'),SPOP_mutation = c('1'='#EEA393'),
      BRCA2_mutation = c('1'='#EEA393'), chromoplexy = c('chromoplexy'='#A3B7F9'),
      chromothripsis = c('chromothripsis'='#D35FB7'),
      TMB = brewer.pal(4, "Blues"),
      ets_all = c('ets_all'='darkgreen'))


heatmap <- pheatmap(normalized_matrix,
                     color = colorRampPalette(c("white", "orange"))(50),
                     annotation_col = column_to_rownames(heatmap_annotation_df, "sample"),
                     cutree_cols = num_clusters,
                     cluster_rows = FALSE,
                     annotation_legend = F,
                     cellheight=15, cellwidth = 10,
                     legend = F,
                     treeheight_col = 0,
                     clustering_distance_cols = "euclidean",
                     annotation_colors = my_colour,
                     clustering_method = cls.method)

ggsave(filename=heatmap_output_file, plot=heatmap, width=20, height=5, units="in") #for short version

