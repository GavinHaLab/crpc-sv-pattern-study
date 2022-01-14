##This will draw figures 2-B

library(plyr)

options(scipen=999)
options(stringsAsFactors = FALSE)


# pvalue threshold for fishers test
thres <- 0.05 

# ##################################################
# # meta data
# ##################################################

gene_file <- "../../../metadata/list.gene.fishhook"

gene_list_prostate <- read.table(gene_file, header=F,
                                 as.is=T)[, 1]


sample_list_mcrpc <- read.table("../../../metadata/list.sample.mcrpc", header=F,
                                as.is=T)[, 1]
sample_list_primary <- read.table("../../../metadata/list.sample.primary", header=F,
                                as.is=T)[, 1]

n_primary <- length(sample_list_primary)
n_mcrpc <- length(sample_list_mcrpc)


##################################################
# load calls
##################################################
d_cn_primary <- read.table("../../../metadata/comut_primary_cn.tsv", header=F, as.is=T)
d_sv_primary <- read.table("../../../metadata/comut_primary_sv.tsv", header=F, as.is=T)
d_snv_primary <- read.table("../../../metadata/comut_primary_snv.tsv", header=F, as.is=T)

d_cn_mcrpc <- read.table("../../../metadata/comut_mcrpc_cn.csv", sep="\t",header=T, as.is=T)
d_sv_mcrpc <- read.table("../../../metadata/comut_mcrpc_sv.csv", sep="\t",header=T, as.is=T)
d_snv_mcrpc <- read.table("../../../metadata/comut_mcrpc_snv.csv", sep="\t",header=T, as.is=T)

d_cn_mcrpc <- subset(d_cn_mcrpc, sample %in% sample_list_mcrpc)
d_sv_mcrpc <- subset(d_sv_mcrpc, sample %in% sample_list_mcrpc)
d_snv_mcrpc <- subset(d_snv_mcrpc, sample %in% sample_list_mcrpc)

d_cn_mcrpc <- subset(d_cn_mcrpc, category %in% gene_list_prostate)
d_sv_mcrpc <- subset(d_sv_mcrpc, category %in% gene_list_prostate)
d_snv_mcrpc <- subset(d_snv_mcrpc, category %in% gene_list_prostate)

d_cn_primary <- subset(d_cn_primary, V2 %in% gene_list_prostate)
d_sv_primary <- subset(d_sv_primary, V2 %in% gene_list_prostate)
d_snv_primary <- subset(d_snv_primary, V2 %in% gene_list_prostate)


##################################################
# generate status flags
##################################################
status_template_mcrpc <- matrix(data=F, nrow=length(gene_list_prostate),
                                ncol=length(sample_list_mcrpc),
                                dimnames=list(gene_list_prostate,
                                              sample_list_mcrpc))
status_template_primary <- matrix(data=F, nrow=length(gene_list_prostate),
                                  ncol=length(sample_list_primary),
                                  dimnames=list(gene_list_prostate,
                                                sample_list_primary))
has_cn_mcrpc <- status_template_mcrpc
has_sv_mcrpc <- status_template_mcrpc
has_snv_mcrpc <- status_template_mcrpc
has_cn_primary <- status_template_primary
has_sv_primary <- status_template_primary
has_snv_primary <- status_template_primary

for(i in 1:nrow(d_cn_mcrpc)) { has_cn_mcrpc[d_cn_mcrpc[i, 2], d_cn_mcrpc[i, 1]] <- T }
for(i in 1:nrow(d_sv_mcrpc)) { has_sv_mcrpc[d_sv_mcrpc[i, 2], d_sv_mcrpc[i, 1]] <- T }
for(i in 1:nrow(d_snv_mcrpc)) { has_snv_mcrpc[d_snv_mcrpc[i, 2], d_snv_mcrpc[i, 1]] <- T }
for(i in 1:nrow(d_cn_primary)) { has_cn_primary[d_cn_primary[i, 2], d_cn_primary[i, 1]] <- T }
for(i in 1:nrow(d_sv_primary)) { has_sv_primary[d_sv_primary[i, 2], d_sv_primary[i, 1]] <- T }
for(i in 1:nrow(d_snv_primary)) { has_snv_primary[d_snv_primary[i, 2], d_snv_primary[i, 1]] <- T }

count_data <- data.frame(
  gene=gene_list_prostate,
  altered_mcrpc=rowSums(has_cn_mcrpc | has_sv_mcrpc | has_snv_mcrpc, na.rm=T),
  cn_mcrpc=rowSums(has_cn_mcrpc, na.rm=T),
  sv_mcrpc=rowSums(has_sv_mcrpc, na.rm=T),
  snv_mcrpc=rowSums(has_snv_mcrpc, na.rm=T),
  altered_primary=rowSums(has_cn_primary | has_sv_primary | has_snv_primary, na.rm=T),
  cn_primary=rowSums(has_cn_primary, na.rm=T),
  sv_primary=rowSums(has_sv_primary, na.rm=T),
  snv_primary=rowSums(has_snv_primary, na.rm=T)
  )


####################
# fisher's exact test for altered genes between mcrpc and primary
####################

for(type in c("altered","sv")){

  ## 1) find genes enriched more in mcrpc compared to primary

  mcrpc_vs_primaryMat <- lapply(1:nrow(count_data), function(x){
    c(count_data[x, sprintf("%s_mcrpc",type)], n_mcrpc-count_data[x, sprintf("%s_mcrpc",type)], 
      count_data[x, sprintf("%s_primary",type)], n_primary-count_data[x, sprintf("%s_primary",type)])
  })
  mcrpc_vs_primaryMat <- matrix(unlist(mcrpc_vs_primaryMat), ncol=4, byrow=T)
  colnames(mcrpc_vs_primaryMat) <- c("mcrpc", "primary", type, sprintf("Not_%s",type))
  rownames(mcrpc_vs_primaryMat) <- gene_list_prostate
  
  #Pearson's Chi-squared test
  fisher.pval <- aaply(.data=mcrpc_vs_primaryMat, .margins=1, .fun=function(x){
    xtab = t(matrix(x, nrow=2, ncol=2, byrow=F, dimnames=list(names(x)[3:4], names(x)[1:2])))
    fisher.test(xtab, alternative="greater")$p.value
  }, .progress=progress_text(style=1))

  
  sig_genes <- sort(fisher.pval)[sort(fisher.pval) <= thres]
  count_data[match(names(sig_genes),count_data$gene),sprintf("sig_mcrpc_%s",type)] <-  sig_genes

  fisher.pval_out <- as.data.frame(sort(fisher.pval))
  write.table(fisher.pval_out, sprintf("mcrpc_enriched_%s_genes_fisher.test.pval.txt",type), sep="\t", col.names=F, quote=F)

## 2) find genes enriched more in primary compared to mcrpc

  primary_vs_mcrpcMat <- lapply(1:nrow(count_data), function(x){
    c(count_data[x, sprintf("%s_primary",type)], n_primary-count_data[x, sprintf("%s_primary",type)], 
      count_data[x, sprintf("%s_mcrpc",type)], n_mcrpc-count_data[x, sprintf("%s_mcrpc",type)])
  })
  primary_vs_mcrpcMat <- matrix(unlist(primary_vs_mcrpcMat), ncol=4, byrow=T)
  colnames(primary_vs_mcrpcMat) = c("primary", "mcrpc", type, sprintf("Not_%s",type))
  rownames(primary_vs_mcrpcMat) = gene_list_prostate

  fisher.pval <- aaply(.data=primary_vs_mcrpcMat, .margins=1, .fun=function(x){
    xtab = t(matrix(x, nrow=2, ncol=2, byrow=F, dimnames=list(names(x)[3:4], names(x)[1:2])))
    fisher.test(xtab, alternative="greater")$p.value
  }, .progress=progress_text(style=1))

  fisher.qval <- p.adjust(fisher.pval, method = "fdr")
  
  sig_genes <- sort(fisher.pval)[sort(fisher.pval) <= thres]
  count_data[match(names(sig_genes), count_data$gene), sprintf("sig_primary_%s",type)] <- sig_genes

  fisher.pval_out <- as.data.frame(sort(fisher.pval))
  write.table(fisher.pval_out, sprintf("primary_enriched_%s_genes_fisher.test.pval.txt",type), sep="\t", col.names=F, quote=F)

}

source("draw.R")
