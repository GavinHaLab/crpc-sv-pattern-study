#create bed for each breakpoints and gene interval files

library(data.table)
library(dplyr)

args <- commandArgs(TRUE)

bedpe_file <- args[1] 
sample_file <- args[2]
gene_file <- args[3] 
out_dir <- "./results/prepare_bed_and_intervals/"

sv <- fread(bedpe_file) 
sample_list <- fread(sample_file)[[1]] 
genes <- fread(gene_file) 
genes <- genes[chromosome %in% paste0("chr",c(1:22,"X")),]

#filter SV by sample_list
sv <- sv[Sample %in% sample_list,]


#create svID to track down
sv[,svID := 1:nrow(sv)]
write.table(sv, paste0(out_dir, gsub(".bedpe","_with_svID.bedpe",bedpe_file)), sep="\t", row.names=F, col.names=T, quote=F)

sv_bp1 <- sv[,c("chrom1","start1","end1","Sample","CN_overlap_type","SPAN","svID")]
sv_bp2 <- sv[,c("chrom2","start2","end2","Sample","CN_overlap_type","SPAN","svID")]
setnames(sv_bp1,c("chrom1","start1","end1"),c("chrom","start","end"))
setnames(sv_bp2,c("chrom2","start2","end2"),c("chrom","start","end"))

rbind_sv = rbind(sv_bp1, sv_bp2)
write.table(rbind_sv, paste0(out_dir,"all_sv_events_breakpoints_only.txt"), sep="\t", row.names=F, col.names=F, quote=F)

#create intervals for gene body 
#write.table(genes[,c("chromosome","start","end","symbol")], "crpc_genes.txt",sep="\t",row.names=F,col.names=F,quote=F)
write.table(genes[, c("chromosome","start","end","symbol")], paste0(out_dir,"gene_body_intervals.txt"), sep="\t", row.names=F, col.names=F, quote=F)

#add 1mb buffer for up (gene left)/dn (gene right) stream
buffer <- 1e6

#for "EnhancerUp" (geneLeft)
genes[, EnhancerUp_end := start-1]
genes[, EnhancerUp_start :=  as.integer(EnhancerUp_end - buffer)]
genes[EnhancerUp_start < 0, EnhancerUp_start := 1]

#for "EnhancerDown" (geneRight)
genes[, EnhancerDown_start := end+1]
genes[, EnhancerDown_end :=  as.integer(EnhancerDown_start + buffer)]
#genes[EnhancerDown_end < 0, EnhancerDown_end := 1]


#replace End position when it's greater than chromosomal length
chrEnds <- readRDS("chrEnds.rds")

genes_EnhancerUp <- genes[,c("chromosome","EnhancerUp_start","EnhancerUp_end","symbol")]
genes_EnhancerUp <- genes_EnhancerUp %>% 
					merge(chrEnds, by.x="chromosome", by.y="Chromosome")

genes_EnhancerUp[V1 < EnhancerUp_end, EnhancerUp_end := V1]
genes_EnhancerUp[, symbol := paste0(symbol,"_Up")]
write.table(genes_EnhancerUp[,c("chromosome","EnhancerUp_start","EnhancerUp_end","symbol")], paste0(out_dir,"genes_upstream_intervals.txt"), sep="\t", row.names=F, col.names=F, quote=F)

genes_EnhancerDown <- genes[,c("chromosome","EnhancerDown_start","EnhancerDown_end","symbol")]
genes_EnhancerDown <- genes_EnhancerDown %>% 
					merge(chrEnds, by.x="chromosome", by.y="Chromosome")

genes_EnhancerDown[V1 < EnhancerDown_end, EnhancerDown_end := V1]
genes_EnhancerDown[, symbol := paste0(symbol,"_Dn")]

write.table(genes_EnhancerDown[,c("chromosome","EnhancerDown_start","EnhancerDown_end","symbol")], paste0(out_dir,"genes_dnstream_intervals.txt"), sep="\t", row.names=F, col.names=F, quote=F)


