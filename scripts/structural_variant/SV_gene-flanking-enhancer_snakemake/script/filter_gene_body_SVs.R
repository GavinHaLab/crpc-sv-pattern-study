library(data.table)

args <- commandArgs(TRUE)

bedpe_file <- args[1] 
gene_body_sv_file <- args[2] 
out_dir <- "./results/intersect_files/"

#create bed for each breakpoints 
original_sv <- fread(bedpe_file) #input sv bedpe
gene_body_sv <- fread(gene_body_sv_file)

svID_for_gene_body <- unique(gene_body_sv$V7)

message("unique gene body events: ", length(svID_for_gene_body))

#remove rows with sv ID (mapped to genes)
remaining_sv <- original_sv[-svID_for_gene_body,]

message("Remaining events: ", nrow(remaining_sv))

#divide into two groups - short events and long events becuase we need to treat them separately!

#1) SPAN greater than 10MB and Trans
maxSPAN <- 10000000 #10MB

long_events_id <- remaining_sv[SPAN >= maxSPAN | SPAN ==-1 , which = TRUE]
short_sv <- remaining_sv[-long_events_id, ] 
long_sv <- remaining_sv[long_events_id, ] 

message("short_sv events: ", nrow(short_sv))
message("long_sv events: ", nrow(long_sv))

#this file is going to be used for creating mapping table (and keep records for partner genes)
write.table(long_sv, paste0(out_dir,"long_sv_minus_gene_body_transection_full.txt"), sep="\t", row.names=F, col.names=T, quote=F) 

#this file is going to be used for Enhancer analysis
write.table(short_sv, paste0(out_dir,"short_sv_minus_gene_body_transection_full.txt"), sep="\t", row.names=F, col.names=T, quote=F)

#re-construct bed for short sv events for bedtools
short_sv_bp1 <- short_sv[,c("chrom1","start1","end1","Sample","CN_overlap_type","SPAN","svID")]
short_sv_bp2 <- short_sv[,c("chrom2","start2","end2","Sample","CN_overlap_type","SPAN","svID")]
setnames(short_sv_bp1,c("chrom1","start1","end1"),c("chrom","start","end"))
setnames(short_sv_bp2,c("chrom2","start2","end2"),c("chrom","start","end"))

rbind_short_sv <- rbind(short_sv_bp1, short_sv_bp2)
write.table(rbind_short_sv, paste0(out_dir,"short_sv_minus_gene_body_transection.txt"),sep="\t",row.names=F,col.names=F,quote=F)
