# Categorize all SV enhancer translocation + large SV (balanced/unblanced) events

# Define "flanking" SV
# Define "H3k27ac" nearby parnter breakpoints (within 1mb) but carefully consider directionalities for TRANS
# Record distance to gene (for bp1) and distance from peak (for bp2 = partner)

library(data.table)
library(dplyr)
library(GenomicRanges)
library(purrr)

args <- commandArgs(TRUE)


infile <- args[1] 
infile_H3K27ac <- args[2]
geneFile <- args[3] 
chrEnd_rds <- args[4]
H3K27ac_peaks <- args[5]

#load all translocation events
all_events = fread(infile)
chrEnds = readRDS(chrEnd_rds)
gene_dat = read.table(geneFile, sep="\t", head=T)

# #remove 'AREnhancer'
# gene_dat = gene_dat[!grepl("AREnhancer", gene_dat$symbol, fixed=TRUE),]
setnames(gene_dat,c("chromosome","start","end","symbol"),c("chr_gene","start_gene","end_gene","gene"))
gene_list = as.character(gene_dat[,"gene"])

#Restrict to genes of interest
all_events_of_interest = all_events[gene_flanking %in% gene_list,]
all_events_of_interest_with_gene_info = all_events_of_interest %>% merge(gene_dat, by.x="gene_flanking", by.y="gene")

#table(all_events_of_interest_with_gene_info$chrom == all_events_of_interest_with_gene_info$chr_gene)

#Add distance to gene
all_events_of_interest_with_gene_info[annotType=="EnhancerUp", distance_to_gene := start_gene - start]
all_events_of_interest_with_gene_info[annotType=="EnhancerDown", distance_to_gene := start - end_gene]


###############
# add H3k27ac info for partner gene 
# Define buffer region for this analysis
# if partner is fwd : set "new_start" to start_partner+1, "new_end" to start_partner + 100kb
# if partner is rev : set "new_start" to start_partner - 100kb, "new_end" to start_partner-1
###############

#partner_buffer = 100000 #100kb
partner_buffer = 1000000 #1MB

all_events_of_interest_with_gene_info[orient_partner=="fwd", partner_buffer_start := start_partner+1]
all_events_of_interest_with_gene_info[orient_partner=="fwd", partner_buffer_end := partner_buffer_start + partner_buffer]

all_events_of_interest_with_gene_info[orient_partner=="rev", partner_buffer_end := start_partner-1]
all_events_of_interest_with_gene_info[orient_partner=="rev", partner_buffer_start := partner_buffer_end - partner_buffer]

#Replace negative values to 1
all_events_of_interest_with_gene_info[partner_buffer_start < 0, partner_buffer_start := 1]

#If partner_buffer_end is greater than chr Length, replace this with max Length ("CTDP1" case)
all_events_of_interest_with_gene_info_with_chrEnd = all_events_of_interest_with_gene_info %>% merge(chrEnds, by.x="chrom_partner", by.y="Chromosome")
all_events_of_interest_with_gene_info_with_chrEnd[V1 < partner_buffer_end, partner_buffer_end := V1]

#confirm
#all_events_of_interest_with_gene_info_with_chrEnd[partner_buffer_end>V1]

## findoverlap between met peaks and partner buffer region
partner_buffer = all_events_of_interest_with_gene_info_with_chrEnd[,c("chrom_partner", "partner_buffer_start", "partner_buffer_end")]
setnames(partner_buffer, c("chr","start","end"))
partner_buffer.gr <- as(partner_buffer, "GRanges")

met_H3K27ac_peaks = read.table(H3K27ac_peaks, sep="\t")[,1:3]
colnames(met_H3K27ac_peaks) = c("chr","start","end")
met_H3K27ac_peaks.gr = makeGRangesFromDataFrame(met_H3K27ac_peaks)

peaks_hits <- as.data.table(findOverlaps(subject = partner_buffer.gr, query = met_H3K27ac_peaks.gr, type="within"))
peaks_hits_collapse <- peaks_hits %>% 
		     group_by(subjectHits) %>% 
		     mutate(met_peaks = paste0(queryHits, collapse = "_"))  %>%
		     select(subjectHits,met_peaks) %>%
		     unique()

# Annotate peaks as "partner_contain_H3K27ac" and add peak info if it exists.
all_events_of_interest_with_gene_info_with_chrEnd[peaks_hits_collapse$subjectHits, partner_contain_H3K27ac := "yes"]
all_events_of_interest_with_gene_info_with_chrEnd[peaks_hits_collapse$subjectHits, partner_H3K27ac_info := peaks_hits_collapse$met_peaks]

# Add distance info to nearest peak from partner_breakpoint
# If orient_partner is "rev", pick the "end" of the last peak to calculate the distance
# If orient_partner is "fwd", pick the "start" of the first peak to calculate the distance

partner_rev_idx =  which(all_events_of_interest_with_gene_info_with_chrEnd$orient_partner == "rev")
partner_fwd_idx =  which(all_events_of_interest_with_gene_info_with_chrEnd$orient_partner == "fwd")

last_peak_for_rev =  unlist(map(strsplit(all_events_of_interest_with_gene_info_with_chrEnd[partner_rev_idx,"partner_H3K27ac_info"][[1]],"_"),last))
first_peak_for_fwd =  unlist(lapply(strsplit(all_events_of_interest_with_gene_info_with_chrEnd[partner_fwd_idx,"partner_H3K27ac_info"][[1]],"_"),'[[',1))

all_events_of_interest_with_gene_info_with_chrEnd[partner_rev_idx, nearest_peak := last_peak_for_rev]
all_events_of_interest_with_gene_info_with_chrEnd[partner_fwd_idx, nearest_peak := first_peak_for_fwd]

nearest_peak_pos_list <- lapply(1:nrow(all_events_of_interest_with_gene_info_with_chrEnd), FUN = function(x) 
						if(all_events_of_interest_with_gene_info_with_chrEnd[x,"orient_partner"]=="rev"){
							met_H3K27ac_peaks$end[as.integer(all_events_of_interest_with_gene_info_with_chrEnd[x, "nearest_peak"][[1]])]
						}else{
							met_H3K27ac_peaks$start[as.integer(all_events_of_interest_with_gene_info_with_chrEnd[x, "nearest_peak"][[1]])]
						}
					)
# Record nearest peak and distance
all_events_of_interest_with_gene_info_with_chrEnd[ ,nearest_peak_pos := unlist(nearest_peak_pos_list)]
all_events_of_interest_with_gene_info_with_chrEnd[ ,distance_to_nearest_peak := abs(nearest_peak_pos-start_partner)]

###################
# Add flanking information (definition : SV arc should go toward gene body)
# If annotType == "EnhancerUp" (geneLeft) -> orient == "fwd"
# If annotType == "EnhancerDown" (geneRight) -> orient == "rev"
###################
all_events_of_interest_with_gene_info_with_chrEnd[annotType=="EnhancerUp" & orient == "fwd", flanking := "flankUp"]
all_events_of_interest_with_gene_info_with_chrEnd[annotType=="EnhancerDown" & orient == "rev", flanking := "flankDown"]
all_events_of_interest_with_gene_info_with_chrEnd[,ID := paste0(Sample,"_",chrom,"_",start,"_",gene_flanking,"_",svID)] 


###################
# Add Enhancer information from bp1
# need it because we might want genes without enhancer nearby and gained enhancer through partners
# All directions can be possible as long as SV arcs within enhancer region contains H3K27ac peaks.
# Separate into case1 and case2 further (towardGene and awayGene)
###################

h3k27ac_only = fread(infile_H3K27ac)

#selected translocation events only
h3k27ac_events = h3k27ac_only[, c("Sample","chr","start","end","orient","EnhancerUp.single","EnhancerDown.single","svID")]

h3k27ac_events[orient == "fwd", pos := start]
h3k27ac_events[orient == "rev", pos := end]
h3k27ac_events[EnhancerUp.single != "", gene := EnhancerUp.single]
h3k27ac_events[EnhancerDown.single != "", gene := EnhancerDown.single]

h3k27ac_events[EnhancerUp.single != "" & orient == "fwd", enhancer := "towardGene_enhancerUp"]
h3k27ac_events[EnhancerUp.single != "" & orient == "rev", enhancer := "awayGene_enhancerUp"]

h3k27ac_events[EnhancerDown.single != "" & orient == "fwd", enhancer := "awayGene_EnhancerDown"]
h3k27ac_events[EnhancerDown.single != "" & orient == "rev", enhancer := "towardGene_EnhancerDown"]

h3k27ac_events[,ID := paste0(Sample,"_",chr,"_",pos,"_",gene,"_",svID)]

###################
# Merge dataframe all_events_of_interest_with_gene_info with h3k27ac_events
###################

merge_events = all_events_of_interest_with_gene_info_with_chrEnd %>% merge(h3k27ac_events, by.x="ID", by.y="ID", all.x=T, all.y=T)
write.table(merge_events, "./results/longSV_analysis/all_longSV_flanking_events_with_enhancer_single_summary_full_version.txt", sep="\t", quote=F, row.names=F)

#Reorder and then output into new file 

arranged_dat = merge_events[,c("Sample.x","chrom","start.x","orient.x","gene_flanking","SPAN","svID.x","annotType","distance_to_gene","CN_overlap_type","flanking","enhancer","chrom_partner","start_partner","orient_partner","gene_partner","partner_buffer_start","partner_buffer_end","partner_contain_H3K27ac","partner_H3K27ac_info","nearest_peak","distance_to_nearest_peak")]
setnames(arranged_dat, c("Sample.x","start.x","orient.x","gene_flanking","enhancer","svID.x","annotType","chrom_partner","start_partner","orient_partner","gene_partner"), c("Sample","start","orient","flanking_gene_1mb","enhancer_nearby","svID","up_dn","partner_chrom","partner_start","partner_orient","partner_gene_1mb"))

arranged_dat[,up_dn := gsub("Enhancer","",up_dn)]
arranged_dat[,distance_sum := distance_to_nearest_peak + distance_to_gene]

write.table(arranged_dat, "./results/longSV_analysis/all_longSV_flanking_events_with_enhancer_single_summary_short_version.txt", sep="\t", quote=F, row.names=F)


###################
# Create comut files by filtering information
# Create this for each up/dn
###################

# 1) flanking only 
flanking_all = arranged_dat[!is.na(flanking),]
flanking_up_uniq = unique(flanking_all[up_dn=="Up",c("Sample","flanking_gene_1mb","CN_overlap_type")])
flanking_dn_uniq = unique(flanking_all[up_dn=="Down",c("Sample","flanking_gene_1mb","CN_overlap_type")])
setnames(flanking_up_uniq, c("Sample","flanking_gene_1mb","CN_overlap_type"),c("sample","category","value"))
setnames(flanking_dn_uniq, c("Sample","flanking_gene_1mb","CN_overlap_type"),c("sample","category","value"))

if(nrow(flanking_up_uniq)!=0) write.csv(flanking_up_uniq, "./results/comut_files/comut_SV_enhUpMatrix_longSV_single.csv", quote=F)
if(nrow(flanking_dn_uniq)!=0) write.csv(flanking_dn_uniq, "./results/comut_files/comut_SV_enhDownMatrix_longSV_single.csv", quote=F)

# 2) flanking + H3K27ac in partner
# translocation events where bp1 in in flanking region (without enhancer mark) and partner has enhancer mark near by 1mb buffer region 
enhancer_all = arranged_dat[!is.na(flanking) & is.na(enhancer_nearby) & partner_contain_H3K27ac=="yes" ,]
enhancer_up_uniq = unique(enhancer_all[up_dn=="Up",c("Sample","flanking_gene_1mb","CN_overlap_type")])
enhancer_dn_uniq = unique(enhancer_all[up_dn=="Down",c("Sample","flanking_gene_1mb","CN_overlap_type")])
setnames(enhancer_up_uniq, c("Sample","flanking_gene_1mb","CN_overlap_type"),c("sample","category","value"))
setnames(enhancer_dn_uniq, c("Sample","flanking_gene_1mb","CN_overlap_type"),c("sample","category","value"))

if(nrow(enhancer_up_uniq)!=0) write.csv(enhancer_up_uniq, "./results/comut_files/comut_SV_enhUpMatrix_H3K27ac_longSV_single.csv", quote=F)
if(nrow(enhancer_dn_uniq)!=0) write.csv(enhancer_dn_uniq, "./results/comut_files/comut_SV_enhDownMatrix_H3K27ac_longSV_single.csv", quote=F)






