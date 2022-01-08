#################################
#run dndscv to cobined mCRPC cohort
#################################
options(scipen=999, stringsAsFactors = FALSE)

library(dndscv)
require(liftOver)
library(rtracklayer)
library(GenomicRanges)
library(data.table)

#load data
crpc10x_snv = fread("../../../metadata/crpc10X.lolopicker.stats_calls.annovar.hg38_multianno.cdk12_rescue.txt")
crpc10x_indel = fread("../../../metadata/crpc10X.strelka.passed.somatic.indels.annovar.hg38_multianno.cdk12_rescue.txt")
wcdt_snv_indel = fread("../../../metadata/wcdt.somatic_pass.annovar.hg38_multianno.txt")

crpc10x_snv_subset = crpc10x_snv[,c("sample","Chr","Start","Ref","Alt")]
crpc10x_indel_subset = crpc10x_indel[,c("sample","Chr","Start","Ref","Alt")]
wcdt_snv_indel_subset = wcdt_snv_indel[,c("sample","Chr","Start","Ref","Alt")]

combined_mutations = rbind(crpc10x_snv_subset, crpc10x_indel_subset, wcdt_snv_indel_subset)


#################################
#Driver discovery (positive selection) in cancer exomes/genomes
#no covariates available for hg38 yet
#################################
#To reduce the risk of false positives and increase the signal to noise ratio, this example will only consider mutations in Cancer Gene Census genes (v81).
data("cancergenes_cgc81", package="dndscv")

#Use unique indel sites instead of the total number of indels (it tends to be more robust)
genes2remove = c("C2orf44", "CASC5", "CDKN2A.p14arf", "CDKN2A.p16INK4a", "FAM46C", 
	"KIAA1598", "LHFP", "MDS2", "MKL1", "MLLT4", "WHSC1", "WHSC1L1")


dndsout = dndscv(combined_mutations, refdb="../../../metadata/RefCDS_human_GRCh38.p12.rda", cv=NULL, gene_list=setdiff(known_cancergenes, genes2remove))

sel_cv = dndsout$sel_cv
print(head(sel_cv), digits = 3)
write.table(sel_cv,"dndscv_results.txt", sep="\t", row.names=F, col.names=T, quote=F)

signif_q_genes = sel_cv[sel_cv$qglobal_cv < 0.1,]
signif_p_genes = sel_cv[sel_cv$pglobal_cv < 0.05,]

write.table(signif_q_genes,"signif_genes_dndscv_q.1.txt", sep="\t",row.names=F,col.names=T,quote=F)
write.table(signif_p_genes,"signif_genes_dndscv_p.05.txt", sep="\t",row.names=F,col.names=T,quote=F)

print(dndsout$globaldnds)
#     name       mle     cilow   cihigh
#wmis wmis 0.9903614 0.9225230 1.063188
#wnon wnon 1.1271962 0.9583578 1.325780
#wspl wspl 1.0198626 0.8242946 1.261830
#wtru wtru 1.0870550 0.9500630 1.243800
#wall wall 0.9994117 0.9318956 1.071820

head(dndsout$annotmuts)



