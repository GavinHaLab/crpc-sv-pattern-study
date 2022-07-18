#Parse SV events for sample-gene pair produced from script/makeSVmatrixByGenes_with_H3K27ac_enhancer_for_intraSV.R

library(dplyr)
library(data.table)
library(tidyr)
options(stringsAsFactors = FALSE, scipen=999)

args <- commandArgs(TRUE)


prefix <- args[1] 
genefile <- args[2]
out_dir <- args[3]

files = list.files("results/shortSV_enh_analysis", pattern=".*Matrix_H3K27ac.*txt",full.names=TRUE)


#load gene list
gene_dat = read.table(genefile, sep="\t", head=T)
gene_list = gene_dat$symbol

for(f in files){
	print(f)

	d = read.table(f, sep='\t', head=T, quote="", fill=FALSE, check.names=FALSE)
	match_id = as.vector(match(gene_list, colnames(d)))
	match_id = match_id[!is.na(match_id)]
	names(match_id) = gene_list
	sv_comut_list = list()

	samples = rownames(d)

	for(i in 1:length(match_id)){
		id = match_id[i]
		#print(id)
		gene = names(id)
		print(gene)
		sv = d[,id]
		idx_avail = which(sv !="NaN")

		if( length(idx_avail)!=0 )	{

			sv_avail <- tibble(
			  sample = samples[idx_avail],
			  svNum = sv[idx_avail]
			)

			sv_avail_separate <- as.data.table(separate_rows(sv_avail, svNum, convert = TRUE,sep="_"))
			
			#assign events
			sv_avail_separate[svNum==1,value:="TandemDup"]
			sv_avail_separate[svNum==2,value:="Deletion"]
			sv_avail_separate[svNum==3,value:="Inversion-Unbal"]
			sv_avail_separate[svNum==4,value:="Inversion-Bal"]
			sv_avail_separate[svNum==5,value:="Inversion-FoldBack"]
			sv_avail_separate[svNum==6,value:="Trans-Unbal"]
			sv_avail_separate[svNum==7,value:="Trans-Bal"]
			sv_avail_separate[svNum==8,value:="Unbalanced"]
			sv_avail_separate[svNum==9,value:="Balanced"]
			sv_avail_separate[svNum==10,value:="LR-shortDup"]

			complete_sv_avail_separate = unique(sv_avail_separate[complete.cases(sv_avail_separate)])
			complete_sv_avail_separate$category = gene

			sv_comut_list[[gene]] = unique(complete_sv_avail_separate[,c("sample","category","value")])
		}

	}

	if( length(sv_comut_list)!=0 ){
		out_comut = do.call(rbind, sv_comut_list)

		out_f = gsub(".txt",".csv", gsub("results/shortSV_enh_analysis/Lucap_cohort_","comut_",f))

		write.csv(out_comut, paste0(out_dir, out_f), quote=F)
	}
}
write.csv("Done!", "results/shortSV_enh_analysis/shortSV_enh_analysis_complete.txt", quote=F)

