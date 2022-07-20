
library(data.table)

bind_list = list()
out_dir = "results/comut_files/"

if(file.exists("./results/comut_files/comut_intersect_sv_and_gene_body.csv")){
	sv_gene = read.csv("./results/comut_files/comut_intersect_sv_and_gene_body.csv")
	sv_gene$source = "gene"
	bind_list[["gene"]] = sv_gene

}

#for short SV events
if(file.exists("./results/comut_files/comut_intersect_sv_and_flanking_up.csv")){
	sv_flankLeft = read.csv("./results/comut_files/comut_intersect_sv_and_flanking_up.csv")
	sv_flankLeft$source = "svShort_flankLeft"
	bind_list[["svShort_flankLeft"]] = sv_flankLeft
}

if(file.exists("./results/comut_files/comut_intersect_sv_and_flanking_dn.csv")){
	sv_flankRight = read.csv("./results/comut_files/comut_intersect_sv_and_flanking_dn.csv")
	sv_flankRight$source = "svShort_flankRight"
	bind_list[["sv_flankRight"]] = sv_flankRight

}
if(file.exists("./results/comut_files/comut_shortSV_enhUpMatrix_H3K27ac_single.csv")){
	sv_enhancerLeft = read.csv("./results/comut_files/comut_shortSV_enhUpMatrix_H3K27ac_single.csv")
	sv_enhancerLeft$source = "svShort_enhancerLeft"
	bind_list[["sv_enhancerLeft"]] = sv_enhancerLeft
}
if(file.exists("./results/comut_files/comut_shortSV_enhDownMatrix_H3K27ac_single.csv")){
	sv_enhancerRight = read.csv("./results/comut_files/comut_shortSV_enhDownMatrix_H3K27ac_single.csv")
	sv_enhancerRight$source = "svShort_enhancerRight"
	bind_list[["sv_enhancerRight"]] = sv_enhancerRight
}
#for long events like TRANS
if(file.exists("./results/comut_files/comut_SV_enhUpMatrix_longSV_single.csv")){
	svlong_flankLeft = read.csv("./results/comut_files/comut_SV_enhUpMatrix_longSV_single.csv")
	svlong_flankLeft$source = "svLong_flankLeft"
	bind_list[["svlong_flankLeft"]] = svlong_flankLeft
}
if(file.exists("./results/comut_files/comut_SV_enhDownMatrix_longSV_single.csv")){
	svlong_flankRight = read.csv("./results/comut_files/comut_SV_enhDownMatrix_longSV_single.csv")
	svlong_flankRight$source = "svLong_flankRight"
	bind_list[["svlong_flankRight"]] = svlong_flankRight
}
if(file.exists("./results/comut_files/comut_SV_enhUpMatrix_H3K27ac_longSV_single.csv")){
	svlong_enhancerLeft = read.csv("./results/comut_files/comut_SV_enhUpMatrix_H3K27ac_longSV_single.csv")
	svlong_enhancerLeft$source = "svLong_enhancerLeft"
	bind_list[["svlong_enhancerLeft"]] = svlong_enhancerLeft

}
if(file.exists("./results/comut_files/comut_SV_enhDownMatrix_H3K27ac_longSV_single.csv")){
	svlong_enhancerRight = read.csv("./results/comut_files/comut_SV_enhDownMatrix_H3K27ac_longSV_single.csv")
	svlong_enhancerRight$source = "svLong_enhancerRight"
	bind_list[["svlong_enhancerRight"]] = svlong_enhancerRight

}

combined_all = rbindlist(bind_list)[,-1]
table(combined_all$source)

# replace event name 
# Inversion-Bal, Inversion-FoldBack and Inversion-Unbal to "Inversion"
# Trans-Bal and Trans-Unbal to Trans
combined_all$value = gsub("Inversion-.*","Inversion",combined_all$value)
combined_all$value = gsub("Trans-.*","Trans",combined_all$value)


sample_gene_pairs = paste(combined_all$sample ,combined_all$category, sep="_")

combined_all$sample_gene_pair = sample_gene_pairs
write.csv(combined_all, paste0(out_dir, "all_comut_SV_combined_single.csv"), quote=F)

dup_events_to_check = names(sort(table(sample_gene_pairs)[table(sample_gene_pairs)!=1]))

duplicates = lapply(dup_events_to_check, FUN = function(x) combined_all[combined_all$sample_gene==x,])
bind_duplicates = rbindlist(duplicates)

write.table(bind_duplicates, paste0(out_dir, "duplicated_gene_sample_pair_SV_events.txt"), sep="\t", col.names=T, row.names=F, quote=F)
