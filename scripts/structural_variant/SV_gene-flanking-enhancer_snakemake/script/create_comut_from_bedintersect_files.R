
args <- commandArgs(TRUE)


intersect_bed_dir <- args[1] 
intersect_files = list.files(intersect_bed_dir, pattern="intersect.*", recursive=T, full.names=T)

for (f in intersect_files){
	print(f)
	#dat = read.table(f, sep='\t', head=F)
	dat = tryCatch(read.table(f,  sep='\t', head=F), error=function(e) NULL)
	if (!is.null(dat)){
		uniq_dat = unique(dat[,c("V4","V11","V5")])
		colnames(uniq_dat) = c("sample","category","value")

		uniq_dat$category = gsub("_Up|_Dn","",uniq_dat$category)
		#print(head(uniq_dat))
		out_f  = gsub(".txt",".csv", paste0("comut_",basename(f)))
		write.csv(uniq_dat, sprintf("results/comut_files/%s",out_f), quote=F)
	}
}