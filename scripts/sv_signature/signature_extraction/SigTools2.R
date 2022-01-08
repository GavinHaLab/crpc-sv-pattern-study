#' runSigTools.R
#' authors: Anna Hoge, Gavin Ha
#' institution: Fred Hutchinson Cancer Research Center
#' contact: <gha@fredhutch.org>
#' date: June 7, 2020

#' requires R-3.6.0+

library(optparse)
option_list <- list(
  make_option(c("--runID"), type = "character", default = "test", help = "ID for run. Default: [%default]"),
  make_option(c("--inDir"), type = "character", help = "Input directory containing bedpe files. Required."),
  make_option(c("--ext"), type = "character", default = ".bedpe", help = "Common filename extension of bedpe files. Required."),
  make_option(c("--sampleListFile"), type = "character", default = NULL, help = "File containing list of sample IDs to consider. Required as of now."),
  make_option(c("--runExtraction"), type = "logical", default = FALSE, help = "TRUE/FALSE to run SV Extraction to predict new signatures. Default: [%default]"),
  make_option(c("--maxNumSignatures"), type = "numeric", default = 15, help = "Maximum number of signatures to consider. Will output SV extraction results for 3 to maxNumSignatures. Default: [%default]"),
  make_option(c("--signatureFile"), type = "character", default = NULL, help = "File containing published PCAWG signatures found in signature.tools.lib/data. If provided, will output results for fit to these signatures based on \"tumorType\" Optional"),
  make_option(c("--selectSignatures"), type = "character"),
  make_option(c("--refSignatureFile", type = "character", default = NULL, help = "File containing published reference signatures found in signature.tools.lib paper supplementary data.")),
  make_option(c("--minSPAN"), type = "numeric", default = 10000, help = "Minimum span (base pairs) to consider for SV events."),
  make_option(c("--excludeTools"), type = "character", default = NULL, help = "Exclude SVs based on string in \"Tool\" column"),
  make_option(c("--excludeSVsupport"), type = "character", default = NULL, help = "Exclude SV with support in \"support\" column."),
  make_option(c("--excludeSVClass"), type = "character", default = "c(\"Unknown-ShortSVwithCN\")", help = "Exclude SV class in \"CN_overlap_type\" column.."),
  make_option(c("--outDir"), type="character", default="./results", help = "Output Directory. Default: [%default]"),
  make_option(c("--libdir"), type = "character", default = NULL, help = "Script library path. "),
  make_option(c("--cores"), type="numeric", default = 1, help = "Number of cores to use for EM. Default: [%default]")
)
parseobj <- OptionParser(option_list=option_list)
opt <- parse_args(parseobj)
print(opt)
options(scipen = 0, stringsAsFactors = FALSE, bitmapType = 'cairo')

#library(signature.tools.lib)
library(data.table)
library(foreach)
library(tidyverse)

libdir <- opt$libdir
if (!is.null(libdir) && libdir != "None"){
	R_files <- list.files(libdir, pattern = ".R", full.names = TRUE)
	for (i in 1:length(R_files)){
		source(file.path(R_files[i]))
	}
	load(list.files(libdir, pattern = ".rda", full.names = TRUE)[1])
} 


numCores <- opt$cores
runID <- opt$runID
inDir <- opt$inDir
ext <- opt$ext
runExtraction <- opt$runExtraction
maxNumSignatures <- opt$maxNumSignatures
select_signatures <- eval(parse(text = opt$selectSignatures))
sampleListFile <- opt$sampleListFile
signature_file <- opt$signatureFile
ref_signature_file <- opt$refSignatureFile
outDir <- opt$outDir
dir.create(outDir)
outImage <- paste0(outDir, "/", runID, "_svSig.RData")
minSPAN <- opt$minSPAN
excludeSVClass <- eval(parse(text = opt$excludeSVClass)) #c("Unknown-ShortSVwithCN")
excludeTools <- eval(parse(text = opt$excludeTools)) #LONGRANGER perhaps?
excludeSVsupport <- eval(parse(text = opt$excludeSVsupport))  #BX perhaps?

nboots <- 20
clusteringMethod <- "MC"

save.image(outImage)

####################################
## load inputs ##

## load bedpe filenames
bedpe_files <- list.files(inDir, pattern = ext, full.names = TRUE)
ids <- gsub(ext, "", basename(bedpe_files))
names(bedpe_files) <- ids

# load sample list
samples <- fread(sampleListFile)
bedpe_files <- bedpe_files[samples$Sample]
sample_names <- names(bedpe_files)

# load signature file if provided
#select which signatures to use from file (recommended to use just those of closest
#tumor type; using many signatures not recommended)
signatures_df <- NULL
if (!is.null(signature_file)){
	signatures_df <- read.table(signature_file, sep="\t", header = T, as.is = T, check.names = FALSE)
	signatures_df <- signatures_df[, select_signatures]
	known_sigs_id <- paste0(gsub(".tsv", "", basename(signature_file)), "_selectsigs")
}

#load reference signature file
ref_signatures_df <- read.table(ref_signature_file, sep="\t", header = T, as.is = T, check.names = FALSE, row.names = 1)
ref_sigs_id <- "ref_sigs_rearr"

filterSvabaSVs <- function(x, excludeSVClass = c("Unknown-ShortSVwithCN"), columnName_SVclass = "CN_overlap_type", 
	excludeTools = NULL, columnName_tools = "Tool", excludeSVsupport = NULL, minSPAN = 10000){
	x <- copy(x)
	numSVs <- nrow(x)
	ind_SPAN <- 1:nrow(x) 
    ###fixed error messages on 11/3/2020 --Anna
	if (!is.null(minSPAN)){
        #indices of SVs that DO pass SPAN filter
		ind_SPAN <- x[SPAN == -1 | SPAN >= minSPAN, which = TRUE]
		message("Excluding ", numSVs - length(ind_SPAN), " SVs shorter than ", minSPAN)
	}
	ind_SVclass <- 1:nrow(x)
	if (!is.null(excludeSVClass)){
        #indices of SVs that DO pass SV class filter
		ind_SVclass <- x[!get(columnName_SVclass) %in% excludeSVClass, which = TRUE]
		message("Excluding ", numSVs - length(ind_SVclass), " SVs in class ", paste(excludeSVClass, collapse=", "))
	}
	ind_tools <- 1:nrow(x)
	if (!is.null(excludeTools)){
        #indices of SVs that DO pass tools filter
		ind_tools <- x[!get(columnName_tools) %in% excludeTools, which = TRUE]
		message("Excluding ", numSVs - length(ind_tools), " SVs predicted by ", paste(excludeTools, collapse=", "))
	}
	ind_support <- 1:nrow(x)
	if (!is.null(excludeSVsupport)){
        #indices of SVs that DO pass SV support filter
		ind_support <- x[grep(excludeSVsupport, support, invert = TRUE), which = TRUE]
	}
	ind <- Reduce(intersect, list(ind_SPAN, ind_SVclass, ind_tools, ind_support)) 
	message("Keeping ", length(ind), " / ", nrow(x), " events.")
	return(copy(x[ind]))
}

#####################################################################
## function to convert SVABA strand convention to BRASS convention ##
# The column "svclass" should correspond to (Sanger BRASS convention): (strand1/strand2)
#      inversion (+/-), if mates on the same chromosome
#      inversion (-/+), if mates on the same chromosome
#      deletion (+/+), if mates on the same chromosome
#      tandem-duplication (-/-), if mates on the same chromosome
#      translocation, if mates are on different chromosomes
convertSVstrandInfo <- function(x){
	# strand1 strand2 type
	# "diff chr"      "InterChr"
	# "+"     "-"     "Deletion"
	# "+"     "+"     "Inversion"
	# "-"     "-"     "Inversion"
	# "-"     "+"     "Duplication"
	x <- copy(x)
	x[, svclass := ""]
	x[type == "Deletion", svclass := "deletion"]
	x[type == "Duplication", svclass := "tandem-duplication"]
	x[type == "Inversion", svclass := "inversion"]
	x[type == "InterChr", svclass := "translocation"]
	# remove strand1 and strand2 columns since convention is different than BRASS
	x[, strand1 := NULL]; x[, strand2 := NULL]
	return(copy(x))
}

#####################################################################
## function to find closest known signature with similarity ##
#returns the list of signatures identified and the corresponding similarity. For example,
#res$RS.Breast560 = c("R1","R3","R5")
#res$cos.sim = c(0.94,0.85,0.7)
#means that Rearrangement (R) signatures 1, 3 and 5 were found, while
#the corresponding similarities to those signatures are 0.94, 0.85 and 0.7
findClosestRearrSigs_withSimilarity <- function(sigs, known_sigs, id){
  #compute cos sim matrix
  cos_sim_df <- data.frame()
  for (s in colnames(sigs)){
    for(a in colnames(known_sigs)){
      cos_sim_df[s,a] <- cos.sim(sigs[,s], known_sigs[,a])
    }
  }
  max.sim <- apply(cos_sim_df, 1, max)
  closestRS.known_sigs <- apply(cos_sim_df, 1, which.max)
  closestRS.known_sigs <- colnames(known_sigs)[closestRS.known_sigs] #paste0("R", closestRS.known_sigs)
  res <- list()
  res[["sig"]] <- names(max.sim)
  res[[id]] <- closestRS.known_sigs
  res[["cos.sim"]] <- max.sim
  return(data.frame(res, check.names = FALSE))
}


####################################
#RUN
#step 1: extract sv features for each sample and combine them all into a big df
message("Step 1: Extracting SV features from samples ...")
sv_catalogue_list <- list()
for (i in 1:length(bedpe_files)){
    sample_id <- names(bedpe_files[i])
    message(sample_id)
    bedpe_df <- fread(bedpe_files[i], check.names = FALSE, na.strings = c(".")) ###added na.strings 12/11/2020
    setnames(bedpe_df, c("Sample"), c("sample"))
    ## filter bedpe file 
	bedpe_df <- filterSvabaSVs(bedpe_df, excludeSVClass = excludeSVClass, 
		columnName_SVclass = "CN_overlap_type", excludeTools = excludeTools, 
		columnName_tools = "Tool", excludeSVsupport = excludeSVsupport, minSPAN = minSPAN)
	bedpe_df <- convertSVstrandInfo(bedpe_df)
    bedpe_catalogue <- bedpeToRearrCatalogue(as(bedpe_df, "data.frame"))
    #bedpe_catalogue has 2 parts:
    #bedpe_catalogue$rearr_catalogue -- 32 rows with counts of various sv features, like clustered_del_1-10kb or non-clustered_trans
    #bedpe_catalogue$annotated_bedpe -- this is your bedpe df as before but with is.clustered and svclass columns added if they weren't there before
    sv_catalogue_list[[sample_id]] <- bedpe_catalogue$rearr_catalogue
}
sv_catalogues <- do.call(cbind, sv_catalogue_list)

save.image(outImage)

#step 2: plot sv features for inspection
message("Step 2: Plotting rearrangement signature catalogues.")
outplot_file <- paste0(outDir, "/", runID, "_svSig_catalogues.jpg")
plotRearrSignatures(signature_data_matrix = sv_catalogues, output_file = outplot_file) #or NULL

#save sv features as text file for inspection
message("Saving rearrangement signature catalogues as text file.")
sv_catalogues_out_file <- paste0(outDir, "/", runID, "_svSig_catalogues.txt")
sv_catalogues_df <- sv_catalogues %>%
                        as.data.frame() %>%
                        rownames_to_column()
names(sv_catalogues_df)[1] <- "feature"
sv_catalogues_df %>%
    fwrite(sv_catalogues_out_file, sep = "\t")

save.image(outImage)

#step 3: fit signatures to samples
#refer online to look at the signature definitions:
#https://signal.mutationalsignatures.com/explore/cancer/17?mutationType=REARRANGEMENT
if (!is.null(signature_file)){
	message("Step 3: Fitting to known signatures with Bootstrap")
	outDir_fit <- paste0(outDir, "/signature_fit/")
	signature_fits <- SignatureFit_withBootstrap_Analysis(outdir = outDir_fit, cat = sv_catalogues, signature_data_matrix = signatures_df, type_of_mutations = "rearr", nboot = 100, nparallel = numCores)
	sv_exposures <- signature_fits$E_median_filtered
	save.image(outImage)
}

#step 4: try extracting signatures from data directly
#"According to our in silico analysis, only while using Lee–KLD in combination with filtering for the best NMF runs did we recover all ten signatures with a cosine similarity robustly above 0.9."
#"In the ‘filter’ approach, for each of 20 bootstrap catalogs, we perform 50 NMF runs and select only the runs that have optimal objective function value within a RTOL of 0.1% from the best. To limit the computation time of the clustering, we limit the number of signatures that need to be clustered by selecting at most ten random best NMF runs for each bootstrap catalog."
#The following from website https://github.com/Nik-Zainal-Group/signature.tools.lib#faq
#"In general, we advise to use 20 bootstraps, 200 repeats, the clustering with matching algorithm (CM), the KLD objective function (nmfmethod brunet) and RTOL=0.001."
if (runExtraction){
	message("Step 4: Extracting signatures.")
	outDir_extractSigs <- paste0(outDir, "/signature_extraction/")
	signature_extraction <- SignatureExtraction(cat = sv_catalogues, outFilePath = outDir_extractSigs, blacklist = c(), nrepeats = 200, nboots = nboots, clusteringMethod = clusteringMethod, useMaxMatching = TRUE, filterBestOfEachBootstrap = TRUE, filterBest_RTOL = 0.001, filterBest_nmaxtokeep = 10, nparallel = numCores, parallel = TRUE, nsig = c(3:maxNumSignatures), project = runID, type_of_extraction = "rearr", nmfmethod = "brunet", removeDuplicatesInCatalogue = FALSE, normaliseCatalogue = FALSE, plotCatalogue = TRUE, plotResultsFromAllClusteringMethods = FALSE)
	
	save.image(outImage)

	
	known_sigs_id <- paste0(gsub(".tsv", "", basename(signature_file)), "_selectsigs")	
	## for each number of signature solution...
	for (ns in 5:maxNumSignatures){			
		sigs_extracted_NS_outDir <- paste0(outDir_extractSigs, "/round_1/sig_", ns, "/")
        sigs_extracted_tsvfile <- paste0(sigs_extracted_NS_outDir, "/Sigs_plot_", runID, "_ns", ns, "_nboots", nboots, ".tsv")
		#sigs_extracted_tsvfile <- paste0(sigs_extracted_NS_outDir, "/Sigs_plot_", runID, "_ns", ns, "_nboots", nboots, "_", clusteringMethod, ".tsv")
		sigs_extracted <- read.delim(sigs_extracted_tsvfile, as.is = TRUE)
	
		#step 5: Fit current samples to extracted signatures (from )
		message("Step 5: Fitting to extracted signatures sig_", ns)
		sig_extracted_fits <- SignatureFit_withBootstrap_Analysis(outdir = sigs_extracted_NS_outDir, cat = sv_catalogues, signature_data_matrix = sigs_extracted, type_of_mutations = "rearr", nboot = 100, nparallel = numCores)

		#step 6: compute cosine similarity between extracted signatures and provided PCAWG signatures "signatures_df"
		#message("Step 6: Computing Cosine Similarity with known signatures ", known_sigs_id)
		if (!is.null(signature_file)){
			sig_cosSim_withKnownSigs <- findClosestRearrSigs_withSimilarity(sigs = sigs_extracted, known_sigs = signatures_df, id = known_sigs_id)
			outCosSim_file <- paste0(sigs_extracted_NS_outDir, "/Sigs_rearr_", known_sigs_id, "_similar_", runID, "_ns", ns, "_nboots", nboots, "_", clusteringMethod, ".tsv")
			write.table(sig_cosSim_withKnownSigs, file = outCosSim_file, row.names = FALSE, col.names = TRUE, quote = FALSE, sep = "\t")
		}
        
        #compute cosine similarity between extracted signatures and provided reference signatures
        sig_cosSim_withRefSigs <- findClosestRearrSigs_withSimilarity(sigs = sigs_extracted, known_sigs = ref_signatures_df, id = ref_sigs_id)
        outCosSim_refSigsFile <- paste0(sigs_extracted_NS_outDir, "/Sigs_rearr_", ref_sigs_id, "_similar_", runID, "_ns", ns, "_nboots", nboots, "_", clusteringMethod, ".tsv")
        write.table(sig_cosSim_withRefSigs, file = outCosSim_refSigsFile, row.names = FALSE, col.names = TRUE, quote = FALSE, sep = "\t")

	}	

	save.image(outImage)
    
    #move rdata files into a folder
    sig_extraction_rdata_folder <- paste0(outDir_extractSigs, "rdata/")
    dir.create(sig_extraction_rdata_folder)
    sig_extraction_rdata_files <- basename(Sys.glob(paste0(outDir_extractSigs, "*.Rdata")))
    for(file in sig_extraction_rdata_files){
        file.rename(from = paste0(outDir_extractSigs, file), paste0(sig_extraction_rdata_folder, file))
    }
}




