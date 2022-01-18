library(argparse)
parser <- ArgumentParser(description="Run fishHook using given parameters.")
parser$add_argument("-b", "--bin", type="integer", default=100,
  help="Bin size in kb. Default: %(default)d", metavar="10/100/1000")
parser$add_argument("-c", "--use-covariate", action="store_true", dest="use_cov",
  default=F, help="Turn on to use all covariates")
parser$add_argument("--no-centromere", action="store_true", dest="use_cen",
  default=F, help="Turn on to filter out centromeric SVs")
parser$add_argument("-e", "--frac-eligible", type="double",
  dest="frac_eligible", default=0.75,
  help="Minimum fraction of eligible region per bin. Default: %(default)f")
parser$add_argument("-i", "--input", default="none", dest="file_sv",
  help="Input file for SV. Format must be bedpe.")
parser$add_argument("-f", "--sv-filter", required=T, default="none", dest="sv_filter",
  help="Apply type filter for SVs. Available filters: unknown, 1kb, 10kb, duplication, invdel, translocation")
parser$add_argument("--filter-cmd", required=F, default="none", dest="filter_cmd",
  help="Custom command for filtering SVs. Must use \"-f custom\" at the same time.")
parser$add_argument("--out-list", required=T, dest="outlist",
  help="File name for list of loci")
parser$add_argument("--out-plot", required=T, dest="outplot",
  help="File name for qqplot")
args <- parser$parse_args()
if (!args$sv_filter %in% c("none", "unknown", "1kb", "10kb", "duplication",
                           "invdel", "translocation", "custom")) {
  message(sprintf("Invalid parameter: --sv-filter %s", args$sv_filter))
  quit()
}

########## covariate file names ########## 
# 1. nucleotide composition
file_nt_comp <- sprintf("covariate/nt_composition.hg38.rds")

# 2. replication timing
file_rep_timing <- "covariate/reptiming_LNCaP.hg38.smoothed.bed"

# 3. DHS
file_dhs <- "covariate/dhs_LNCaP_rep1.hg38.bedGraph"

# 4. retrotransposons
file_rmsk_line <- "covariate/retrotransposon.line.hg38.bed"
file_rmsk_sine <- "covariate/retrotransposon.sine.hg38.bed"
file_rmsk_ltr <- "covariate/retrotransposon.ltr.hg38.bed"
file_rmsk_dna <- "covariate/retrotransposon.dna.hg38.bed"
file_rmsk_simple <- "covariate/retrotransposon.simple.hg38.bed"

# 5. chromatin states
file_chrom_hmm <- "covariate/chromHMM_LNCaP.hg38.bed"

# 6. common fragile sites
file_cfs <- "covariate/common_fragile_sites.hg38.txt"

############ other input files ###########
file_eligible <- "covariate/um75-hs38DH.complement.bed"
file_sv <- args$file_sv
file_cytoband <- "covariate/hg38.cytoband.bed"
########################################## 

suppressMessages(library(fishHook))
suppressMessages(library(rtracklayer))
suppressMessages(library(BSgenome.Hsapiens.NCBI.GRCh38))
message("Loading genome")
genome_hg38 <- BSgenome.Hsapiens.NCBI.GRCh38
genome_hg38_info <- seqinfo(genome_hg38)[seqlevels(genome_hg38)[1:23]]
chr_window <- tileGenome(genome_hg38_info, tilewidth=args$bin*1000,
                         cut.last.tile.in.chrom=T)
eligible_reg <- import.bed(file_eligible, seqinfo=genome_hg38_info)

# SV info
# load SV now
message("Loading SV")
sv_data_raw <- read.table(file_sv, header=T, as.is=T, sep="\t")
sv_data_raw$UID <- seq(1, nrow(sv_data_raw))
sv_data <- sv_data_raw
# additional filtering based on SV types
if (args$sv_filter == "unknown") {
  sv_data <- subset(sv_data_raw, CN_overlap_type != "Unknown-ShortSVwithCN")
} else if (args$sv_filter == "1kb") {
  sv_data <- subset(sv_data_raw, CN_overlap_type != "Unknown-ShortSVwithCN" &
                    (SPAN == -1 | SPAN >= 1000))
} else if (args$sv_filter == "10kb") {
  sv_data <- subset(sv_data_raw, CN_overlap_type != "Unknown-ShortSVwithCN" &
                    (SPAN == -1 | SPAN >= 10000))
} else if (args$sv_filter == "duplication") {
  sv_data <- subset(sv_data_raw, type == "Duplication" & SPAN >= 1000)
} else if (args$sv_filter == "invdel") {
  sv_data <- subset(sv_data_raw, (type == "Deletion" | type == "Inversion")
                                & SPAN >= 1000)
} else if (args$sv_filter == "translocation") {
  sv_data <- subset(sv_data_raw, SPAN == -1)
} else if (args$sv_filter == "custom") {
  message(sprintf("Filtering SV based on custom conditions [n=%d]", nrow(sv_data)))
  eval(parse(text=args$filter_cmd))
}
bp1_data <- data.frame(start=sv_data$start1,
                       end=sv_data$start1+1,
                       strand=rep('*', nrow(sv_data)),
                       seqnames=gsub("chr", "", sv_data$chrom1))
bp2_data <- data.frame(start=sv_data$start2,
                       end=sv_data$start2+1,
                       strand=rep('*', nrow(sv_data)),
                       seqnames=gsub("chr", "", sv_data$chrom2))
sv_gr <- dt2gr(rbind(bp1_data, bp2_data))
sv_gr <- sortSeqlevels(sv_gr)
seqinfo(sv_gr) <- genome_hg38_info
sv_gr <- trim(sv_gr)
# (optional) filtering out centromeric regions
if (args$use_cen) {
  cytoband <- read.table(file_cytoband, sep="\t", as.is=T,
                         col.names=c("chr", "start", "end", "band", "stain"))
  centromere <- subset(cytoband, stain=="acen" & chr!="chrY", select=-stain)
  cen_gr <- dt2gr(data.frame(start=centromere$start, end=centromere$end,
                             strand=rep('*', nrow(centromere)),
                             seqnames=gsub("chr", "", centromere$chr)))
  cen_gr <- sortSeqlevels(cen_gr)
  seqinfo(cen_gr) <- genome_hg38_info
  sv_gr <- sv_gr[!(sv_gr %over% cen_gr)]
}
message(sprintf("In total %d SVs are used.", length(sv_gr)))

# load covariates
if (args$use_cov) {
  message("Loading covariates")
  data_nt_comp <- readRDS(file_nt_comp)
  # Manually compute GC content, as the method in fishHook tutorial seems to
  # provide two separated covariates rather than G+C
  data_nt_comp$GC_content <- data_nt_comp$G + data_nt_comp$C
  data_rep_timing <- import.bedGraph(file_rep_timing, seqinfo=genome_hg38_info)
  data_dhs <- import.bedGraph(file_dhs, seqinfo=genome_hg38_info)
  data_rmsk_line <- import.bed(file_rmsk_line, seqinfo=genome_hg38_info)
  data_rmsk_sine <- import.bed(file_rmsk_sine, seqinfo=genome_hg38_info)
  data_rmsk_ltr <- import.bed(file_rmsk_ltr, seqinfo=genome_hg38_info)
  data_rmsk_dna <- import.bed(file_rmsk_dna, seqinfo=genome_hg38_info)
  data_rmsk_simple <- import.bed(file_rmsk_simple, seqinfo=genome_hg38_info)
  data_chrom_hmm <- import.bed(file_chrom_hmm, seqinfo=genome_hg38_info)
  # using all repressive states for heterochromatin
  data_chrom_het <- data_chrom_hmm %Q% 
    (name %in% c("Quies", "Het", "ReprPC", "ReprPCWk"))
  data_cfs <- import.bed(file_cfs, seqinfo=genome_hg38_info)
  
  cov_nt_gc <-  Cov(data_nt_comp, name="GC_content", field="GC_content", type="numeric")
  cov_nt_cpg <- Cov(data_nt_comp, name="CpG", field="CG", type="numeric")
  cov_nt_tpc <- Cov(data_nt_comp, name="TpC", field="TC", type="numeric")
  cov_rep_timing <- Cov(data_rep_timing, name="Replication_timing", field="score", type="numeric")
  cov_dhs <- Cov(data_dhs, name="DHS", type="interval")
  cov_rmsk_line <- Cov(data_rmsk_line, name="LINE", type="interval")
  cov_rmsk_sine <- Cov(data_rmsk_sine, name="SINE", type="interval")
  cov_rmsk_ltr <- Cov(data_rmsk_ltr, name="LTR", type="interval")
  cov_rmsk_dna <- Cov(data_rmsk_dna, name="DNA_transposon", type="interval")
  cov_rmsk_simple <- Cov(data_rmsk_simple, name="Simple_repeat", type="interval")
  cov_hetero <- Cov(data_chrom_het, name="Heterochromatin", type="interval")
  cov_fragile <- Cov(data_cfs, name="Fragile_site", type="interval")
}

########## Run fishHook ##########
message("Starting fishHook")
fh <- Fish(hypotheses=chr_window, events=sv_gr, eligible=eligible_reg,
           mc.cores=8)

output_result <- function(fh_model, qqplot_file, list_file, cov_list=NULL,
                          plot_size=6) {
  if (length(cov_list) > 0) {
    fh_model$covariates <- cov_list
  }
  fh_model$score()
  pdf(qqplot_file, height=plot_size, width=plot_size)
  par(mar=c(4,4.5,1,1))
  fh_model$qqp(plotly=FALSE)
  dev.off()
  write.table(fh_model$res %Q% order(p), list_file, sep="\t",
              quote=F, row.names=F)
}

if (args$use_cov) {
  list_cov <- c(
                cov_nt_gc,
                cov_nt_cpg, cov_nt_tpc,
                cov_rmsk_line,
                cov_rmsk_sine, cov_rmsk_ltr, cov_rmsk_dna, cov_rmsk_simple,
                cov_rep_timing, cov_dhs, cov_hetero, cov_fragile
                )
} else {
  list_cov <- NULL
}

fh <- fh[which(fh$data$frac.eligible >= args$frac_eligible), ]
output_result(fh, args$outplot, args$outlist, list_cov)
