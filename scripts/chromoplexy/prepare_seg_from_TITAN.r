library(argparse)
parser <- ArgumentParser(description="Prepare a SEG file from TITAN result.")
parser$add_argument("-i", "--in", dest="input", required=T,
  help="Input file from TITAN combined with ichorCNA no SNP")
parser$add_argument("-o", "--out", dest="output", required=T,
  help="Output SEG file")
parser$add_argument("-n", "--zero-neutral", action="store_true", dest="zero_neut",
  default=F, help="Turn on to set neutral calls to 0 logR")
parser$add_argument("--keep-homd", action="store_true", dest="keep_homd",
  default=F, help="Turn on to keep HOMD segments")
parser$add_argument("--max-homd-len", type="integer", default=100000,
  help="Max segment length for HOMD. Only works when NOT setting --keep-homd. Default: %(default)d",
  dest="max_homd_len")
parser$add_argument("--max-homd-snp", type="integer", default=40,
  help="Max SNP length for HOMD. Only works when NOT setting --keep-homd. Default: %(default)d",
  dest="max_homd_snp")
parser$add_argument("--keep-short", action="store_true", dest="keep_short",
  default=F, help="Turn on to keep short segments")
parser$add_argument("--min-len", type="integer", default=50000,
  help="Min segment length to be kept. Only works when NOT setting --keep-short. Default: %(default)d",
  dest="min_len")
parser$add_argument("--min-snp", type="integer", default=4,
  help="Min SNP length to be kept. Only works when NOT setting --keep-short. Default: %(default)d",
  dest="min_snp")
parser$add_argument("--no-subclone", action="store_true", dest="no_subclone",
  default=F, help="Turn on to exclude subclone segments")
args <- parser$parse_args()

library(data.table)
segInput <- fread(args$input)
ind.chrX <- grepl("X", segInput$Chromosome)
segInput[!ind.chrX, Corrected_logR := log2(logR_Copy_Number / 2)]
segInput[ind.chrX,  Corrected_logR := log2(logR_Copy_Number / 1)]
segInput[Corrected_logR < -1.5, Corrected_logR := -1.5]

if (args$zero_neut) {
  #segInput[!ind.chrX & Corrected_Copy_Number == 2, Corrected_logR := 0]
  #segInput[ind.chrX  & Corrected_Copy_Number == 1, Corrected_logR := 0]
  segInput[!ind.chrX & Corrected_Call == "NEUT", Corrected_logR := 0]
  segInput[ind.chrX  & Corrected_Call == "NEUT", Corrected_logR := 0]
}
ind <- 1:nrow(segInput)
if (!args$keep_homd){
  exclude.homd.ind <- segInput[Corrected_Call == "HOMD"
    & (Length.snp. < args$max_homd_snp
       | (End-Start) < args$max_homd_len),
    which = TRUE]
  ind <- setdiff(ind, exclude.homd.ind)
}
if (!args$keep_short){
  exclude.size.ind <- segInput[Length.snp. < args$min_snp
    | (End-Start) < args$min_len, which = TRUE]
  ind <- setdiff(ind, exclude.size.ind)
}
if (args$no_subclone){
  exclude.subclone.ind <- segInput[Cellular_Prevalence < 1, which = TRUE]
  ind <- setdiff(ind, exclude.subclone.ind)
}
segInput <- segInput[ind]
igvSegs <- unique(segInput[, .(Sample, Chromosome, Start, End, Length.snp., Corrected_logR)])
fwrite(igvSegs, file=args$output, col.names=T, row.names=F, quote=F, sep="\t")
