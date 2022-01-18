args <- commandArgs(trailingOnly=T)
suppressMessages(library(fishHook))
suppressMessages(library(BSgenome.Hsapiens.NCBI.GRCh38))
genome_hg38 <- BSgenome.Hsapiens.NCBI.GRCh38
genome_hg38_info <- seqinfo(genome_hg38)[seqlevels(genome_hg38)[1:23]]

# *NOTE*: BIN_SIZE and BIN_STEP in kb
BIN_SIZE <- as.integer(args[1])
# use non-overlapping bins
OUT_FILE <- sprintf("nt_composition.tile%dkb.hg38.rds", BIN_SIZE)
chr_window <- tileGenome(genome_hg38_info, tilewidth=BIN_SIZE*1e3,
                         cut.last.tile.in.chrom=T)
seq <- getSeq(genome_hg38, chr_window)
comp <- sapply(1:3, oligonucleotideFrequency, x=seq, as.prob=T)
mcols(chr_window) <- comp
saveRDS(chr_window, OUT_FILE)
