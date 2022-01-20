
# Join consecutives segments with the same "Corrected_Copy_Number"
# Required as copy number input for shatterseek analysis (run_ShatterSeek.Rmd)

library(dplyr)

gap_threshold <- 1

input_rds <- "../../metadata/CRPC10X42-WCDT101_copy_number.rds"
output_prefix <- gsub(".rds","",basename(input_rds))

# copy number input file; requires 4 columns ('Sample','Chromosome','Start','End','Corrected_Copy_Number')
CNdata <- readRDS(input_rds)

centromere <- read.table("../../metadata/GRCh38.GCA_000001405.2_centromere_acen.txt", header = TRUE, sep='\t')

## merge segments based on same copy number
joined_CNdata <- CNdata[1,]  # init starting sample
currSample <- joined_CNdata[,'Sample']
currEnd <- joined_CNdata[,'End']
currChr <- joined_CNdata[,'Chromosome']
currCN <- joined_CNdata[,'Corrected_Copy_Number']

for (i in 2:nrow(CNdata)){
	if (CNdata[i,'Sample'] == currSample & 
	    CNdata[i,'Chromosome'] == currChr & 
	    CNdata[i,'Corrected_Copy_Number']== currCN & CNdata[i,'Start']-currEnd <= gap_threshold){
		currEnd <- CNdata[i,'End']
	} else{
		joined_CNdata[nrow(joined_CNdata),'End'] <- currEnd
		joined_CNdata <- rbind(joined_CNdata, CNdata[i,])
		currSample <- CNdata[i,'Sample']
		currEnd <- CNdata[i,'End']
		currChr <- CNdata[i,'Chromosome']
		currCN <- CNdata[i,'Corrected_Copy_Number']
	}
	if (i %% 1000 == 0) {
	  print(i)
	}
}
joined_CNdata[nrow(joined_CNdata),'End'] <- currEnd

##remove centromere
keepInd <- !logical(length = length(joined_CNdata$Chromosome))
for (c in centromere$Chr){
	ind <- which((joined_CNdata[joined_CNdata$Chromosome==c,'Start'] >= centromere[centromere$Chr==c, "Start"]) & (joined_CNdata[joined_CNdata$Chromosome==c,'End'] <= centromere[centromere$Chr==c, "End"]))
	keepInd[ind] <- FALSE
}

joined_CNdata <- joined_CNdata[which(keepInd), ]
write.table(joined_CNdata, sprintf('%s.joined.txt', output_prefix), sep='\t', row.names=F, col.names=T, quote=F)

