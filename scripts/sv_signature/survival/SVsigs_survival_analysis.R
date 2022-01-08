
###imports
library(tidyverse)
library(data.table)
library(knitr)
library(survival)
library(ggfortify)
library(survminer)
library(RColorBrewer)


svsigs_annotated_file <- "../sv_cluster_enrichment/outputs/sample_by_cluster.txt"
annotated_file <- "../sv_cluster_enrichment/outputs/compiled_annotation_data_by_sample.txt"

survival_data_file <- "../../../metadata/WCDT_survival_data.csv"

###run
#get SVsig data
svsigs_df <- fread(svsigs_annotated_file)
#rename sample
svsigs_df$sample <- unlist(lapply(strsplit(svsigs_df$sample,"-"), function(x) paste(x[1], x[2], sep="-")))
svsigs_df$cluster <- paste0("cluster", svsigs_df$cluster)
filt_svsigs_df <- svsigs_df %>% filter(cluster == "cluster3" | cluster == "cluster5" | cluster == "cluster6")

#get survival data 
survival_df <- fread(survival_data_file)

#annotation data 
annot_df <- fread(annotated_file)

mutation_categories <- c("sample","SPOP_mutation", "TP53_mutation", "BRCA2_mutation")
annot_mut_df <- annot_df[,..mutation_categories]

#rename sample
annot_mut_df$sample <- unlist(lapply(strsplit(annot_mut_df$sample,"-"), function(x) paste(x[1], x[2], sep="-")))

SPOP_alt = annot_mut_df[SPOP_mutation==1,"sample"]
TP53_alt = annot_mut_df[TP53_mutation==1,"sample"]
BRCA2_alt = annot_mut_df[BRCA2_mutation==1,"sample"]

SPOP_df = cbind(SPOP_alt,"SPOP")
TP53_df = cbind(TP53_alt,"TP53")
BRCA2_df = cbind(BRCA2_alt,"BRCA2")


colnames(filt_svsigs_df)[2] = "V2"

######################################################
##choose one of categories to plot
######################################################

#all_events_df = rbind(filt_svsigs_df, SPOP_df,TP53_df, BRCA2_df) #all
#all_events_df = rbind(SPOP_df, TP53_df, BRCA2_df) #mutation only
all_events_df = rbind(filt_svsigs_df) #SV cluster only

#merge dataframes
df <- all_events_df %>%
          merge(survival_df, by.x = "sample", by.y = "Patient.ID")
colnames(df)[2] = "class"

kaplan_meyer_fit <- survfit(Surv(OS.mCRPC.days, Death.Event) ~ class, data = df)

#Log-Rank test comparing survival curves: survdiff()
log_rank_test <- survdiff(Surv(OS.mCRPC.days, Death.Event) ~ class, data = df)

outname="kaplan_meyer_curve_ggsurvplot_with_pval.pdf"

cols = c("class=BRCA2"="#E7298A","class=no_SPOP_TP53_BRCA2"="#66A61E", "class=cluster3"="#1B9E77","class=SPOP"="#A6761D",
  "class=cluster5"="#E6AB02","class=cluster6"="#D95F02","class=TP53"="#666666")

survp <- ggsurvplot(kaplan_meyer_fit,
            pval = TRUE, ## show p-value of log-rank test.
            conf.int = FALSE, # point estimaes of survival curves.
            #risk.table = TRUE, # Add risk table
            conf.int.style = "step",  # customize style of confidence intervals
  		      break.time.by = 500,     # break X axis in time intervals by 200. 
            risk.table.col = "strata", # Change risk table color by groups
            linetype = "strata", # Change line type by groups
            #ncensor.plot = TRUE,  
            ggtheme = theme_bw(), # Change ggplot2 theme
            font.legend = c(10),
            #ggtheme = theme_classic2(base_size=8),
            palette = cols) +
  		      labs(x="Time (days)", y="Overall Survival")


pdf.options(reset = TRUE, onefile = FALSE)
  pdf(file=outname,height=5, width=5)
  print(survp)
dev.off()


res <- pairwise_survdiff(Surv(OS.mCRPC.days, Death.Event) ~ class,
    data = df)
    symnum(res$p.value, cutpoints = c(0, 0.0001, 0.001, 0.01, 0.05, 0.1, 1),
    symbols = c("****", "***", "**", "*", "+", " "),
   abbr.colnames = FALSE, na = "")

