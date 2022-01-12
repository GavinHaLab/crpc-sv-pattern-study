
###imports
library(tidyverse)
library(data.table)
library(knitr)
library(survival)
library(ggfortify)
library(survminer)
library(RColorBrewer)


svsigs_annotation_file <- "../sv_cluster_enrichment/outputs/compiled_annotation_data_by_sample.txt"
survival_data_file <- "../../../metadata/WCDT_survival_data.csv"
outname <- "kaplan_meyer_curve_ggsurvplot_with_pval.pdf"


#get survival data 
survival_df <- fread(survival_data_file)

#get sv signautre cluster and alteration information
svsigs_df <- fread(svsigs_annotation_file)

#rename samples to match with the sample names from survival data
svsigs_df$sample <- unlist(lapply(strsplit(svsigs_df$sample,"-"), function(x) paste(x[1], x[2], sep="-")))
svsigs_df$cluster <- paste0("cluster", svsigs_df$cluster)

svsigs_cluster356 <- svsigs_df %>% 
                    filter(cluster == "cluster3" | cluster == "cluster5" | cluster == "cluster6") %>%
                    select(sample,cluster)

mutation_categories <- c("sample","SPOP_mutation", "TP53_mutation", "BRCA2_mutation")
svsigs_alt_df <- svsigs_df[,..mutation_categories]
gather_svsigs_alt <- gather(svsigs_alt_df, detail, value, -sample)
gather_svsigs_alt_complete <- gather_svsigs_alt[complete.cases(gather_svsigs_alt),] %>% select(-value)


######################################################
##choose one of categories to plot survival curve
######################################################

# 1) comparison between mutation types (SPOP, TP53, BRCA2)
# events_df <-  as.data.table(gather_svsigs_alt_complete) 

# 2) comparison between SV clusters
events_df <- svsigs_cluster356 

#merge dataframes
df <- events_df %>%
          merge(survival_df, by.x = "sample", by.y = "Patient.ID")
colnames(df)[2] = "class"

kaplan_meyer_fit <- survfit(Surv(OS.mCRPC.days, Death.Event) ~ class, data = df)

#Log-Rank test comparing survival curves: survdiff()
log_rank_test <- survdiff(Surv(OS.mCRPC.days, Death.Event) ~ class, data = df)


cols <- c("class=BRCA2_mutation"="#E7298A", "class=cluster3"="#1B9E77","class=SPOP_mutation"="#A6761D",
  "class=cluster5"="#E6AB02","class=cluster6"="#D95F02","class=TP53_mutation"="#666666")

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
            font.legend = c(8),
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

