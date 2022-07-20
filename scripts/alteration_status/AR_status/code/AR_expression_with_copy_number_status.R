
# This will produce main figure 3-B 
# Also fit the ANCOVA model using AR expression as the response variable, AR amplification status as the predictor variable, and ploidy, purity as the covariate.

library(ggplot2)
library(ggpubr)
library(multcomp)
library(car)
library(tidyr)

options(stringsAsFactors=F)

n_fun <- function(x){
  return(data.frame(y = max(x), label = paste0("n = ",length(x))))
}

# This is output from compareCN_2regions_v2.R
ARenh_CN <- read.table("CRPC10X42_WCDT_AR-enh_compare_AR-vs-Enh.txt", sep="\t", head=T)
log10_tpm <- read.table("../../../metadata/CRPC10X-WCDT_rnaMatrix_log10tpm_commonGenes_batchEffectRemoved.tsv", sep="\t", head=T, check.names=F)

AR_exp <- log10_tpm[log10_tpm$hgnc_gs=="AR", 2:ncol(log10_tpm)]

#wide to long
AR_exp_gather <- gather(AR_exp, Sample, AR_logTPM)

#merge ARenh copy nubmer and AR expression information
ARenh_CN_exp_combined <- AR_exp_gather %>% merge(ARenh_CN, by.x="Sample", by.y="Sample")

# Reorder the groups order
ARenh_CN_exp_combined$call <- factor(ARenh_CN_exp_combined$call , levels=c("No amplification","Selective AR","Selective Enhancer","Coamplification"))

cols <- c("No amplification"="grey", "Coamplification"="red", "Selective Enhancer"="#E69F00", "Selective AR"="#009E73")

#wilcox test between 4 groups
#wilcox_results <- compare_means(AR_logTPM ~ call, data = ARenh_CN_exp_combined, method = "wilcox.test")

pdf("boxplot_for_AR_expression_CN.pdf", width=6, height=4, useDingbats=FALSE)
    ggboxplot(ARenh_CN_exp_combined, x = "call", y = "AR_logTPM",
          color = "call",  add = "jitter", palette =c("grey", "#009E73", "#E69F00","red"), add.params = list(size = 1.5))+
    theme_bw() +
    stat_summary(fun.data = n_fun, geom = "text", position = position_dodge(width = 1), size = 3) +
    ylab("AR expression [log10(TPM+1)]")+
    theme(axis.text.y=element_text(size=12), axis.text.x=element_text(size=8), axis.title.x = element_blank())
dev.off()


#perform ANCOVA anslysis
#https://www.statology.org/ancova-in-r/
leveneTest(AR_logTPM ~ call, data=ARenh_CN_exp_combined)

ancova_model <- aov(AR_logTPM ~ ploidy + purity + call, data = ARenh_CN_exp_combined)
Anova(ancova_model, type="III") 
# Anova Table (Type III tests)

# Response: AR_logTPM
#             Sum Sq  Df F value    Pr(>F)
# (Intercept) 16.898   1 46.4761 4.411e-10 ***
# ploidy       0.295   1  0.8108    0.3697
# purity       0.344   1  0.9457    0.3328
# call        19.982   3 18.3192 8.495e-10 ***
# Residuals   42.177 116

postHocs <- glht(ancova_model, linfct = mcp(call = "Tukey"))
summary(postHocs)

#view the confidence intervals associated with the multiple comparisons
confint(postHocs)
