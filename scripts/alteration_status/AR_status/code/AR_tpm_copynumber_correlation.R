library(ggplot2)
library(grid)
library(gridExtra)
library(reshape2)
# library(sjPlot)
# library(sjmisc)
library(ggrepel)
options(stringsAsFactors = FALSE)

#load batch corrected combined cohort TPM
log_tpm = read.table("../../../metadata/CRPC10X-WCDT_rnaMatrix_log10tpm_commonGenes_batchEffectRemoved.tsv", 
	sep="\t", head=T, check.names=F)
AR_log_exp = log_tpm[log_tpm$hgnc_gs=="AR",2:ncol(log_tpm)]

#load copy number
ARenh_CN = read.table("CRPC10X42_WCDT_AR-enh_compare_AR-vs-Enh.txt",sep="\t",head=T)

samples = ARenh_CN[,"Sample"]

match_id = match(samples, colnames(AR_log_exp))
del_id = which(is.na(match_id))

ARenh_CN_tpm_avail = ARenh_CN[-del_id,]
ARenh_CN_tpm_avail$tpm = unlist(AR_log_exp[na.exclude(match_id)])

###### x-axis AR copy number y-axis Enhancer copy number dot as expression value
#for lower copies
gp <- ggplot(ARenh_CN_tpm_avail, aes(x=AR.norm, y=Enhancer.norm, fill=tpm)) +
    geom_abline(slope=1, intercept=0, color="red") +
    geom_point(shape=21, size=4, alpha=0.75) +
    scale_fill_distiller(palette = "Spectral")+
    ylim(c(0.5, 3.55)) + xlim(c(0.5,3.55)) + ##for low_copy
    labs(fill = "log10(TPM+1)") +
    theme_bw() + coord_fixed() +
    theme(axis.text.y=element_text(size=20), axis.text.x=element_text(size=16),
          axis.title = element_blank(), plot.title=element_text(size=24, face="bold"))

ggsave(gp, file="AR-vs-Enh_with_logTPM_low_copy.pdf", width=5, height=5)

gp <- ggplot(ARenh_CN_tpm_avail, aes(x=AR.norm, y=Enhancer.norm, fill=tpm)) +
    geom_abline(slope=1, intercept=0, color="red") +
    geom_point(shape=21, alpha=0.75, aes(size = exp(tpm))) + 
    geom_text_repel(aes(label=Sample), color="black", size=2) +
    scale_fill_distiller(palette = "Spectral")+
    ylim(c(0.5, 3.55)) + xlim(c(0.5,3.55)) + ##for low_copy
    labs(fill = "log10(TPM+1)") +
    theme_bw() + coord_fixed() +
    theme(axis.text.y=element_text(size=20), axis.text.x=element_text(size=16),
          axis.title = element_blank(), plot.title=element_text(size=24, face="bold"))

ggsave(gp, file="AR-vs-Enh_with_logTPM_low_copy_with_sample_names.pdf", width=5, height=5)

#for all copies
gp <- ggplot(ARenh_CN_tpm_avail, aes(x=AR.norm, y=Enhancer.norm, fill=tpm)) +
    geom_abline(slope=1, intercept=0, color="red") +
    geom_point(shape=21, size=4, alpha=0.75) +
    scale_fill_distiller(palette = "Spectral")+
    scale_x_continuous(limits=c(0.5,280), breaks=c(1,2,8,32,128), labels=c(1,2,8,32,128), trans="log2") + 
    scale_y_continuous(limits=c(0.5,280), breaks=c(1,2,8,32,128), labels=c(1,2,8,32,128), trans="log2") +    
    theme_bw() + coord_fixed() +
    labs(fill = "log10(TPM+1)") +
    theme(axis.text.y=element_text(size=20), axis.text.x=element_text(size=16),
          axis.title = element_blank(), plot.title=element_text(size=24, face="bold"))

ggsave(gp, file="AR-vs-Enh_with_logTPM.pdf", width=5, height=5)

gp <- ggplot(ARenh_CN_tpm_avail, aes(x=AR.norm, y=Enhancer.norm, fill=tpm)) +
    geom_abline(slope=1, intercept=0, color="red") +
    geom_point(shape=21, size=4, alpha=0.75) +
    geom_text_repel(aes(label=Sample), color="black", size=2) +
    scale_fill_distiller(palette = "Spectral")+
    scale_x_continuous(limits=c(0.5,280), breaks=c(1,2,8,32,128), labels=c(1,2,8,32,128), trans="log2") + 
    scale_y_continuous(limits=c(0.5,280), breaks=c(1,2,8,32,128), labels=c(1,2,8,32,128), trans="log2") +    
    theme_bw() + coord_fixed() +
    labs(fill = "log10(TPM+1)") +
    theme(axis.text.y=element_text(size=20), axis.text.x=element_text(size=16),
          axis.title = element_blank(), plot.title=element_text(size=24, face="bold"))

ggsave(gp, file="AR-vs-Enh_with_logTPM_with_sample_names.pdf", width=5, height=5)

######################
#run multiple linear regression afterwards
######################
# #fit <- lm(tpm ~ ploidy + purity + AR.norm + Enhancer.norm + AR.norm*Enhancer.norm, data=ARenh_CN_tpm_avail)
# fit <- lm(tpm ~ ploidy + purity + AR + Enhancer + AR*Enhancer, data=ARenh_CN_tpm_avail)
# summary(fit) # show results

# # fit <- glm(tpm ~ ploidy + purity + AR + Enhancer + AR*Enhancer, data=ARenh_CN_tpm_avail)
# # summary(fit) # show results

# # fit model with interaction
# pdf("model_fit_estimates.pdf", width=6, height=4)
#     plot_model(fit, sort.est = TRUE,show.values = TRUE, title="log10(TPM)")
# dev.off()

# pdf("model_interaction_terms_copy_0_4.pdf", width=6, height=4)
#     #plot_model(fit, type = "pred", terms = c("AR", "Enhancer [0,2,4]"))
#     plot_model(fit, type = "pred", terms = c("Enhancer [0,2,4]", "AR [0,2,4]"), title="Predicted values of log10(TPM)")
#     #plot_model(fit, type = "pred", terms = c("AR [0,2,4]", "Enhancer [0,2,4]"))
# dev.off()


# pdf("model_interaction_terms_expanded_copies.pdf", width=6, height=4)
#     #plot_model(fit, terms = c("AR", "Enhancer"),type = "int", mdrt.values = "minmax" )
#     #plot_model(fit, terms = c("AR", "Enhancer"), type = "pred")
#     plot_model(fit, type = "pred", terms = c("AR", "Enhancer [30, 50, 70]"), title="Predicted values of log10(TPM)")
# dev.off()

