library(ggplot2)
library(ggrepel)
col_pts_border <- 1
col_pts_fill <- "#99999988"

dist_to_zero <- function(x, y) {
  sqrt(x**2 + y**2)
}

dist_to_eq <- function(x, y) {
  abs(x - y)
}

topn_genes <- function(data_x, data_y, all_labels, n=10) {
  labels <- rep(NA, nrow(count_data))
  rank1 <- rank(-dist_to_zero(data_x, data_y), ties.method="random")
  rank2 <- rank(-dist_to_eq(data_x, data_y), ties.method="random")
  idx <- which(rank1 <= n | rank2 <= n)
  labels[idx] <- all_labels[idx]
  labels
}


p_met_vs_pri <- ggplot(
  count_data,
  aes(
      x=sv_mcrpc/n_mcrpc*100,
      y=sv_primary/n_primary*100
      )
  ) +

  xlim(0, 75) +
  ylim(0, 75) +
  labs(
       #title="Primary vs metastatic",
       x="mCRPC samples (%)",
       y="primary samples (%)"
       ) +

   geom_point(
             data=count_data[is.na(count_data$sig_mcrpc_sv) & is.na(count_data$sig_primary_sv) ,], 
             color=col_pts_border,
             fill=col_pts_fill,
             shape=21,
             size =1
             ) +  

  geom_point(data=subset(count_data[!is.na(count_data$sig_mcrpc_sv),],sig_mcrpc_sv <= thres), fill=col_pts_fill, shape=21, size=1.5) +
  geom_point(data=subset(count_data[!is.na(count_data$sig_primary_sv),],sig_primary_sv <= thres), fill=col_pts_fill, shape=21, size=1.5) +
  geom_abline(slope=1, intercept=0, linetype="dashed", colour="#00000080") +

  geom_text_repel(data=subset(count_data[!is.na(count_data$sig_mcrpc_sv)|!is.na(count_data$sig_primary_sv),],sig_primary_sv <= thres | sig_mcrpc_sv <= thres),    
      segment.alpha=0.2, min.segment.length= 0, #Use min.segment.length = 0 to draw all line segments, no matter how short they are.
    aes(label=gene), size=2.5, box.padding=unit(0.35, "lines"), point.padding=unit(0.1, "lines"), nudge_y=2, max.overlaps = 15) +

    theme_classic() +
    theme(axis.title.x = element_text(size=12),
    axis.title.y = element_text(size=12))


 
pdf(sprintf("fig_met-vs-primary_sv_pval_%.2f.pdf", thres), height=2.5, width=2.5)
  print(p_met_vs_pri)
dev.off()
