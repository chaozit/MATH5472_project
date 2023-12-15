library(SCAVENGE)
library(chromVAR)
library(gchromVAR)
library(SummarizedExperiment)
library(data.table)
library(BiocParallel)
library(BSgenome.Hsapiens.UCSC.hg19)
library(dplyr)
library(igraph)
library(ggplot2)
library(RColorBrewer)
library(ggridges)
library(dplyr)
library(grid)
library(reshape2)
library(RColorBrewer)
library(ggrepel)

all_traits <- c("Lym", "WBC", "Neutro", "Mono", "Eosino", "Baso", "Plt", "RBC", "MCH", "MCHC", "MCV", "Hb")
method <- "xmap"
# method <- "ukb"
# method <- "bbj"
out_all <- data.frame()
for (i in 1:length(all_traits)) {
  trait <- all_traits[i]
  trait_mat <- fread(paste0("/home/share/mingxuan/fine_mapping/analysis/results/scavenge_output_", trait, "_", method, ".txt"))
  out_all <- rbind(out_all, data.frame(trait_mat, trait = all_traits[i]))
  cat(i, "-th trait finished.\n")
}
out_all$BioClassification<-factor(out_all$BioClassification,
                                    levels = c("Late.Eryth", "Early.Eryth", "Early.Baso", "CD14.Mono.1", "GMP.Neut", "cDC", "GMP", "CD14.Mono.2", "HSC", "CMP.LMPP", "CLP.1", "CLP.2", "Plasma", "B", "Pre.B", "pDC", "CD8.EM", "CD8.CM", "CD4.M", "CD8.N", "NK", "CD4.N1", "CD4.N2"))
     


class_avg <- trait_mat %>%
  group_by(BioClassification) %>%
  summarise(
    UMAP1 = median(UMAP1),
    UMAP2 = median(UMAP2)
  )

colourCount = length(unique(trait_mat$BioClassification))
mycolors <- colorRampPalette(brewer.pal(12, "Paired"))(colourCount)


for (i in 1:2) {
  trait_name <- all_traits[i]
  p_umap <- ggplot(data = out_all[out_all$trait == trait_name,], aes(UMAP1, UMAP2, color = TRS)) +
    geom_point(size = 1, na.rm = TRUE, alpha = 0.6) +
    scale_color_gradientn(colors = viridis) +
    scale_alpha() +
    theme_classic() +
    theme(legend.title = element_text(size = 18),
          legend.text = element_text(size = 14),
          legend.position = c(0.8, 0.8),
          axis.title = element_text(size = 18)) +
    xlab("UMAP 1") +
    ylab("UMAP 2") +
    annotate("text", x = -9, y = 13, label = trait_name, size = 8)
  print(p_umap)
  #size: 750*650

  p_bar <- out_all[out_all$trait == trait_name,] %>%
    group_by(BioClassification, trait) %>%
    summarise(enriched_cell = mean(true_cell_top_idx)) %>%
    ggplot(aes(x = BioClassification, y = enriched_cell, fill = BioClassification)) +
    geom_bar(stat = "identity", position = position_dodge()) +
    geom_text(aes(y = enriched_cell, label = scales::percent(enriched_cell, accuracy = 0.1)), position = position_dodge(width = 0.8)) +
    theme_bw() +
    scale_y_continuous(labels = scales::percent_format(scale = 100), limits = c(-1, 1.05), breaks = c(0, 0.25, 0.5, 0.75, 1)) +
    ylab("Trait-enriched cell proportion") +
    xlab(element_blank()) +
    theme(
      panel.border = element_blank(),
      plot.margin = unit(rep(0, 4), "cm"),
      legend.position = "none",
      axis.text.x = element_text(size = 13),
      axis.text.y = element_text(size = 13),
      axis.title.y = element_text(size = 18)
    ) +
    labs(fill = "Cell type") +
    coord_polar(start = 0) +
    annotate("text", x = 0, y = -1, label = trait_name, size = 8)
  #size: 1000*950
  print(p_bar)

}
