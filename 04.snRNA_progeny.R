###############################################################
# Project: scRNA-seq Cardiomyocyte PROGENy Pathway Analysis
# Author: maihuanzhuo
# Date: 2026-03-27
# Description:
# This script performs PROGENy pathway activity analysis on
# cardiomyocytes from Angiotensin and BPF_Angiotensin groups.

############################
## 1. Load packages
############################
if(!require(Seurat))BiocManager::install("Seurat") 
if(!require(ggplot2))BiocManager::install("ggplot2")
if(!require(dplyr))BiocManager::install("dplyr")
if(!require(stringr))BiocManager::install("stringr")
if(!require(tidyverse))BiocManager::install("tidyverse")
if(!require(patchwork))BiocManager::install("patchwork")
if(!require(qs))BiocManager::install("qs")
if(!require(progeny))BiocManager::install("progeny")
if(!require(limma))BiocManager::install("limma")

setwd("E:/Omics_Data/scRNA/results_scanpy")

############################
## 2. Load data
############################

sce.all.filt <- qs::qread("./results/2-DA/sce.all.filt.qs")
sce.all.filt
cm <- subset(sce.all.filt, subset = cell_type == "Cardiomyocytes")
# 
metadata <- sce.all.filt@meta.data %>% 
  dplyr::select(sample, group) %>% 
  rownames_to_column(var = "barcode")
# 
metadata_cm <- cm@meta.data %>% 
  dplyr::select(sample, group) %>% 
  dplyr::filter(group == "Angiotensin" | group == "BPF_Angiotensin") %>% 
  mutate(group = factor(group, levels = c("Angiotensin", "BPF_Angiotensin"))) %>% 
  rownames_to_column(var = "barcode") %>% 
  dplyr::select(-barcode) %>% 
  distinct(sample, .keep_all = TRUE)

############################
## 3. PROGENy scoring
############################

progeny_score <- progeny(cm, top = 500, organism = "Mouse", scale = F, 
                         perm = 1, return_assay = TRUE)
progeny <- ScaleData(progeny_score, assay = "progeny")

############################
## 4. Sample-level aggregation
############################
# Extract raw PROGENy pathway scores from the assay matrix
progeny_score_mat <- progeny@assays$progeny$data %>% 
  t() %>% 
  as.data.frame() %>% 
  rownames_to_column(var = "barcode") %>% 
  left_join(metadata, by = "barcode") %>% 
  dplyr::filter(group == "Angiotensin" | group == "BPF_Angiotensin") %>% 
  group_by(sample, group) %>%
  # group_by(group) %>%
  summarise(across(where(is.numeric), mean), .groups = "drop") %>% 
  dplyr::select(-group) %>%
  column_to_rownames(var = "sample") %>%
  # column_to_rownames(var = "group") %>% 
  t() %>% 
  as.matrix()
############################
## 5. Heatmap visualization
############################

matrix_scale <- t(scale(t(progeny_score_mat)))
# 
annot_df <- data.frame(sample = colnames(matrix_scale), 
                       Group = c(rep("Ang",5), rep("Ang+BPF_50", 5)))
annot_df$Group <- factor(annot_df$Group, levels = c("Ang", "Ang+BPF_50"))
col_fun <- circlize::colorRamp2(c(-2, 0, 2), rev(RColorBrewer::brewer.pal(11, "RdBu"))[c(1,6,11)])
# annotation
col_ha <- HeatmapAnnotation(Group = factor(annot_df$Group, levels = c("Ang", "Ang+BPF_50")),
                            col = list(Group = c("Ang+BPF_50" = "#F7A944", "Ang" = "#ACACAC")),
                            gp = gpar(col = "black", lwd = 2), 
                            border = T,
                            show_legend = F,
                            annotation_name_gp = gpar(fontface = "bold"))
Legend_1 <- Legend(labels = c("Ang", "Ang+BPF_50"),
                   legend_gp = gpar(fill = c("#ACACAC","#F7A944"), col = "black", lwd = 1.5), 
                   title = "Group", border = T, row_gap = unit(0.5, "mm"))
Legend_2 <- Legend(col_fun = col_fun, title = "scale_expression", at = c(-2, -1, 0, 1, 2), border = T,
                   legend_height = unit(4, "cm"), legend_width = unit(2, "cm"))
Legend_list <- packLegend(Legend_1, Legend_2)
# 
labs <- rownames(matrix_scale)
# 
ha <-  rowAnnotation(foo = anno_mark(at = 1:14,
                                     labels = labs,
                                     labels_gp = gpar(fontsize = 10, fontface = 'bold'),
                                     link_width = unit(3, "mm"),
                                     link_gp = gpar(col = "black", lwd = 2, lty = 1)))
# plot
ht <- Heatmap(matrix_scale, 
              name = "scale_expression",
              col = col_fun,
              show_row_names = F, 
              show_column_names = F,
              cluster_rows = TRUE,
              cluster_columns = F,
              show_row_dend = F,
              top_annotation = col_ha,
              right_annotation = ha, 
              row_title = NULL, 
              column_title_gp = gpar(fontface = "bold", col = "black"),
              cluster_column_slices = F, 
              
              border = TRUE,
              border_gp = gpar(col = "black", lwd = 2),
              column_dend_gp = gpar(lwd = 2),
              show_heatmap_legend = F,
              heatmap_legend_param = list(border = T, legend_gp = gpar(col = "black", lwd = 1.5),
                                          legend_height = unit(4, "cm"),
                                          legend_width = unit(2, "cm"))
)
ht
draw(ht, annotation_legend_list = Legend_list)
pdf("./complexheatmap.pdf", width = 5, height = 4)
draw(ht, annotation_legend_list = Legend_list)
dev.off()

############################
## 6. Differential analysis (limma)
############################
# Extract scaled PROGENy pathway scores from the assay matrix
progeny_score_mat_scale <- progeny@assays$progeny$scale.data %>% 
  t() %>% 
  as.data.frame() %>% 
  rownames_to_column(var = "barcode") %>% 
  left_join(metadata, by = "barcode") %>% 
  dplyr::filter(group == "Angiotensin" | group == "BPF_Angiotensin") %>% 
  group_by(sample, group) %>%
  # group_by(group) %>%
  summarise(across(where(is.numeric), mean), .groups = "drop") %>% 
  dplyr::select(-group) %>%
  column_to_rownames(var = "sample") %>%
  # column_to_rownames(var = "group") %>% 
  t() %>% 
  as.matrix()

design <- model.matrix(~ 0 + factor(metadata_cm$group))
colnames(design) <- levels(factor(metadata_cm$group))
rownames(design) <- colnames(progeny_score_mat_scale)
# 
contrast.matrix <- makeContrasts(BPF_Angiotensin - Angiotensin, levels = design)
# 
fit <- lmFit(progeny_score_mat_scale, design)
fit2 <- contrasts.fit(fit, contrast.matrix)
fit2 <- eBayes(fit2)
diff <- topTable(fit2, coef = 1, n = Inf, adjust.method = "BH", sort.by = "logFC")
head(diff)
write.csv(diff, "./progeny_diff_stat.csv")

############################
## 7. Barplot visualization
############################

diff <- read.csv("./progeny_diff_stat.csv")
# 
df_plot <- diff %>% 
  dplyr::mutate(
    pathway_full = dplyr::recode(
      pathway,
      "p53"      = "p53 Signaling Pathway",
      "VEGF"     = "Vascular Endothelial Growth Factor Signaling Pathway",
      "PI3K"     = "PI3K Signaling Pathway",
      "TGFb"     = "TGF-β Signaling Pathway",
      "Androgen" = "Androgen Receptor Signaling Pathway",
      "WNT"      = "Wnt Signaling Pathway",
      "EGFR"     = "Epidermal Growth Factor Receptor Signaling Pathway",
      "MAPK"     = "MAPK Signaling Pathway",
      "Estrogen" = "Estrogen Receptor Signaling Pathway",
      "Trail"    = "TNF-Related Apoptosis Inducing Ligand Signaling Pathway",
      "TNFa"     = "TNF-α Signaling Pathway",
      "JAK-STAT" = "JAK-STAT Signaling Pathway",
      "NFkB"     = "NFkB Signaling Pathway",
      "Hypoxia"  = "HIF-1 Signaling Pathway"
    )
  )
# 
p <- ggplot(df_plot, aes(x = reorder(pathway, t), y = t, fill = logFC > 0)) +
  geom_col() + 
  geom_text(data = df_plot %>% dplyr::filter(t > 0), 
            aes(x = reorder(pathway, logFC), y = -0.1, label = pathway_full), 
            hjust = 1, size = 3.3, fontface = "bold") +
  geom_text(data = df_plot %>% dplyr::filter(t < 0), 
            aes(x = reorder(pathway, logFC), y = 0.1, label = pathway_full),
            hjust = 0, size = 3.3, fontface = "bold") +
  geom_segment(aes(y = -0.1, yend = -2.5, x = 15, xend = 15),
               arrow = arrow(length = unit(0.2, "cm"), type = "closed"),
               size = 0.5) +
  geom_segment(aes(y = 0.05, yend = 2.5, x = 15, xend = 15),
               arrow = arrow(length = unit(0.2, "cm"), type = "closed"),
               size = 0.5) +
  annotate("text", x = 15, y = -3.3, label = "Down", fontface = "bold", size = 5, color = "#364C75") +
  annotate("text", x = 15, y = 3.3, label = "UP", fontface = "bold", size = 5, color = "#BA0404") +
  scale_fill_manual(values = c("#364C75","#BA0404"))+
  ylim(-5, 5) +
  labs(x = "", y = "t value of PROGENy score", title = "Ang+BPF_50 vs. Ang") +
  ggprism::theme_prism(base_size = 10, base_line_size = 0.8, base_rect_size = 0.8) +
  theme(legend.position = "none", 
        axis.title.x = ggtext::element_markdown(size = 12),
        axis.title.y = ggtext::element_markdown(size = 12),
        axis.line.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.text.y = element_blank()) + 
  coord_flip(clip = "off")
ggsave(width = 8, height = 5, dpi = 600, p,
       filename = "E:/Omics_Data/scRNA/results/progeny/bar_progeny_score.pdf")

