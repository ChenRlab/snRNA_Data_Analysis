# devtools::install_github("Bioconductor/MatrixGenerics")
# devtools::install_github("const-ae/sparseMatrixStats")
# devtools::install_github("neurorestore/Augur")
library(Augur)
library(Seurat)
library(viridis)

# 
augur <- calculate_auc(sce,
                       cell_type_col = "cell_type",
                       label_col = "group", 
                       n_threads = 8,
                       n_subsamples = 50,
                       subsample_size = 20,
                       folds = 3,
                       var_quantile = 0.5, 
                       feature_perc = 0.5, 
                       augur_mode = c("default", "velocity", "permute"), 
                       classifier = c("rf", "lr"), 
                       rf_params = list(trees = 100, mtry = 2, min_n = NULL, importance = "accuracy"),
                       lr_params = list(mixture = 1, penalty = "auto"))
augur$AUC
qs::qsave(augur, file = "./results/2-DA/Augur/augur_con_BPF.qs")
# 
plot_augur <- augur$AUC %>% 
  as.data.frame() %>% 
  arrange(desc(auc)) %>% 
  mutate(cell_type = factor(cell_type, levels = rev(unique(cell_type))))
# 
p_augur <- ggplot(data = plot_augur, x = auc, y = cell_type) + 
  ggforce::geom_link(aes(x = 0.5, xend = auc, y = cell_type, yend = cell_type, color = cell_type,
                         alpha = after_stat(index), size = after_stat(index)),
                     n = 800, show.legend = F) +
  geom_point(aes(x = auc, y = cell_type, color = cell_type, size = auc),
             fill = "white", shape = 21) + 
  # geom_text(aes(x = auc + 0.013, y = cell_type, label = sprintf("%.3f", auc)), size = 3, fontface = "bold") +
  # geom_text(aes(x = auc + 0.018, y = cell_type, label = sprintf("%.3f", auc)), size = 3, fontface = "bold") +
  geom_text(aes(x = auc + 0.010, y = cell_type, label = sprintf("%.3f", auc)), size = 3, fontface = "bold") +
  geom_vline(xintercept = 0.5, linetype = "dashed", color = "grey", size = 0.8) + 
  scale_color_manual(values = row_name_col) +
  xlab("AUC") + ylab("Cell type") + 
  ggprism::theme_prism(base_size = 12, border = T) +
  theme(legend.position = "none")
p_augur
ggsave(p_augur, width = 6, height = 5, dpi = 600, file = "./results/Augur/p_augur.pdf")

