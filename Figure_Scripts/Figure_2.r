# Loading necessary packages
pacman::p_load(data.table, stringr, pbapply, argparse, ggplot2, patchwork, disgenet2r, metap, ggrepel, svglite)

# Loading IBAS results
all_simulation_files <- list.files("./Simulations_All", pattern = ".csv", full.names = TRUE, recursive = TRUE)
all_simulation_file_data <- data.table(file_names = all_simulation_files,
                                      dataset = str_split(basename(dirname(all_simulation_files)), "_", simplify = TRUE)[,2],
                                      noise_variance = basename(dirname(dirname(dirname(all_simulation_files)))),
                                      replicate = basename(dirname(dirname(all_simulation_files))),
                                      weight_type = str_split(basename(dirname(all_simulation_files)), "_", simplify = TRUE)[,3])

# Loading PrediXcan results
all_predixcan_files <- list.files("./Simulations_All_PrediXcan_AssociationResults", pattern = ".txt", full.names = TRUE, recursive = TRUE)
all_predixcan_file_data <- data.table(file_names = all_predixcan_files,
                                      dataset = basename(dirname(all_predixcan_files)),
                                      noise_variance = basename(dirname(dirname(dirname(all_predixcan_files)))),
                                      replicate = basename(dirname(dirname(all_predixcan_files))))

all_datasets <- unique(all_simulation_file_data$dataset)

for(curr_dataset in all_datasets){
  # curr_dataset <- all_datasets[1] # For testing

  curr_sim_results <- rbindlist(pblapply(all_simulation_file_data$file_names[all_simulation_file_data$dataset == curr_dataset], function(x){
      skat_results <- fread(x)
      skat_results$noise_variance <- basename(dirname(dirname(dirname(x))))
      skat_results$replicate <- basename(dirname(dirname(x)))
      skat_results$weight_type <- ifelse(str_split(basename(dirname(x)), "_")[[1]][3] == "b", "Coefficient-based", "p-value-based")
      skat_results$dataset <- str_split(basename(dirname(x)), "_")[[1]][2]
      return(skat_results)
  }))
  curr_sim_results$method <- ifelse(curr_sim_results$method == "pca", "IBAS-PCA", 
                                        ifelse(curr_sim_results$method == "tsne", "IBAS-tSNE",
                                          ifelse(curr_sim_results$method == "umap", "IBAS-UMAP", curr_sim_results$method)))
  curr_sim_results$method <- ifelse(curr_sim_results$weight_type == "p-value-based", paste0(curr_sim_results$method, " (p)"), 
                                        ifelse(curr_sim_results$weight_type == "Coefficient-based", paste0(curr_sim_results$method, " (c)"), curr_sim_results$method))
  curr_sim_results$weight_type <- NULL

  curr_predixcan_data_list <- unique(all_predixcan_file_data[all_predixcan_file_data$dataset == curr_dataset, .(noise_variance,replicate)])
  curr_predixcan_results <- rbindlist(pblapply(seq(1, nrow(curr_predixcan_data_list)), function(x){
      curr_results <- rbindlist(lapply(all_predixcan_file_data$file_names[all_predixcan_file_data$dataset == curr_dataset & all_predixcan_file_data$noise_variance == curr_predixcan_data_list$noise_variance[x] & all_predixcan_file_data$replicate == curr_predixcan_data_list$replicate[x]], fread))
      curr_results <- curr_results[is.na(status)]
      curr_results$noise_variance <- curr_predixcan_data_list$noise_variance[x]
      curr_results$replicate <- curr_predixcan_data_list$replicate[x]
      return(curr_results)
  }))

  curr_sim_meta_results <- curr_sim_results[component == "meta"]

  # Identify genes with at least one P.value < 0.05 for each method and noise_variance combination
  all_sim_meta_significant_genes <- curr_sim_meta_results[P.value < 0.05, .(SetID), by = .(method, noise_variance)]

  # Now, non-equi join with the original data.table to include only those genes, for each method and noise_variance
  all_sim_meta_filtered_dt <- curr_sim_meta_results[all_sim_meta_significant_genes, on = .(SetID, method, noise_variance), allow.cartesian = TRUE]

  variance_dt <- all_sim_meta_filtered_dt[, .(Pvalue_variance = var(P.value, na.rm = TRUE)), by = .(SetID, pathway, method, noise_variance)]
  variance_dt <- variance_dt[!is.na(Pvalue_variance)]
  variance_ibas_sim_results <- variance_dt[, .(SetID, method, noise_variance, Pvalue_variance)]

  # Identify genes with at least one P.value < 0.05 for each method and noise_variance combination
  all_predixcan_sim_significant_genes <- curr_predixcan_results[pvalue < 0.05, .(gene), by = .(noise_variance)]

  # Now, non-equi join with the original data.table to include only those genes, for each method and noise_variance
  all_predixcan_sim_filtered_dt <- curr_predixcan_results[all_predixcan_sim_significant_genes, on = .(gene, noise_variance), allow.cartesian = TRUE]

  variance_predixcan_sim_results <- all_predixcan_sim_filtered_dt[, .(Pvalue_variance = var(pvalue, na.rm = TRUE)), by = .(gene, noise_variance)]
  variance_predixcan_sim_results <- variance_predixcan_sim_results[!is.na(Pvalue_variance)]
  variance_predixcan_sim_results$method <- "PrediXcan"
  colnames(variance_predixcan_sim_results) <- c("SetID", "noise_variance", "Pvalue_variance", "method")

  combined_variance <- rbindlist(list(variance_ibas_sim_results, variance_predixcan_sim_results), use.names = TRUE)

  method_cols <- c("IBAS-PCA (p)" = "#fe9ea4ff", "IBAS-tSNE (p)" = "#fec7a7ff", "IBAS-UMAP (p)" = "#a7daf8ff",
                    "IBAS-PCA (c)" = "#fa4b57", "IBAS-tSNE (c)" = "#f69054", "IBAS-UMAP (c)" = "#58befa", 
                    "PrediXcan" = "#d6a0c4ff")

  p_var_plot <- ggplot(combined_variance, aes(x= method, y = Pvalue_variance, fill = method)) +
    geom_boxplot(outlier.size = 0.7) +
    scale_fill_manual(values = method_cols) +
    facet_wrap(~noise_variance) +
    theme_bw() +
    labs(y = "Variance of p-values", title = "Variance of noise", x = "", fill = "Method") +
    theme(text = element_text(size = 9),
          plot.title  = element_text(size = 9, hjust = 0.5),
          axis.ticks.x = element_blank(),
          axis.text.x = element_blank(),
          legend.position = "bottom",
          legend.direction = "horizontal")
    

  # ggsave(plot = p_var_plot, filename = "Simulations_Pvalue_Variance.png", width = 180, height = 100, units = "mm")

  all_sim_meta_results_sig <- curr_sim_meta_results[P.value < 0.05]
  all_predixcan_sim_results_sig <- curr_predixcan_results[pvalue < 0.05]

  # print_log("Calculating the number of genes and pathways discovered for each method and noise variance")
  # Get unique combinations of SetID (gene), method, noise_variance, and replicate
  unique_dt <- unique(all_sim_meta_results_sig, by = c("SetID", "method", "noise_variance", "replicate"))

  # Count the number of unique replicates each gene (SetID) appears in for each combination of method and noise_variance
  counts <- unique_dt[, .N, by = .(SetID, method, noise_variance)]
  counts$N <- counts$N - 1

  # Calculate the percentage
  counts[, percentage := (N / 9) * 100]  # replace 9 with the known number of replicates (-1)

  unique_predixcan <- unique(all_predixcan_sim_results_sig, by = c("gene", "noise_variance", "replicate"))
  predixcan_counts <- unique_predixcan[, .N, by = .(gene, noise_variance)]
  predixcan_counts$N <- predixcan_counts$N - 1
  predixcan_counts[, percentage := (N / 9) * 100]
  predixcan_counts$method <- "PrediXcan"
  colnames(predixcan_counts) <- c("SetID", "noise_variance", "N", "percentage", "method")

  recurrence_data <- rbindlist(list(counts, predixcan_counts), use.names = TRUE)

  recurrence_plot <- ggplot(recurrence_data, aes(x = method, y = percentage, fill = method)) +
    geom_boxplot(outlier.size = 0.7) +
    scale_fill_manual(values = method_cols) +
    ylim(0, 100) +
    theme_bw() +
    facet_wrap(~noise_variance) +
    theme(text = element_text(size = 9),
          axis.title.x = element_blank(),
          axis.text.x = element_text(angle = 45, hjust = 1),
          strip.background = element_blank(),
          strip.text = element_blank(),
          legend.position = "bottom",
          legend.direction = "horizontal") +
    labs(y = "Recurrence Percentage (%)", x = "Variance of noise", fill = "Method")

  # ggsave(plot = recurrence_plot, filename = "Simulations_Recurrence.pdf", path = plot_dir, width = 180, height = 100, units = "mm")

  combined_plot <- p_var_plot + recurrence_plot + plot_layout(nrow = 2, guides = "collect") + plot_annotation(tag_levels = "A") & theme(legend.position = "bottom", plot.tag = element_text(face = "bold"))
  ggsave(plot = combined_plot, filename = paste0(curr_dataset, "_Simulations_Pvalue_Variance.pdf"), path = "./Figures", width = 180, height = 200, units = "mm")
  ggsave(plot = combined_plot, filename = paste0(curr_dataset, "_Simulations_Pvalue_Variance.eps"), path = "./Figures", width = 180, height = 200, units = "mm")
  ggsave(plot = combined_plot, filename = paste0(curr_dataset, "_Simulations_Pvalue_Variance.svg"), path = "./Figures", width = 180, height = 200, units = "mm")
}
