# Loading necessary packages
pacman::p_load(data.table, stringr, pbapply, argparse, ggplot2, patchwork, disgenet2r, ggvenn, metap, ggrepel)

skat_results <- "./03_SKATResults"
plot_dir <- "./Figures"

if(!dir.exists(plot_dir)){
    dir.create(plot_dir)
}

all_skat_results <- rbindlist(pblapply(list.files(skat_results, pattern = ".csv", full.names = TRUE, recursive = TRUE), function(x){
    skat_results <- fread(x)
    skat_results$dataset <- str_split(basename(dirname(x)), "_")[[1]][1]
    skat_results$phenotype <- str_split(basename(dirname(x)), "_")[[1]][2]
    skat_results$weight_type <- ifelse(str_split(basename(dirname(x)), "_")[[1]][3] == "b", "Coefficient-based", "p-value-based")
    return(skat_results)
}))

setnames(all_skat_results, "SetID", "Gene")
method_dict <- data.table(method = c("pca", "tsne", "umap"), method_name = c("PCA", "t-SNE", "UMAP"))
all_skat_results <- merge(all_skat_results, method_dict, by.x = "method", by.y = "method", all.x = TRUE)
all_skat_results[, method := NULL]
setnames(all_skat_results, "method_name", "method")

curr_test_count <- nrow(all_skat_results)
all_skat_results <- all_skat_results[N.Marker.Test >= 3]
filtered_test_count <- nrow(all_skat_results)

all_skat_results[, adj_p := p.adjust(P.value, method = "fdr"), by = .(dataset, phenotype, method, component, pathway)]

curr_test_count <- nrow(all_skat_results)
sig_skat_results <- all_skat_results[adj_p < 0.05]
filtered_test_count <- nrow(sig_skat_results)

wtccc_skat_results <- sig_skat_results[dataset == "WTCCC" & component == "meta"]

dim_method_cols <- c("PCA" = "#fe9ea4ff", "t-SNE" = "#fec7a7ff", "UMAP" = "#a7daf8ff")

disgenet_api_key <- get_disgenet_api_key(
  email = "thalagalakossinnagep@ucalgary.ca", 
  password = getPass::getPass("Enter password for DisGeNET API: "))
Sys.setenv(DISGENET_API_KEY = disgenet_api_key)

CUIs <- fread("DiseaseCUIs.tsv", header = FALSE)
CUIs <- unique(CUIs)

wtccc_skat_results$DisGeNET <- "No"
for(curr_phenotype in unique(CUIs$V1)){
  curated <- disease2gene(disease = CUIs$V2[CUIs$V1 == curr_phenotype], 
                          database = "ALL", verbose = TRUE)
  curated <- extract(curated)
  wtccc_skat_results$DisGeNET[wtccc_skat_results$phenotype == curr_phenotype & wtccc_skat_results$Gene %in% curated$gene_symbol] <- "Yes"
}

wtccc_discovery_results <- unique(wtccc_skat_results[, .(phenotype, weight_type, method, Gene, DisGeNET)])
wtccc_discovery_results$DisGeNET <- ifelse(wtccc_discovery_results$DisGeNET == "Yes", "Rediscovered", "Novel")

for(curr_method in unique(wtccc_discovery_results$method)){
  curr_plot_data <- wtccc_discovery_results[method == curr_method]

  curr_plot <- ggplot(curr_plot_data, aes(x = DisGeNET, fill = DisGeNET)) +
    geom_bar(position = "dodge") +
    scale_fill_manual(values = c("Novel" = "#be3857", "Rediscovered" = "#efbf50")) +
    facet_grid(weight_type~phenotype) +
    theme_bw() +
    theme(axis.text.x = element_blank(),
          axis.ticks.x = element_blank(),
          text = element_text(size = 9),
          legend.position = "bottom") +
    labs(x = "DisGeNET Discoveries", y = "Number of Genes", fill = "Gene type") +
    ggtitle(paste0("WTCCC - ", curr_method, " - Discovery"))

  ggsave(plot = curr_plot, filename = paste0(plot_dir, "/WTCCC_", curr_method, "_DisGeNET.pdf"), width = 15, height = 11, units = "cm")
  ggsave(plot = curr_plot, filename = paste0(plot_dir, "/WTCCC_", curr_method, "_DisGeNET.png"), width = 15, height = 11, units = "cm")
  ggsave(plot = curr_plot, filename = paste0(plot_dir, "/WTCCC_", curr_method, "_DisGeNET.eps"), width = 15, height = 11, units = "cm")
}
