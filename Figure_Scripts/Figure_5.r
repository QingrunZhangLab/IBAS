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

pca_p_value_results <- wtccc_skat_results[method == "PCA" & weight_type == "p-value-based"]
pca_p_value_results$DisGeNET <- ifelse(pca_p_value_results$DisGeNET == "Yes", "Rediscovered", "Novel")

summarized_pca_p_value_results <- pca_p_value_results[ , .(Combined_P = ifelse(length(adj_p)<2, sumlog(adj_p)$validp, sumlog(adj_p)$p), Num_Pathways = .N), by = .(Gene, phenotype)]

disgenet_results <- unique(pca_p_value_results[, .(phenotype, Gene, DisGeNET)])
summarized_pca_p_value_results <- summarized_pca_p_value_results[disgenet_results, on = .(phenotype, Gene)]

setorder(summarized_pca_p_value_results, phenotype, Combined_P)
summarized_pca_p_value_results[, gene_label := ifelse(frank(Combined_P, ties.method = "min") <= 10, Gene, NA), by = phenotype]

result_summary <- ggplot(summarized_pca_p_value_results, aes(x = Num_Pathways, y= -log(Combined_P), color = DisGeNET)) + 
                    geom_point(size = 0.5) +
                    facet_wrap(~phenotype, scales = "free_y", strip.position = "top") +
                    labs(x = "Number of Mediated Pathways", y = "-log(Combined P-value)", color = "DisGeNET") +
                    scale_color_manual(values = c("Novel" = "#be3857", "Rediscovered" = "#4aaf3a")) +
                    theme_bw() +
                    theme(text = element_text(size = 9),
                          legend.position = "bottom") + 
                    geom_text_repel(aes(label = gene_label), size = 2, show.legend = FALSE, max.overlaps = Inf)

ggsave(plot = result_summary, filename = paste0(plot_dir, "/WTCCC_PCA_p-value-based_DisGeNET.pdf"), width = 20, height = 11, units = "cm")
ggsave(plot = result_summary, filename = paste0(plot_dir, "/WTCCC_PCA_p-value-based_DisGeNET.eps"), width = 20, height = 11, units = "cm")
ggsave(plot = result_summary, filename = paste0(plot_dir, "/WTCCC_PCA_p-value-based_DisGeNET.png"), width = 20, height = 11, units = "cm")