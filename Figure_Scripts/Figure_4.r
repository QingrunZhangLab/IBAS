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

other_results <- fread("./kTWAS-PrediXcan.csv")

kTWAS_genes <- unique(other_results[kTWAS == "Yes", Gene, Disease])
kTWAS_genes$Method <- "kTWAS"
PrediXcan_genes <- unique(other_results[PrediXcan == "Yes", Gene, Disease])
PrediXcan_genes$Method <- "PrediXcan"

IBAS_genes <- unique(wtccc_skat_results[weight_type == "p-value-based", c("phenotype", "Gene", "method")])
colnames(IBAS_genes) <- c("Disease", "Gene", "Method")
IBAS_genes$Method <- paste0("IBAS-", IBAS_genes$Method)

ibas_method_colors <- dim_method_cols
names(ibas_method_colors) <- paste0("IBAS-", names(ibas_method_colors))

all_method_colors <- c(ibas_method_colors, "kTWAS" = "#a0d1b9ff", "PrediXcan" = "#d6a0c4ff")
all_methods_genes <- rbindlist(list(kTWAS_genes, PrediXcan_genes, IBAS_genes))

for(curr_disease in unique(all_methods_genes$Disease)){
  curr_genes <- all_methods_genes[Disease == curr_disease]
  curr_list <- lapply(unique(curr_genes$Method), FUN = function(x){
    curr_genes[Method == x]$Gene
  })
  names(curr_list) <- unique(curr_genes$Method)
  # colorspace::lighten(cell_type_colours[[curr_cell_type]], amount = 0.8, method = "relative", space = "HCL")
  VennDiagram::venn.diagram(curr_list, main = curr_disease,
                            filename = paste0(plot_dir, "/Method_Comparison_", curr_disease, ".tiff"), imagetype = "tiff",
                            cols = all_method_colors[names(curr_list)], fill = all_method_colors[names(curr_list)], alpha = 0.6,
                            margin = 0.15, cex = 1,
                            width = 10, height = 10, units = "cm")
}