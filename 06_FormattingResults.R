# Uses the pacman package to ensure all required packages are installed and loaded
if (!require("pacman")) install.packages("pacman")
if (!require("devtools")) install.packages("devtools")
library(devtools)
if (!require("disgenet2r")) install_bitbucket("ibi_group/disgenet2r")
pacman::p_load(data.table, 
               ggplot2, 
               disgenet2r,
               optparse,
               ggsci,
               stringr,
               ggvenn,
               gridExtra,
               patchwork)

# Defining arguments to be parse from command line
option_list <- list(
  make_option(c("-r", "--results_file"), action = "store", type = "character", 
              default = "Pathway_SKAT_Results.csv",
              help="Location to load final results. Default location is %default."),
  make_option(c("-P", "--pathway_availability_file"), action = "store", type = "character", 
              default = "./Reference_Pathway_Gene_Availability_Summary.csv",
              help = "Location to load pathway availability summary data. Default location is %default."),
  make_option(c("-c", "--disease_cui"), action = "store", type = "character", 
              default = "C0036341",
              help = "CUI for disease from DisGeNet database. Default is %default (for Schizophrenia)."),
  make_option(c("-o", "--output_results_folder"), action = "store", type = "character", 
              default = "./Pathway_SKAT_Results_Output/",
              help = "Location to save final results. Default location is %default.")
)

# Loading parsed arguments
args <- parse_args(OptionParser(option_list = option_list))

if(!dir.exists(args$output_results_folder)){
  dir.create(args$output_results_folder, recursive = TRUE)
}

disgenet_api_key <- get_disgenet_api_key(
  email = "thalagalakossinnagep@ucalgary.ca", 
  password = "e8wt3h35$")
Sys.setenv(DISGENET_API_KEY = disgenet_api_key)

results <- fread(args$results_file)

pathway_data <- fread(args$pathway_availability_file)
req_data <- pathway_data[, c("PathwayID", "PathwayClass", "PathwayName")]

results <- merge(results, req_data, by.x = "pathway", by.y = "PathwayID")
results <- results[, c("PathwayClass", "PathwayName", "pathway", "SetID", "P.value",
                       "N.Marker.All", "N.Marker.Test", "num_snp_sets", "method_dim")]
colnames(results) <- c("PathwayClass", "PathwayName", "PathwayID", "Gene", "p.value", 
                       "N.Marker.All", "N.Marker.Test", "GenesInPathway", "Dimensionality Reduction")

gene_counts <- results[, .(gene_counts = length(Gene)), by = .(PathwayID, `Dimensionality Reduction`)]
results <- merge(results, gene_counts, by = c("PathwayID", "Dimensionality Reduction"))

results$adj.P.value <- 1
for(x in unique(results$PathwayID)){
  for(y in unique(results$`Dimensionality Reduction`)){
    num_comparisons <- results$gene_counts[results$PathwayID == x & results$`Dimensionality Reduction` == y]
    results$adj.P.value[results$PathwayID == x & results$`Dimensionality Reduction` == y] <- p.adjust(results$p.value[results$PathwayID == x & results$PathwayID == x & results$`Dimensionality Reduction` == y],
                                                                           method = "fdr")
  }
}

results$is.Sig <- ifelse(results$adj.P.value < 0.05, "Yes", "No")

disgenet_api_key <- get_disgenet_api_key(
  email = "thalagalakossinnagep@ucalgary.ca", 
  password = "e8wt3h35$")
Sys.setenv(DISGENET_API_KEY = disgenet_api_key)

results$Curated <- "No"
results$All <- "No"

curated <- disease2gene(disease = "C0005586", 
                        database = "CURATED", verbose = TRUE)
curated <- extract(curated)
all_assoc <- disease2gene(disease = "C0005586",
                          database = "ALL", verbose = TRUE)
all_assoc <- extract(all_assoc)

results$Curated[results$Gene %in% curated$gene_symbol] <- "Yes"
results$All[results$Gene %in% all_assoc$gene_symbol] <- "Yes"

fwrite(results, paste0(args$output_results_folder, "/Formatted_Results.csv"))
