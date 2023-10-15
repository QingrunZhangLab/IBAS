# Loading necessary packages
pacman::p_load(data.table, uwot, Rtsne, ggplot2, argparse)

# Function to print log messages
print_log <- function(message){
  print(paste0("[", Sys.time(), "] ", message))
}

# Arguments for the script
parser <- argparse::ArgumentParser()

parser$add_argument("--expr_data", default = "../Brain_Frontal_Cortex_BA9.v8.residuals.matrix.csv", help = "Expression data file")
parser$add_argument("--pathway_data", default = "./GeneSets.csv", help = "Pathway data file (1st column pathway names, 2nd column genes in pathway)")
parser$add_argument("--gene_data", default = "../../IBAS-Stability-Internal/gencode.v26.GRCh38.genes.gtf", help = "Gene data file (gtf)")
parser$add_argument("--seed", default = 1234, help = "Seed for reproducibility")
parser$add_argument("--output_dir", default = "./Brain_Frontal_Cortex_01_DimReduc", help = "Output directory")

print_log("Parsing arguments")
args <- parser$parse_args()

print_log("Arguments parsed:")
print_log(paste0("Expression data: ", args$expr_data))
print_log(paste0("Pathway data: ", args$pathway_data))
print_log(paste0("Gene data: ", args$gene_data))
print_log(paste0("Seed: ", args$seed))
print_log(paste0("Output directory: ", args$output_dir))

print_log("Setting seed")
set.seed(args$seed)

print_log("Checking if output directory exists and creating if not")
if(!dir.exists(args$output_dir)){
  dir.create(args$output_dir, recursive = TRUE)
}

print_log("Checking if expression data exists")
if(!file.exists(args$expr_data)){
  stop("Expression data file does not exist")
}
print_log("Expression data exists")

print_log("Checking if pathway data exists")
if(!file.exists(args$pathway_data)){
  stop("Pathway data file does not exist")
}
print_log("Pathway data exists")

print_log("Checking if gene data exists")
if(!file.exists(args$gene_data)){
  stop("Gene data file does not exist")
}
print_log("Gene data exists")

print_log("Reading in gene data")
gencode <- as.data.table(stringr::str_split(fread(args$gene_data)$V9, ";", simplify = TRUE)[, c(1,4)])
colnames(gencode) <- c("geneID", "geneName")
gencode$geneID <- stringr::str_remove(gencode$geneID, "gene_id ")
gencode$geneName <- stringr::str_remove(gencode$geneName, "gene_name ")
gencode$geneID <- gsub('[\"]', '', gencode$geneID)
gencode$geneName <- gsub('[\"]', '', gencode$geneName)
gencode$geneName <- trimws(gencode$geneName, which = "both")
gencode <- unique(gencode)

print_log("Reading in pathway data")
pathway_data <- fread(args$pathway_data)
if(ncol(pathway_data) < 2){
  stop("Pathway data file must have at least 2 columns")
}
pathway_data <- pathway_data[, c(1,2)]
colnames(pathway_data) <- c("Pathway", "Genes")
# Detect string delimiter in Genes column
if(sum(stringr::str_count(pathway_data$Genes, ";")) == 0){
  stop("Genes column must be delimited by ;")
}
# Split genes column by delimiter and create named vector using names in GeneSet
pathway_genes <- stringr::str_split(pathway_data$Genes, ";")
names(pathway_genes) <- pathway_data$Pathway

# Creating unique list of genes
all_genes <- unique(unlist(pathway_genes))
# Identifying corresponding gene IDs
all_gene_ids <- gencode[gencode$geneName %in% all_genes, geneID]

print_log("Reading in expression data")
expr_data <- fread(args$expr_data)
print_log("Assuming matrix is of the form samples (rows) x genes (columns)")
print_log("Converting to matrix of gene expression")
expr_matrix <- t(as.matrix(expr_data[,-c(1)]))
colnames(expr_matrix) <- expr_data[[1]]

# print_log("Selecting only genes that are in the pathways provided")
# expr_matrix <- expr_matrix[,all_gene_ids]

for(i in seq(1, length(names(pathway_genes)))){
    curr_pathway <- names(pathway_genes)[i]
    print_log(paste0("Working on pathway ", curr_pathway, " (pathway ", i, " of ", length(names(pathway_genes)), ")"))
    curr_gene_ids <- gencode[gencode$geneName %in% pathway_genes[[curr_pathway]], geneID]
    curr_gene_ids <- curr_gene_ids[curr_gene_ids %in% colnames(expr_matrix)]

    curr_expr_matrix <- expr_matrix[,curr_gene_ids]
    print_log("Removing genes with zero variance")
    curr_expr_matrix <- curr_expr_matrix[, apply(curr_expr_matrix, 2, var) > 0]

    pathway_umap <- umap(curr_expr_matrix,
                        n_components = min(10, ncol(curr_expr_matrix)),
                        metric = "correlation",
                        n_epochs = 1000, verbose = FALSE)
    umap_df <- as.data.table(pathway_umap)
    umap_df <- cbind(data.table(ID = rownames(curr_expr_matrix)), umap_df)


    pathway_tsne <- Rtsne(curr_expr_matrix,
                         verbose = FALSE)
    tsne_df <- as.data.table(pathway_tsne$Y)
    tsne_df <- cbind(data.table(ID = rownames(curr_expr_matrix)), tsne_df)

    pathway_pca <- prcomp(curr_expr_matrix, center = TRUE, scale. = TRUE)
    pca_df <- as.data.table(pathway_pca$x[, seq(1,min(10, ncol(curr_expr_matrix)))])
    pca_df <- cbind(data.table(ID = rownames(curr_expr_matrix)), pca_df)

    print_log("Saving UMAP, tSNE and PCA embeddings")
    fwrite(umap_df, file.path(args$output_dir, paste0(curr_pathway, "_umap.tsv")), sep = "\t")
    fwrite(tsne_df, file.path(args$output_dir, paste0(curr_pathway, "_tsne.tsv")), sep = "\t")
    fwrite(pca_df, file.path(args$output_dir, paste0(curr_pathway, "_pca.tsv")), sep = "\t")
}
