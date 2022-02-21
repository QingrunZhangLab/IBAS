# Uses the pacman package to ensure all required packages are installed and loaded
if (!require("pacman")) install.packages("pacman", repos = "https://cloud.r-project.org/")
pacman::p_load(data.table, 
               stringr,
               optparse,
               Rtsne)

# Defining arguments to be parse from command line
option_list <- list(
  make_option(c("-p", "--pathway_file"), action = "store", type = "character", 
              default = "./PathwayData.csv",
              help = "Location of pathway database. Default location is %default."),
  make_option(c("-e", "--expression_file"), action = "store", type = "character", 
              default = "../GTExExpressionData/Whole_Blood.v8.residuals.csv",
              help = "Location of expression data. Default location is %default."),
  make_option(c("-a", "--annotation_file"), action = "store", type = "character", 
              default = "../gencode.v26.GRCh38.genes.gtf",
              help = "Location of annotations (gtf). Default location is %default."),
  make_option(c("-s", "--samples_file"), action = "store", type = "character", 
              default = "../GTExSampleIDs/CommonIDs.txt",
              help = "Sample file containing sample IDs to keep. Default location is %default."),
  make_option(c("-d", "--dimensionality_reduction"), action = "store", type = "character", 
              default = "PCA,tSNE",
              help = "Generate PCs [PCA] and/or tSNEs [tSNEs]. Separate multiple inputs by ','s. Default is %default."),
  make_option(c("-P", "--pathway_availability_file"), action = "store", type = "character", 
              default = "./Reference_Pathway_Gene_Availability_Summary.csv",
              help = "Location to store pathway availability summary. Default location is %default."),
  make_option(c("-n", "--dim_num"), action = "store", type = "numeric", 
              default = "2",
              help = "Number of dimensions to reduce data to. Default is %default.")
)

# Loading parsed arguments
args <- parse_args(OptionParser(option_list = option_list))

pathways <- c()
pathways_downloaded = TRUE
if(file.exists(args$pathway_file)){
  pathways <- fread(args$pathway_file)
  if(all(colnames(pathways) !=c ("PathwayID", "PathwayClass", "PathwayName", "All_Genes")))pathways_downloaded = FALSE
  if(nrow(pathways) < 1) pathways_downloaded = FALSE
}else{
  pathways_downloaded = FALSE
}

if(!pathways_downloaded){
  stop(paste0("Unable to locate pathway file at ", args$pathway_file, 
               ". Please use Script 01 to download pathways."))
}else{
  print(paste0("Pathways found and loaded from ", args$pathway_file))
}

if(!file.exists(args$expression_file)){
  stop(paste0("Expression data: ", args$expression_file, " does not exist."))
}

full_data <- fread(args$expression_file)

if(!file.exists(args$samples_file)){
  stop(paste0("Samples data: ", args$samples_file, " does not exist."))
}

CommonIDs <- fread(args$samples_file, header = FALSE)

# Filtering and preparing expression data
sel_IDs <- colnames(full_data)[colnames(full_data) %in% CommonIDs$V1]
sel_IDs <- sel_IDs[order(sel_IDs, decreasing = TRUE)]
cols <- c("#chr", "start", "end", "gene_id", sel_IDs)
full_data <- full_data[, ..cols]

if(!file.exists(args$annotation_file)){
  stop(paste0("Annotation data: ", args$annotation_file, " does not exist."))
}

annotations <- fread(args$annotation_file, header = FALSE)

mappings <- as.data.table(str_split(annotations$V9, ";", simplify = TRUE)[, c(1, 4)])
colnames(mappings) <- c("gene_id", "gene_name")
mappings$gene_id <- str_replace(mappings$gene_id, "gene_id ", "")
mappings$gene_id <- str_replace_all(mappings$gene_id, '[\"]', "")
mappings$gene_name <- str_replace(mappings$gene_name, "gene_name ", "")
mappings$gene_name <- str_replace_all(mappings$gene_name, '[\"]', "")
mappings$gene_name <- str_trim(mappings$gene_name)

mappings <- unique(mappings)

full_data$gene_name <- mappings$gene_name[match(full_data$gene_id, mappings$gene_id)]

availability_summary <- data.table(PathwayID = character(),
                                   PathwayClass = character(),
                                   PathwayName = character(),
                                   TotalGenes = numeric(),
                                   AllGenes = character(),
                                   AvailableGeneNum = numeric(),
                                   AvailableGenes = character(),
                                   SelectedTranscriptNum = numeric(),
                                   SelectedTranscripts = character())

dim_methods <- str_split(args$dimensionality_reduction, ",")[[1]]

if(!dir.exists("ReducedPathways")){
  dir.create("ReducedPathways")
}

PC_dir <- "ReducedPathways/PCs"
PC_var_dir <- "ReducedPathways/PC_Vars"
tSNE_dir <- "ReducedPathways/tSNEs"

if("tSNE" %in% dim_methods){
  if(!dir.exists(tSNE_dir)){
    dir.create(tSNE_dir)
  }
}

if("PCA" %in% dim_methods){
  if(!dir.exists(PC_dir)){
    dir.create(PC_dir)
  }
  if(!dir.exists(PC_var_dir)){
    dir.create(PC_var_dir)
  }
}

for(i in 1:nrow(pathways)){
  
  print(paste0("Working on pathway ", i, " of ", nrow(pathways), "..."))
  pathway_genes <- str_split(pathways$All_Genes[i], ";")[[1]]
  
  available_genes <- pathway_genes[pathway_genes %in% full_data$gene_name]
  
  if(length(available_genes) < 2){
    print("Insufficient genes in pathway! Skipping...")
    next
  }
  
  selected_transcripts <- full_data$gene_id[full_data$gene_name %in% available_genes]
  
  availability_summary <- rbindlist(list(availability_summary,
                                         data.table(PathwayID = pathways$PathwayID[i],
                                                    PathwayClass = pathways$PathwayClass[i],
                                                    PathwayName = pathways$PathwayName[i],
                                                    TotalGenes = length(pathway_genes),
                                                    AllGenes = paste0(pathway_genes, collapse = ";"),
                                                    AvailableGeneNum = length(available_genes),
                                                    AvailableGenes = paste0(available_genes, collapse = ";"),
                                                    SelectedTranscriptNum = length(selected_transcripts),
                                                    SelectedTranscripts = paste0(selected_transcripts, collapse = ";"))))
  
  pathway_data <- full_data[full_data$gene_id %in% selected_transcripts, -c("#chr", "start", "end", "gene_name")]
  pathway_data <- as.data.frame(pathway_data)
  row.names(pathway_data) <- pathway_data$gene_id
  pathway_data$gene_id <- NULL
  pathway_data <- t(as.matrix(pathway_data))
  
  if("tSNE" %in% dim_methods){
    if(file.exists(paste0(tSNE_dir, "/", pathways$PathwayID[i], "_tSNEs.tsv"))){
      print(paste0("Pathway ", pathways$PathwayID[i], " tSNE file exists. Skipping..."))
    }else{
      print("Calculating tSNEs...")
      pathway_tsne <- Rtsne(pathway_data, dims = args$dim_num, num_threads = 2)
      
      pathway_tsne <- as.data.table(cbind(row.names(pathway_data), 
                                          pathway_tsne$Y[,1:args$dim_num]))
      colnames(pathway_tsne) <- c("ID", paste0("tSNE", 1:args$dim_num))
      
      fwrite(pathway_tsne, paste0(tSNE_dir, "/", pathways$PathwayID[i], "_tSNEs.tsv"), sep = "\t")
    }
  }
  
  if("PCA" %in% dim_methods){
    if(file.exists(paste0(PC_dir, "/", pathways$PathwayID[i], "_PCs.tsv"))){
      print(paste0("Pathway ", pathways$PathwayID[i], " PC file exists. Skipping..."))
    }else{
      print("Calculating PCAs...")
      pathway_pca <- prcomp(pathway_data, center = TRUE, scale = TRUE)
      
      pathway_pcs <- as.data.table(cbind(row.names(pathway_pca$x), pathway_pca$x[,1:min(args$dim_num, ncol(pathway_pca$x))]))
      colnames(pathway_pcs)[1] <- "ID"
      colnames(pathway_pcs)[2:ncol(pathway_pcs)] <- paste0("PC", 1:(ncol(pathway_pcs)-1))
      
      impor <- t(summary(pathway_pca)$importance[,1:min(args$dim_num, ncol(pathway_pca$x))])
      
      fwrite(pathway_pcs, paste0(PC_dir, "/", pathways$PathwayID[i], "_PCs.tsv"), sep = "\t")
      fwrite(impor, paste0(PC_var_dir, "/", pathways$PathwayID[i], "_Vars.tsv"), sep = "\t")
    }
  }
}

fwrite(availability_summary, args$pathway_availability_file)