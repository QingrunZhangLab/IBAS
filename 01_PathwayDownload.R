# Uses the pacman package to ensure all required packages are installed and loaded
if (!require("pacman")) install.packages("pacman")
pacman::p_load(data.table, 
               stringr,
               optparse,
               KEGGREST)

# Defining arguments to be parse from command line
option_list <- list(
  make_option(c("-o", "--output_file"), action = "store", type = "character", 
              default = "./PathwayData.csv",
              help = "Location of pathway database. Default location is %default. 
              If data is not found, it will be downloaded.")
)

# Loading parsed arguments
args <- parse_args(OptionParser(option_list = option_list))

pathways <- c()
pathways_downloaded = TRUE
if(file.exists(args$output_file)){
  pathways <- fread(args$output_file)
  if(all(colnames(pathways) != c("PathwayID", "PathwayClass", "PathwayName", "All_Genes")))pathways_downloaded = FALSE
  if(nrow(pathways) < 1) pathways_downloaded = FALSE
}else{
  pathways_downloaded = FALSE
}

if(!pathways_downloaded){
  print(paste0("Unable to locate pathway file at ", args$output_file, ". Downloading pathways..."))
  
  link <- keggLink("pathway", "hsa")
  unique_pathways <- unique(link)
  
  pathways <- data.table(PathwayID = character(),
                         PathwayClass = character(),
                         PathwayName = character(),
                         All_Genes = character())
  
  for(i in 1:length(unique_pathways)){
    query <- keggGet(unique_pathways[i])
    
    print(paste0("Working on ", query[[1]]$NAME, ". Pathway ", i, " of ", length(unique_pathways), "..."))
    
    genes <- query[[1]]$GENE
    if(is.null(genes))next
    
    gene_names <- genes[seq(2, length(genes), 2)]
    gene_names <- str_split(gene_names, ";", simplify = TRUE)[, 1]
    
    pathways <- rbindlist(list(pathways, data.table(PathwayID = query[[1]]$ENTRY,
                                                    PathwayClass = query[[1]]$CLASS,
                                                    PathwayName = query[[1]]$NAME,
                                                    All_Genes = paste0(gene_names, collapse = ";"))))
  }
  
  if(!dir.exists(dirname(args$output_file))){
    dir.create(dirname(args$output_file), recursive = TRUE)
  }
  
  fwrite(pathways, args$output_file)
}else{
  print(paste0("Pathways found and loaded from ", args$output_file))
}

print(paste0(nrow(pathways), " pathways can be found at ", args$output_file))