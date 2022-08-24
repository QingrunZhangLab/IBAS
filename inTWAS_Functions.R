# Uses the pacman package to ensure all required packages are installed and loaded
if (!require("pacman")) install.packages("pacman")
pacman::p_load(data.table, 
               stringr,
               optparse,
               KEGGREST,
               KEGGgraph)

downloadPathways <- function(downloadLocation = "./PathwayData.csv"){
  
  pathways <- c()

  link <- keggLink("pathway", "hsa")
  uniquePathways <- unique(link)
  
  pathways <- c()
  
  for(i in 1:length(uniquePathways)){
    query <- keggGet(uniquePathways[i])
    
    print(paste0("Working on ", query[[1]]$NAME, ". Pathway ", i, " of ", length(uniquePathways), "..."))
    
    tmpFile <- tempfile(fileext = ".xml")
    tryCatch(retrieveKGML(uniquePathways[i], organism = "hsa", destfile = tmpFile, method = "libcurl", quiet = TRUE), error = function(e) print("Error?"))
    mapkG <- parseKGML2Graph(tmpFile, expandGenes = TRUE, genesOnly = TRUE)
    medianEdges <- median(lengths(edges(mapkG)))
    
    genes <- query[[1]]$GENE
    if(is.null(genes))next
    
    gene_names <- genes[seq(2, length(genes), 2)]
    gene_names <- str_split(gene_names, ";", simplify = TRUE)[, 1]
    
    pathways <- rbindlist(list(pathways, data.table(PathwayID = query[[1]]$ENTRY,
                                                    PathwayClass = query[[1]]$CLASS,
                                                    PathwayName = query[[1]]$NAME,
                                                    All_Genes = paste0(gene_names, collapse = ";"),
                                                    MedianEdges = medianEdges)))
  }
  
  if(!dir.exists(dirname(downloadLocation))){
    dir.create(dirname(downloadLocation), recursive = TRUE)
  }
  
  fwrite(pathways, downloadLocation)
  
  print(paste0(nrow(pathways), " pathways can be found at ", downloadLocation))
  
}

pathwayLoader <- function(cacheLocation = "./PathwayData.csv", forceDownload = FALSE){
  
  if(forceDownload | !file.exists(cacheLocation)){
    downloadPathways(cacheLocation)
  }
  
  return(fread(cacheLocation))
  
}

pathwayReducer <- function(exprData, method = "tSNE", numDims = 2, pathwayData){
  
  
  
}