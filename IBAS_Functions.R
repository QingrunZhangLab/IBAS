# Uses the pacman package to ensure all required packages are installed and loaded
if (!require("pacman")) install.packages("pacman")
pacman::p_load(data.table, 
               stringr,
               optparse,
               KEGGREST,
               KEGGgraph,
               Rtsne,
               umap)

downloadPathways <- function(downloadLocation = "./PathwayData.csv"){
  
  pathways <- c()

  link <- keggLink("pathway", "hsa")
  uniquePathways <- unique(link)
  
  pathways <- c()
  
  for(i in 1:length(uniquePathways)){
    query <- keggGet(uniquePathways[i])
    
    print(paste0("Working on ", query[[1]]$NAME, ". Pathway ", i, " of ", length(uniquePathways), "..."))
    
    tmpFile <- "./tmpfile.xml"
    file.create(tmpFile)
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
                                                    AllGenes = paste0(gene_names, collapse = ";"),
                                                    MedianEdges = medianEdges)))
  }
  
  if(!dir.exists(dirname(downloadLocation))){
    dir.create(dirname(downloadLocation), recursive = TRUE)
  }
  
  fwrite(pathways, downloadLocation)
  
  print(paste0(nrow(pathways), " pathways can be found at ", downloadLocation))
  
}

pathwaysToGeneSets <- function(pathwayLocation = "./PathwayData.csv", geneSetLocation = "./GeneSets.csv"){
  
  pathways <- fread(pathwayLocation)
  
  geneSets <- pathways[, c("PathwayID", "AllGenes", "MedianEdges")]
  colnames(geneSets) <- c("GeneSet", "Genes", "Edges")
  
  fwrite(geneSets, geneSetLocation)
  
}

# exprData <- expr.data
# sampleIDs <- sample.ids
# method = "tSNE"
# numDims <- 2
# pathwayData <- pathways
# geneLookupTable <- gencode
# outputFolder = "./ReducedPathways"
# seed = 2222
geneSetReducer <- function(exprData, sampleIDs, geneSets, 
                           geneLookupTable = NULL, method = "tSNE", numDims = 2,
                           edgeThreshold = 3,
                           outputFolder = "./ReducedPathways",
                           seed = 2222){
  
  geneSetData <- geneSets[Edges >= edgeThreshold]
  
  if(!dir.exists(outputFolder)) dir.create(outputFolder, recursive = TRUE)
  
  for(i in seq(1, nrow(geneSetData))){
    geneSetName <- geneSetData$GeneSet[i]
    geneSetGenes <- str_split(geneSetData$Genes[i], ";")[[1]]
    
    if(!is.null(geneLookupTable)){
      geneSetGenes <- geneLookupTable[[1]][geneLookupTable[[2]] %in% geneSetGenes]
    }
    
    genesToSelect <- intersect(geneSetGenes, colnames(exprData))
    
    if(length(genesToSelect)<1){
      print(paste0("Insufficient genes for geneset ", geneSetName, ". Skipping reduction."))
      next
    }
    
    geneSetExprData <- as.matrix(exprData[, ..genesToSelect])
    
    if(tolower(method) == "tsne"){
      set.seed(seed)
      geneSetReduced <- Rtsne(geneSetExprData, dims = numDims, perplexity = geneSetData$Edges[i])$Y
    }else if(tolower(method) == "umap"){
      geneSetReduced <- umap(geneSetExprData, n_neighbors = geneSetData$Edges[i], 
                             n_components = numDims, random_state = seed)$layout
    }else{
      stop("No valid method specified!")
    }
    
    geneSetReduced <- data.table(ID = sampleIDs, geneSetReduced)
    fwrite(geneSetReduced, paste0(outputFolder, "/", method, "_", geneSetName, ".tsv"), sep = "\t")
    
  }
  
}

# gencode <- mappings
# pathwayData <- pathways
# cisBP <- 1000000
# outputFolder = "./SNPAvailability"
snpAvailability <- function(gencode, geneSets, cisBP = 1000000, 
                            outputFolder = "./SNPAvailability"){
  
  if(!dir.exists(outputFolder)) dir.create(outputFolder, recursive = TRUE)
  
  for(i in 1:nrow(geneSets)){
    
    print(paste0("Currently working on ", geneSets$GeneSet[i], " (", i, " of ", nrow(geneSets), ")..."))
    if(file.exists(paste0(outputFolder, "/", geneSets$GeneSet[i], ".tsv"))){
      print(paste0("Geneset ", geneSets$GeneSet[i], " SNP file exists. Skipping..."))
      next
    }
    
    curr_genes <- str_split(geneSets$Genes[i], ";")[[1]]
    
    gene_details <- gencode[gencode$gene_name %in% curr_genes,]
    gene_details$start <- pmax(gene_details$start - cisBP, 0)
    gene_details$end <- gene_details$end + cisBP
    
    region_file <- gene_details[, 1:3]
    region_file$POS <- paste0(region_file$start, ":", region_file$end)
    region_file <- region_file[, c(1,4)]
    colnames(region_file) <- c("CHR", "LOC")
    
    print(paste0("Writing geneset ", geneSets$GeneSet[i]))
    fwrite(region_file, paste0(outputFolder, "/", geneSets$GeneSet[i], ".tsv"), sep = "\t")
  }
  
}