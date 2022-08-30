source("inTWAS_Functions.R")

pathways <- pathwayLoader(forceDownload = FALSE)
pathways$MedianEdges <- ceiling(pathways$MedianEdges)

# expr.data <- fread("../../Data/GTEx/GTEx_Analysis_v8_eQTL_expression_matrices/Whole_Blood.v8.normalized_expression.bed")
expr.data <- fread("../../Data/GTEx/Whole_Blood.v8.residuals.csv")
expr.data$start <- NULL
expr.data$end <- NULL
expr.data$`#chr` <- NULL
expr.data <- transpose(expr.data, keep.names = "SampleID", make.names = "gene_id")
sample.ids <- expr.data$SampleID
expr.data$SampleID <- NULL

gencode <- as.data.table(str_split(fread("../../Data/GTEx/gencode.v26.GRCh38.genes.gtf")$V9, ";", simplify = TRUE)[, c(1,4)])
colnames(gencode) <- c("geneID", "geneName")
gencode$geneID <- str_remove(gencode$geneID, "gene_id ")
gencode$geneName <- str_remove(gencode$geneName, "gene_name ")
gencode$geneID <- gsub('[\"]', '', gencode$geneID)
gencode$geneName <- gsub('[\"]', '', gencode$geneName)
gencode$geneName <- trimws(gencode$geneName, which = "both")
gencode <- unique(gencode)

pathwayReducer(exprData = expr.data, sampleIDs = sample.ids, pathwayData = pathways, 
               geneLookupTable = gencode, method = "tSNE", numDims = 2, 
               outputFolder = "./ReducedPathways", seed = 1234)

pathwayReducer(exprData = expr.data, sampleIDs = sample.ids, pathwayData = pathways, 
               geneLookupTable = gencode, method = "UMAP", numDims = 2, 
               outputFolder = "./ReducedPathways", seed = 1234)

gencode <- fread("../../Data/GTEx/gencode.v26.GRCh38.genes.gtf")
mappings <- as.data.table(str_split(gencode$V9, ";", simplify = TRUE)[, c(1, 4)])
colnames(mappings) <- c("gene_id", "gene_name")
mappings <- cbind(gencode, mappings)
mappings$gene_id <- str_replace(mappings$gene_id, "gene_id ", "")
mappings$gene_id <- str_replace_all(mappings$gene_id, '[\"]', "")
mappings$gene_name <- str_replace(mappings$gene_name, "gene_name ", "")
mappings$gene_name <- str_replace_all(mappings$gene_name, '[\"]', "")
mappings$gene_name <- str_trim(mappings$gene_name)

mappings <- mappings[mappings$V3 == "gene", ]
mappings <- mappings[mappings$gene_name %in% setdiff(unique(mappings$gene_name), unique(mappings$gene_name[duplicated(mappings$gene_name)]))]
mappings <- mappings[, c(1, 4, 5, 11)]
colnames(mappings) <- c("chromosome", "start", "end", "gene_name")
mappings$chromosome <- str_remove_all(mappings$chromosome, "chr")
mappings <- mappings[!grepl("\\D", mappings$chromosome), ]

snpAvailability(gencode = mappings, pathwayData = pathways, cisBP = 1000000, outputFolder = "./SNPAvailability")


results.dir <- "WTCCC"
all.files <- list.files(results.dir, full.names = TRUE)

all.results <- c()

for(f in all.files){
  file.name <- basename(f)
  file.name <- str_remove(file.name, ".csv")
  disease <- str_split(file.name, "_")[[1]][1]
  
  curr.data <- fread(f)
  curr.data <- data.table(dataset = results.dir,
                          disease = disease,
                          curr.data)
  
  all.results <- rbindlist(list(all.results, curr.data))
}

results.dir <- "NFBC"
all.files <- list.files(results.dir, full.names = TRUE)

for(f in all.files){
  file.name <- basename(f)
  file.name <- str_remove(file.name, ".csv")
  disease <- str_split(file.name, "_")[[1]][1]
  
  curr.data <- fread(f)
  curr.data <- data.table(dataset = results.dir,
                          disease = disease,
                          curr.data)
  
  all.results <- rbindlist(list(all.results, curr.data))
}

req.pathway.data <- pathways[, c("PathwayID", "PathwayClass", "PathwayName")]
all.results <- merge(all.results, req.pathway.data, by.x = "pathway", by.y = "PathwayID")

all.results <- all.results[, c("dataset", "disease", "method", "method_comp", "PathwayClass", "PathwayName", "pathway", "SetID", "P.value",
                       "N.Marker.All", "N.Marker.Test", "num_snp_sets")]
colnames(all.results) <- c("Dataset", "Disease", "Method", "Component", "PathwayClass", "PathwayName", "PathwayID", "Gene", "p.value", 
                       "N.Marker.All", "N.Marker.Test", "GenesInPathway")

adjusted.results <- c()

all.results$adj.P.value <- 1
for(dataset in unique(all.results$Dataset)){
  for(disease in unique(all.results$Disease)){
    for(method in unique(all.results$Method)){
      for(component in unique(all.results$Component)){
        for(pathway in unique(all.results$PathwayID)){
          curr.selection <- all.results[Dataset==dataset & Disease==disease & Method==method & Component==component & PathwayID == pathway]
          curr.selection$adj.P.value <- p.adjust(curr.selection$p.value, method = "bonf")
          adjusted.results <- rbindlist(list(adjusted.results, curr.selection))
        }
      }
    }
  }
}

adjusted.results$is.Sig <- ifelse(adjusted.results$adj.P.value < 0.05, "Yes", "No")

if (!require("devtools")) install.packages("devtools")
library(devtools)
if (!require("disgenet2r")) install_bitbucket("ibi_group/disgenet2r")
library(disgenet2r)

disgenet_api_key <- get_disgenet_api_key(
  email = "email", 
  password = "password")
Sys.setenv(DISGENET_API_KEY = disgenet_api_key)

CUIs <- fread("DiseaseCUIs.tsv", header = FALSE)
CUIs <- unique(CUIs)

adjusted.results$DisGeNETCurated <- "No"
for(disease in unique(CUIs$V1)){
  curated <- disease2gene(disease = CUIs$V2[CUIs$V1 == disease], 
                          database = "ALL", verbose = TRUE)
  curated <- extract(curated)
  adjusted.results$DisGeNETCurated[adjusted.results$Disease==disease & adjusted.results$Gene %in% curated$gene_symbol] <- "Yes"
}

fwrite(adjusted.results, "WTCCC_NFBC_SKAT_Results.csv")

plot.data <- adjusted.results[is.Sig == "Yes",]
table(plot.data$DisGeNETCurated, plot.data$Disease)
