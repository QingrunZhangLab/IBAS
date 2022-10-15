source("IBAS_Functions.R")

downloadPathways()
pathwaysToGeneSets()

geneSets <- fread("./GeneSets.csv")

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

geneSetReducer(exprData = expr.data, sampleIDs = sample.ids, geneSets = geneSets, 
               geneLookupTable = gencode, method = "tSNE", numDims = 2, 
               outputFolder = "./ReducedGeneSets", seed = 1234)

geneSetReducer(exprData = expr.data, sampleIDs = sample.ids, geneSets = geneSets, 
               geneLookupTable = gencode, method = "UMAP", numDims = 2, 
               outputFolder = "./ReducedGeneSets", seed = 1234)


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

snpAvailability(gencode = mappings, geneSets = geneSets, cisBP = 1000000, outputFolder = "./SNPAvailability")