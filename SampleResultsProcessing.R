
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

pathways <- fread("PathwayData.csv")
req.pathway.data <- pathways[, c("PathwayID", "PathwayClass", "PathwayName")]
all.results <- merge(all.results, req.pathway.data, by.x = "geneSet", by.y = "PathwayID")

all.results <- all.results[, c("dataset", "disease", "method", "method_comp", "PathwayClass", "PathwayName", "geneSet", "SetID", "P.value",
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
  email = "thalagalakossinnagep@ucalgary.ca", 
  password = "e8wt3h35$")
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
