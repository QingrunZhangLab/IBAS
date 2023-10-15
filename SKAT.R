# Loading necessary packages
pacman::p_load(data.table,
                stringr,
                argparse,
                GenomicRanges,
                metap,
                SKAT)

# Function to print log messages
print_log <- function(message){
  print(paste0("[", Sys.time(), "] ", message))
}

# Arguments for the script
parser <- argparse::ArgumentParser()

parser$add_argument("--pathway_data", default = "./GeneSets.csv", help = "Pathway data file (1st column pathway names, 2nd column genes in pathway)")
parser$add_argument("--gene_data", default = "../../IBAS-Stability-Internal/gencode.v26.GRCh38.genes.gtf", help = "Gene data file (gtf)")
parser$add_argument("--genotype_data", default = "../SampleData/ADAPTmap_genotypeTOP_20160222_full", help = "Prefix to bed, bim and fam files to test association. Default location is %default.")
parser$add_argument("--results_prefix", default = "T1D_SKAT", help = "Prefix of location to store final results. A CSV file will be generated. Default location is %default.")
parser$add_argument("--alpha_val", default = "0.05", help = "Level of significance at which to consider SNPs for SKAT. Default value is %default.")
parser$add_argument("--cisBP", default = "1000000", help = "Number of cis-bps to consider. Default value is %default.")
parser$add_argument("--phenotype_type", default = "D", help = "An indicator of the outcome type. 'C' for the continuous outcome and 'D' for the dichotomous outcome. Default value is %default.")
parser$add_argument("--weight_type", default = "b", help = "An indicator of the type of weights to use. 'b' uses the betas of the regression while 'p' uses the p-values. Default value is %default.")
parser$add_argument("--temp_folder", default = "TEMP", help = "Folder in which to temporarily store SKAT files. Default location is %default.")
parser$add_argument("--weight_directory", default = "GeneSetWeights/hsa00010_pca", help = "Folder containing EMMAX results with weights.")
parser$add_argument("--output_dir", default = "./03_SKATResults", help = "Output directory")

print_log("Parsing arguments")
args <- parser$parse_args()

print_log("Arguments parsed:")
print_log(paste0("Pathway data: ", args$pathway_data))
print_log(paste0("Gene data: ", args$gene_data))
print_log(paste0("Genotype data: ", args$genotype_data))
print_log(paste0("Results prefix: ", args$results_prefix))
print_log(paste0("Alpha value: ", args$alpha_val))
print_log(paste0("cisBP: ", args$cisBP))
print_log(paste0("Phenotype type: ", args$phenotype_type))
print_log(paste0("Weight type: ", args$weight_type))
print_log(paste0("Temp folder: ", args$temp_folder))
print_log(paste0("Weight directory: ", args$weight_directory))
print_log(paste0("Output directory: ", args$output_dir))

if(!file.exists(paste0(args$genotype_data, ".bed")) | !file.exists(paste0(args$genotype_data, ".bim")) | !file.exists(paste0(args$genotype_data, ".fam"))){
  stop(paste0("Genotype data: ", args$genotype_data, " does not exist."))
}

if(!file.exists(args$gene_data)){
  stop(paste0("Annotation data: ", args$gene_data, " does not exist."))
}

if(!file.exists(args$pathway_data)){
  stop(paste0("Gene Set data: ", args$pathway_data, " does not exist."))
}

if(!dir.exists(args$temp_folder)){
  dir.create(args$temp_folder, recursive = TRUE)
}

if(!dir.exists(args$output_dir)){
  dir.create(args$output_dir, recursive = TRUE)
}

if(!dir.exists(args$weight_directory)){
  stop(paste0("Weight directory: ", args$weight_directory, " does not exist."))
}

print_log("Reading in pathway data")
pathway_data <- fread(args$pathway_data)

print_log("Reading in gene data")
annotations <- fread(args$gene_data, header = FALSE)

print_log("Creating mappings")
mappings <- as.data.table(str_split(annotations$V9, ";", simplify = TRUE)[, c(1, 4)])
colnames(mappings) <- c("gene_id", "gene_name")
mappings <- cbind(annotations, mappings)
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

File.Bed <- paste0(args$genotype_data, ".bed")
File.Bim <- paste0(args$genotype_data, ".bim")
File.Fam <- paste0(args$genotype_data, ".fam")

print_log("Reading in bim data")
bim_data <- fread(File.Bim)[, c(1,4,2)]
colnames(bim_data) <- c("chr", "loc", "rs_id")
if(0 %in% bim_data$chr) bim_data$chr <- bim_data$chr + 1

weights_to_process <- list.files(args$weight_directory, pattern = ".top", full.names = TRUE)

# Read in all weight files with a new column indicating the number before ".top"
print_log("Reading in weight files...")
all_weights <- rbindlist(lapply(weights_to_process, function(x){
  tmp <- fread(x)
  tmp$component <- as.numeric(gsub(".*_PC([0-9]+).*", "\\1", basename(x)))
  tmp
}))

all_weights <- all_weights[, c("#chr", "location", "component", "coefficient", "pvalue")]

# Collapse weights by averaging coefficients and using sumlog from the metap package for the pvalues
print_log("Collapsing weights...")
collapsed_weights <- all_weights[, list(coefficient = mean(coefficient), pvalue = sumlog(pvalue)$p), by = c("#chr", "location")]
collapsed_weights <- collapsed_weights[complete.cases(collapsed_weights),]

all_results <- c()

curr_pathway <- str_split(basename(args$weight_directory), "_")[[1]][1]
curr_method <- str_split(basename(args$weight_directory), "_")[[1]][2]
curr_genes <- str_split(pathway_data[[2]][pathway_data[[1]] == curr_pathway], ";")[[1]]

gene_details <- mappings[mappings$gene_name %in% curr_genes,]

gene_details$chromosome <- as.integer(gene_details$chromosome)
gene_details$start <- pmax(gene_details$start - as.numeric(args$cisBP), 0) 
gene_details$end <- gene_details$end + as.numeric(args$cisBP)

for(curr_weight_file in weights_to_process){
    curr_component <-  as.numeric(gsub(".*_.*?([0-9]+)\\..*", "\\1", basename(curr_weight_file)))
    print_log(paste0("Currently working on component ", curr_component, " of ", curr_pathway, "..."))

    print_log("Filtering weights...")
    curr_weights <- all_weights[component == curr_component & pvalue <= as.numeric(args$alpha_val),]

    if(nrow(curr_weights) == 0){
      print_log(paste0("No significant weights for component ", curr_component, " of ", curr_pathway, ". Skipping..."))
      next
    }

    print_log("Creating weight data based on pvalue or coefficient...")
    if(args$weight_type == "p"){
        curr_weights <- curr_weights[, c("#chr", "location", "pvalue")]
        colnames(curr_weights) <- c("chr", "loc", "weight")
        curr_weights$weight[curr_weights$weight == 0] <- .Machine$double.eps
        curr_weights$weight <- -log10(curr_weights$weight)
    }else{
        curr_weights <- curr_weights[, c("#chr", "location", "coefficient")]
        colnames(curr_weights) <- c("chr", "loc", "weight")
        curr_weights$weight <- abs(curr_weights$weight)
    }

    print_log("Merging weights with bim data...")
    curr_results <- merge(curr_weights, bim_data, by = c("chr", "loc"))
  
    if(nrow(curr_results)==0){
        print_log(paste0("No SNPs found for component ", curr_component, " of ", curr_pathway, ". Skipping..."))
        next
    }

    curr_SNP_Weights <- curr_results[, c("rs_id", "weight")]
    if(!dir.exists(paste0(args$temp_folder, "/", curr_pathway, "_", curr_method, "_", curr_component))){
        dir.create(paste0(args$temp_folder, "/", curr_pathway, "_", curr_method, "_", curr_component), recursive = TRUE)
    }
    fwrite(curr_SNP_Weights, paste0(args$temp_folder, "/", curr_pathway, "_", curr_method, "_", curr_component, "/Snp.Weights.txt"), col.names = FALSE, sep = "\t")

    print_log("Finding SNPs in genes...")
    isnps <- with(curr_results, GRanges(curr_results$chr, IRanges(loc, width=1, names=rs_id), "*"))
    igenes <- with(gene_details, GRanges(gene_details$chromosome, IRanges(start, end, names=gene_name), "*"))
    olaps <- findOverlaps(isnps, igenes)
  
    output <- cbind(curr_results[queryHits(olaps),], gene_details[subjectHits(olaps),])
    output$distance <- abs((output$end+output$start)/2-output$loc)
    output <- output[order(output$rs_id, output$distance ), ]
    output <- output[ !duplicated(output$rs_id), ]
  
    curr_SNP_set <- output[, c("gene_name", "rs_id")]
    fwrite(curr_SNP_set, paste0(args$temp_folder, "/", curr_pathway, "_", curr_method, "_", curr_component, "/Snp.Sets.txt"), col.names = FALSE, sep = "\t")
  
    print(paste0("Running SKAT for component ", curr_component, " of ", curr_pathway, " (", curr_method, ")..."))
    File.SetID <- paste0(args$temp_folder, "/", curr_pathway, "_", curr_method, "_", curr_component, "/Snp.Sets.txt")
    File.SSD <- paste0(args$temp_folder, "/", curr_pathway, "_", curr_method, "_", curr_component, "/temp.SSD")
    File.Info <-paste0(args$temp_folder, "/", curr_pathway, "_", curr_method, "_", curr_component, "/temp.SSD.Info")
  
    SnpSets <- fread(paste0(args$temp_folder, "/", curr_pathway, "_", curr_method, "_", curr_component, "/Snp.Sets.txt"), header = FALSE)
    obj.SNPWeight <- Read_SNP_WeightFile(paste0(args$temp_folder, "/", curr_pathway, "_", curr_method, "_", curr_component, "/Snp.Weights.txt"))
  
    Generate_SSD_SetID(File.Bed, File.Bim, File.Fam, File.SetID, File.SSD, File.Info)
    FAM <- Read_Plink_FAM(File.Fam, Is.binary=(args$phenotype_type=="D"))
    y <- FAM$Phenotype
  
    SSD.INFO <- Open_SSD(File.SSD, File.Info)
    obj <- SKAT_Null_Model(y ~ 1, out_type=args$phenotype_type)
  
    out <- SKAT.SSD.All(SSD.INFO, obj, obj.SNPWeight = obj.SNPWeight)

    print_log("Formatting results...")
    curr_outputs <- as.data.table(out$results)
  
    curr_outputs$pathway <- curr_pathway
    curr_outputs$method <- curr_method
    curr_outputs$num_snp_sets <- length(unique(SnpSets$V1))
    curr_outputs$component <- curr_component
  
    all_results[[paste0(curr_pathway, "_", curr_method, "_", curr_component)]] <- curr_outputs
}

all_results <- rbindlist(all_results)

# Running meta-analysis
print_log("Running meta-analysis...")

curr_weights <- collapsed_weights[pvalue <= as.numeric(args$alpha_val),]

if(nrow(curr_weights) == 0){
    print_log(paste0("No significant weights for meta-analysis of ", curr_pathway, ". Skipping..."))
    stop()
}

if(args$weight_type == "p"){
    curr_weights <- curr_weights[, c("#chr", "location", "pvalue")]
    colnames(curr_weights) <- c("chr", "loc", "weight")
    curr_weights$weight[curr_weights$weight == 0] <- .Machine$double.eps
    curr_weights$weight <- -log10(curr_weights$weight)
}else{
    curr_weights <- curr_weights[, c("#chr", "location", "coefficient")]
    colnames(curr_weights) <- c("chr", "loc", "weight")
    curr_weights$weight <- abs(curr_weights$weight)
}

curr_results <- merge(curr_weights, bim_data, by = c("chr", "loc"))

if(nrow(curr_results)==0){
    print_log(paste0("No SNPs found for meta-analysis of ", curr_pathway, ". Skipping..."))
    stop()
}

curr_SNP_Weights <- curr_results[, c("rs_id", "weight")]
if(!dir.exists(paste0(args$temp_folder, "/", curr_pathway, "_", curr_method, "_meta"))){
    dir.create(paste0(args$temp_folder, "/", curr_pathway, "_", curr_method, "_meta"), recursive = TRUE)
}
fwrite(curr_SNP_Weights, paste0(args$temp_folder, "/", curr_pathway, "_", curr_method, "_meta", "/Snp.Weights.txt"), col.names = FALSE, sep = "\t")

isnps <- with(curr_results, GRanges(curr_results$chr, IRanges(loc, width=1, names=rs_id), "*"))
igenes <- with(gene_details, GRanges(gene_details$chromosome, IRanges(start, end, names=gene_name), "*"))
olaps <- findOverlaps(isnps, igenes)

output <- cbind(curr_results[queryHits(olaps),], gene_details[subjectHits(olaps),])
output$distance <- abs((output$end+output$start)/2-output$loc)
output <- output[order(output$rs_id, output$distance ), ]
output <- output[ !duplicated(output$rs_id), ]

curr_SNP_set <- output[, c("gene_name", "rs_id")]
fwrite(curr_SNP_set, paste0(args$temp_folder, "/", curr_pathway, "_", curr_method, "_meta", "/Snp.Sets.txt"), col.names = FALSE, sep = "\t")

print(paste0("Running SKAT for meta-analysis of ", curr_pathway, " (", curr_method, ")..."))
File.SetID <- paste0(args$temp_folder, "/", curr_pathway, "_", curr_method, "_meta", "/Snp.Sets.txt")
File.SSD <- paste0(args$temp_folder, "/", curr_pathway, "_", curr_method, "_meta", "/temp.SSD")
File.Info <- paste0(args$temp_folder, "/", curr_pathway, "_", curr_method, "_meta", "/temp.SSD.Info")

SnpSets <- fread(paste0(args$temp_folder, "/", curr_pathway, "_", curr_method, "_meta", "/Snp.Sets.txt"), header = FALSE)
obj.SNPWeight <- Read_SNP_WeightFile(paste0(args$temp_folder, "/", curr_pathway, "_", curr_method, "_meta", "/Snp.Weights.txt"))

Generate_SSD_SetID(File.Bed, File.Bim, File.Fam, File.SetID, File.SSD, File.Info)
FAM <- Read_Plink_FAM(File.Fam, Is.binary=(args$phenotype_type=="D"))
y <- FAM$Phenotype

SSD.INFO <- Open_SSD(File.SSD, File.Info)
obj <- SKAT_Null_Model(y ~ 1, out_type=args$phenotype_type)

out <- SKAT.SSD.All(SSD.INFO, obj, obj.SNPWeight = obj.SNPWeight)

print_log("Formatting results...")
curr_outputs <- as.data.table(out$results)

curr_outputs$pathway <- curr_pathway
curr_outputs$method <- curr_method
curr_outputs$num_snp_sets <- length(unique(SnpSets$V1))
curr_outputs$component <- "meta"

all_results <- rbindlist(list(all_results, curr_outputs))

print_log("Writing results to file...")
if(!dir.exists(paste0(args$output_dir, "/", args$results_prefix))){
    dir.create(paste0(args$output_dir, "/", args$results_prefix), recursive = TRUE)
}
fwrite(all_results, paste0(args$output_dir, "/", args$results_prefix, "/", curr_pathway, "_", curr_method, ".csv"), sep = "\t")
