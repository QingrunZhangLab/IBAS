# Uses the pacman package to ensure all required packages are installed and loaded
if (!require("pacman")) install.packages("pacman")
pacman::p_load(data.table, 
               stringr,
               optparse,
               GenomicRanges,
               SKAT)

# Defining arguments to be parse from command line
option_list <- list(
  make_option(c("-P", "--pathway_availability_file"), action = "store", type = "character", 
              default = "./Reference_Pathway_Gene_Availability_Summary.csv",
              help = "Location to load pathway availability summary. Default location is %default."),
  make_option(c("-g", "--genotype_data"), action = "store", type = "character", 
              default = "../../WTCCC_Split/BD/BD.final",
              help = "Prefix to bed, bim and fam files to test association. Default location is %default."),
  make_option(c("-A", "--annotation_file"), action = "store", type = "character", 
              default = "../gencode.v26.GRCh38.genes.gtf",
              help = "Location of annotations. Default location is %default."),
  make_option(c("-r", "--results_file"), action = "store", type = "character", 
              default = "Pathway_SKAT_Results",
              help = "Prefix of location to store final results. A CSV file will be generated. 
              Default location is %default."),
  make_option(c("-a", "--alpha_val"), action = "store", type = "numeric", 
              default = "0.05",
              help = "Level of significance at which to consider SNPs for SKAT. Default value is %default."),
  make_option(c("-p", "--phenotype_type"), action = "store", type = "character", 
              default = "D",
              help = "An indicator of the outcome type. 'C' for the continuous outcome 
              and 'D' for the dichotomous outcome. Default value is %default."),
  make_option(c("-t", "--temp_folder"), action = "store", type = "character", 
              default = "TEMP",
              help = "Folder in which to temporarily store SKAT files. Default location is %default.")
)

# Loading parsed arguments
args <- parse_args(OptionParser(option_list = option_list))

if(!file.exists(paste0(args$genotype_data, ".bed")) | !file.exists(paste0(args$genotype_data, ".bim")) | !file.exists(paste0(args$genotype_data, ".fam"))){
  stop(paste0("Genotype data: ", args$genotype_file, " does not exist."))
}

if(!file.exists(args$annotation_file)){
  stop(paste0("Annotation data: ", args$annotation_file, " does not exist."))
}

if(!file.exists(args$pathway_availability_file)){
  stop(paste0("Pathway Availability data: ", args$pathway_availability_file, " does not exist."))
}


if(!dir.exists(args$temp_folder)){
  dir.create(args$temp_folder)
}

pathways <- fread(args$pathway_availability_file)

annotations <- fread(args$annotation_file, header = FALSE)

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

bim_data <- fread(File.Bim)[, c(1,4,2)]
colnames(bim_data) <- c("chr", "loc", "rs_id")

PC_dirs <- c()
tSNE_dirs <- c()

if(dir.exists("Reference_PathwayWeights/PCs")){
  PC_dirs <- list.dirs("Reference_PathwayWeights/PCs/", recursive = FALSE, full.names = FALSE)
}

if(dir.exists("Reference_PathwayWeights/tSNEs")){
  tSNE_dirs <- list.dirs("Reference_PathwayWeights/tSNEs/", recursive = FALSE, full.names = FALSE)
}

full_list <- union(PC_dirs, tSNE_dirs)

all_results <- c()

for(d in full_list){
  print(paste0("Currently working on ", d, "..."))
  
  if(!dir.exists(paste0(args$temp_folder, "/", d))){
    dir.create(paste0(args$temp_folder, "/", d), recursive = TRUE)
  }
  
  curr_genes <- str_split(pathways$AllGenes[pathways$PathwayID == d], ";")[[1]]
  
  gene_details <- mappings[mappings$gene_name %in% curr_genes,]
  
  gene_details$chromosome <- as.integer(gene_details$chromosome)
  gene_details$start <- pmax(gene_details$start - 1000000, 0)
  gene_details$end <- gene_details$end + 1000000
  
  if(d %in% PC_dirs){
    pc_list <- list.files(paste0("Reference_PathwayWeights/PCs/", d, "/"))
    pc_list <- str_replace(pc_list, "EMMAX.*_", "")
    pc_list <- str_replace(pc_list, ".top", "")
    nums <- str_replace(pc_list, "PC", "")
    for(pc in pc_list){
      if(!dir.exists(paste0(args$temp_folder, "/", d, "/", pc))){
        dir.create(paste0(args$temp_folder, "/", d, "/", pc))
      }
    }
    for(n in nums){
      curr_results <- fread(paste0("Reference_PathwayWeights/PCs/", d, "/EMMAX.", (as.numeric(n)-1), "_PC", n, ".top"), skip = "#chr")
      curr_results <- curr_results[curr_results$pvalue <= args$alpha_val,]
      if(nrow(curr_results)==0){
        print("Insufficient results for PC", n, " of ", d, " to generate SKAT files...")
        next
      }
      
      curr_results <- curr_results[, 1:3]
      colnames(curr_results) <- c("chr", "loc", "weight")
      min_v <- min(curr_results$weight[curr_results$weight != 0])
      curr_results$weight[curr_results$weight == 0] <- min_v/2
      curr_results$weight <- -log10(curr_results$weight)
      curr_results <- merge(curr_results, bim_data, by = c("chr", "loc"))
      if(nrow(curr_results)==0){
        next
      }
      curr_SNP_Weights <- curr_results[, c("rs_id", "weight")]
      fwrite(curr_SNP_Weights, paste0(args$temp_folder, "/", d, "/PC", n, "/Snp.Weights.txt"), col.names = FALSE, sep = "\t")
      
      isnps <- with(curr_results, GRanges(curr_results$chr, IRanges(loc, width=1, names=rs_id), "*"))
      igenes <- with(gene_details, GRanges(gene_details$chromosome, IRanges(start, end, names=gene_name), "*"))
      olaps <- findOverlaps(isnps, igenes)
      
      output <- cbind(curr_results[queryHits(olaps),], gene_details[subjectHits(olaps),])
      output$distance <- abs((output$end+output$start)/2-output$loc)
      output <- output[order(output$rs_id, output$distance ), ]
      output <- output[ !duplicated(output$rs_id), ]
      
      curr_SNP_set <- output[, c("gene_name", "rs_id")]
      fwrite(curr_SNP_set, paste0(args$temp_folder, "/", d, "/PC", n, "/Snp.Sets.txt"), col.names = FALSE, sep = "\t")
      
      print(paste0("Running SKAT for PC ", n, " of ", d,"..."))
      File.SetID <- paste0(args$temp_folder, "/", d, "/PC", n, "/Snp.Sets.txt")
      File.SSD <- paste0(args$temp_folder, "/", d, "/PC", n, "/temp.SSD")
      File.Info <- paste0(args$temp_folder, "/", d, "/PC", n, "/temp.SSD.Info")
      
      SnpSets <- fread(paste0(args$temp_folder, "/", d, "/PC", n, "/Snp.Sets.txt"), header = FALSE)
      obj.SNPWeight <- Read_SNP_WeightFile(paste0(args$temp_folder, "/", d, "/PC", n, "/Snp.Weights.txt"))
      
      Generate_SSD_SetID(File.Bed, File.Bim, File.Fam, File.SetID, File.SSD, File.Info)
      FAM <- Read_Plink_FAM(File.Fam, Is.binary=TRUE)
      y <- FAM$Phenotype
      
      SSD.INFO <- Open_SSD(File.SSD, File.Info)
      obj <- SKAT_Null_Model(y ~ 1, out_type=args$phenotype_type)
      
      out <- SKAT.SSD.All(SSD.INFO, obj, obj.SNPWeight = obj.SNPWeight)
      
      curr_results <- as.data.table(out$results)
      
      curr_results$pathway <- d
      curr_results$num_snp_sets <- length(unique(SnpSets$V1))
      curr_results$method_dim <- paste0("PC", n)
      
      all_results <- rbindlist(list(all_results, curr_results))
    }
  }
  
  if(d %in% tSNE_dirs){
    tSNE_list <- list.files(paste0("Reference_PathwayWeights/tSNEs/", d, "/"))
    tSNE_list <- str_replace(tSNE_list, "EMMAX.*_", "")
    tSNE_list <- str_replace(tSNE_list, ".top", "")
    nums <- str_replace(tSNE_list, "tSNE", "")
    for(tSNE in tSNE_list){
      if(!dir.exists(paste0(args$temp_folder, "/", d, "/", tSNE))){
        dir.create(paste0(args$temp_folder, "/", d, "/", tSNE))
      }
    }
    for(n in nums){
      curr_results <- fread(paste0("Reference_PathwayWeights/tSNEs/", d, "/EMMAX.", (as.numeric(n)-1), "_tSNE", n, ".top"), skip = "#chr")
      curr_results <- curr_results[curr_results$pvalue <= args$alpha_val,]
      if(nrow(curr_results) == 0){
        print("Insufficient results for tSNE", n, " of ", d, " to generate SKAT files...")
        next
      }
      
      curr_results <- curr_results[, 1:3]
      colnames(curr_results) <- c("chr", "loc", "weight")
      min_v <- min(curr_results$weight[curr_results$weight != 0])
      curr_results$weight[curr_results$weight == 0] <- min_v/2
      curr_results$weight <- -log10(curr_results$weight)
      curr_results <- merge(curr_results, bim_data, by = c("chr", "loc"))
      if(nrow(curr_results) == 0){
        next
      }
      curr_SNP_Weights <- curr_results[, c("rs_id", "weight")]
      fwrite(curr_SNP_Weights, paste0(args$temp_folder, "/", d, "/tSNE", n, "/Snp.Weights.txt"), col.names = FALSE, sep = "\t")
      
      isnps <- with(curr_results, GRanges(curr_results$chr, IRanges(loc, width=1, names=rs_id), "*"))
      igenes <- with(gene_details, GRanges(gene_details$chromosome, IRanges(start, end, names=gene_name), "*"))
      olaps <- findOverlaps(isnps, igenes)
      
      output <- cbind(curr_results[queryHits(olaps),], gene_details[subjectHits(olaps),])
      output$distance <- abs((output$end+output$start)/2-output$loc)
      output <- output[order(output$rs_id, output$distance ), ]
      output <- output[ !duplicated(output$rs_id), ]
      
      curr_SNP_set <- output[, c("gene_name", "rs_id")]
      fwrite(curr_SNP_set, paste0(args$temp_folder, "/", d, "/tSNE", n, "/Snp.Sets.txt"), col.names = FALSE, sep = "\t")
      
      print(paste0("Running SKAT for tSNE ", n, " of ", d,"..."))
      File.SetID <- paste0(args$temp_folder, "/",d, "/tSNE", n, "/Snp.Sets.txt")
      File.SSD <- paste0(args$temp_folder, "/", d, "/tSNE", n, "/temp.SSD")
      File.Info <- paste0(args$temp_folder, "/", d, "/tSNE", n, "/temp.SSD.Info")
      
      SnpSets <- fread(paste0(args$temp_folder, "/", d, "/tSNE", n, "/Snp.Sets.txt"), header = FALSE)
      obj.SNPWeight <- Read_SNP_WeightFile(paste0(args$temp_folder, "/", d, "/tSNE", n, "/Snp.Weights.txt"))
      
      Generate_SSD_SetID(File.Bed, File.Bim, File.Fam, File.SetID, File.SSD, File.Info)
      FAM <- Read_Plink_FAM(File.Fam, Is.binary = TRUE)
      y <- FAM$Phenotype
      
      SSD.INFO <- Open_SSD(File.SSD, File.Info)
      obj <- SKAT_Null_Model(y ~ 1, out_type = args$phenotype_type)
      
      out <- SKAT.SSD.All(SSD.INFO, obj, obj.SNPWeight = obj.SNPWeight)
      
      curr_results <- as.data.table(out$results)
      
      curr_results$pathway <- d
      curr_results$num_snp_sets <- length(unique(SnpSets$V1))
      curr_results$method_dim <- paste0("tSNE", n)
      
      all_results <- rbindlist(list(all_results, curr_results))
    }
  }
}

fwrite(all_results, paste0(args$results_file, ".csv"))