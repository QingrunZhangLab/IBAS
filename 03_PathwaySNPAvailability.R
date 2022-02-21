# Uses the pacman package to ensure all required packages are installed and loaded
if (!require("pacman")) install.packages("pacman", repos = "https://cloud.r-project.org/")
pacman::p_load(data.table,
               stringr,
               optparse)

# Defining arguments to be parse from command line
option_list <- list(
  make_option(c("-a", "--annotation_file"), action = "store", type = "character", 
              default = "../gencode.v26.GRCh38.genes.gtf",
              help = "Location of annotations. Default location is %default."),
  make_option(c("-P", "--pathway_availability_file"), action = "store", type = "character", 
              default = "./Reference_Pathway_Gene_Availability_Summary.csv",
              help = "Location to load pathway availability summary. Default location is %default."),
  make_option(c("-g", "--genotype_file"), action = "store", type = "character", 
              default = "../GTExGenotypeData/GTExFinal.num.csv",
              help = "Location of genotype data. Default location is %default.")
)

# Loading parsed arguments
args <- parse_args(OptionParser(option_list = option_list))

if(!file.exists(args$genotype_file)){
  stop(paste0("Genotype data: ", args$genotype_file, " does not exist."))
}

if(!file.exists(args$annotation_file)){
  stop(paste0("Annotation data: ", args$annotation_file, " does not exist."))
}

if(!file.exists(args$pathway_availability_file)){
  stop(paste0("Pathway Availability data: ", args$pathway_availability_file, " does not exist."))
}

annotations <- fread(args$annotation_file, header = FALSE)

genotypes <- fread(args$genotype_file)

availability_summary <- fread(args$pathway_availability_file)

if(!dir.exists("Reference_SNPFiles")){
  dir.create("Reference_SNPFiles")
}

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

for(i in 1:nrow(availability_summary)){
  
  print(paste0("Currently working on ", availability_summary$PathwayID[i], " (", i, " of ", nrow(availability_summary), ")..."))
  if(file.exists(paste0("Reference_SNPFiles/", availability_summary$PathwayID[i], ".tsv"))){
    print(paste0("Pathway ", pathways$PathwayID[i], " SNP file exists. Skipping..."))
    next
  }
  
  curr_genes <- str_split(availability_summary$AllGenes[i], ";")[[1]]

  gene_details <- mappings[mappings$gene_name %in% curr_genes,]
  gene_details$start <- pmax(gene_details$start - 1000000, 0)
  gene_details$end <- gene_details$end + 1000000
  
  sel_geno <- c()
  
  required_chrs <- unique(gene_details$chromosome)
  required_chrs <- required_chrs[order(as.numeric(required_chrs))]
  
  for(chr in required_chrs){
    print(paste0("Currently working on chromosome ", chr, "..."))
    curr_gene_details <- gene_details[gene_details$chromosome == chr, ]
    
    chr_geno <- genotypes[genotypes$CHR == chr,]
    ranges <- curr_gene_details[, c("start", "end")]
    chr_sel_geno <- chr_geno[chr_geno$LOC %inrange% ranges, c("CHR", "LOC")]
    
    chr_sel_geno <- chr_sel_geno[order(chr_sel_geno$LOC)]
    sel_geno <- rbindlist(list(sel_geno, chr_sel_geno))
  }
  
  print(paste0("Writing pathway genotypes to Reference_SNPFiles/", availability_summary$PathwayID[i], ".tsv"))
  fwrite(sel_geno, paste0("Reference_SNPFiles/", availability_summary$PathwayID[i], ".tsv"), sep = "\t")
}