
# Uses the pacman package to ensure all required packages are installed and loaded
if (!require("pacman")) install.packages("pacman")
pacman::p_load(data.table, 
               stringr,
               optparse,
               MetaSKAT,
               SKAT)

# Defining arguments to be parse from command line
option_list <- list(
  make_option(c("-r", "--results_file"), action = "store", type = "character", 
              default = "Pathway_SummarySKAT_Results",
              help = "Prefix of location to store final results. A CSV file will be generated. 
              Default location is %default."),
  make_option(c("-p", "--phenotype_type"), action = "store", type = "character", 
              default = "D",
              help = "An indicator of the outcome type. 'C' for the continuous outcome 
              and 'D' for the dichotomous outcome. Default value is %default."),
  make_option(c("-s", "--summary_folder"), action = "store", type = "character", 
              default = "SKAT_Summary",
              help = "Folder in which to SKAT summary files for each reduced pathway are stored. 
              Default location is %default.")
)

# Loading parsed arguments
args <- parse_args(OptionParser(option_list = option_list))

if(!dir.exists(args$summary_folder)){
  stop(paste0("Directory ", args$summary_folder, " does not exist!"))
}

MSSD_Files <- list.files(args$summary_folder, pattern = "*.MSSD", full.names = TRUE, recursive = TRUE)

file_count <- length(MSSD_Files)
print(paste0("Detected ", file_count, " MSSD files."))

all_results <- c()

for(i in 1:file_count){
  print(paste0("Currently working on ", f, "(", i, " of ", file_count, ")..."))

  f <- MSSD_Files[i]
  
  curr_dir <- dirname(f)
  dir_split <- str_split(curr_dir, "/", simplify = TRUE)[1, ]
  
  method <- dir_split[3]
  pathway <- dir_split[2]
  
  MInfo.File <- paste0(curr_dir, "/SKAT.MInfo")
  
  weights <- fread(paste0(curr_dir, "/Snp.Weights.txt"))
  
  info_lines <- readLines(MInfo.File)[1:6]
  temp.info <- fread(MInfo.File)
  
  temp.info$MAF <- weights$V2[match(temp.info$SNPID, weights$V1)]
  
  New.MInfo.File <- paste0(curr_dir, "/SKAT.New.MInfo")
  
  writeLines(info_lines, New.MInfo.File)
  
  fwrite(temp.info, New.MInfo.File, append = TRUE, sep = " ")
  
  Curr.Info <- Open_MSSD_File_2Read(f, New.MInfo.File)
  
  curr_results <- MetaSKAT_MSSD_ALL(Curr.Info, combined.weight = FALSE, is.separate = TRUE,
                                    weights.beta = weights$V2, custom.weights = TRUE)
  
  marker_counts <- temp.info[, .N, by = .(SetID)]
  curr_results$N.Marker.All <- marker_counts$N[match(curr_results$SetID, marker_counts$SetID)]
  curr_results$N.Marker.Test <- curr_results$N.Marker.All
  curr_results$pathway <- pathway
  curr_results$num_snp_sets <- nrow(curr_results)
  curr_results$method_dim <- method
  
  all_results <- rbindlist(list(all_results, curr_results))
}

fwrite(all_results, paste0(args$results_file, ".csv"))