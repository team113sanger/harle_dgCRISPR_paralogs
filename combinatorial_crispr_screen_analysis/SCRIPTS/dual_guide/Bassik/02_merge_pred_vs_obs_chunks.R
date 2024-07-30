suppressPackageStartupMessages(library(optparse))
suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(vroom))
suppressPackageStartupMessages(library(data.table))

###############################################################################
#* --                                                                     -- *#
#* --                             OPTIONS                                 -- *#
#* --                                                                     -- *#
###############################################################################

option_list = list(
  make_option(c("-f", "--files"), type = "character",
               help = "files to merge (file of file paths one per line)", metavar = "character"),
  make_option(c("-p", "--prefix"), type = "character",
               help = "output file prefix", metavar = "character"),
  make_option(c("-o", "--out"), type = "character", default = '.',
               help = "TSV output directory [Default: . ]", metavar = "character"),
  make_option(c("-r", "--rds"), type = "character", default = '.',
               help = "RDS output directory [Default: . ]", metavar = "character")
);

opt_parser <- OptionParser(option_list = option_list);
opt <- parse_args(opt_parser);

###############################################################################
#* --                                                                     -- *#
#* --                             VALIDATION                              -- *#
#* --                                                                     -- *#
###############################################################################

# Fold change matrix file
if (is.null(opt$files)) {
  print_help(opt_parser)
  stop("Please provide a comma-separated list of files", call.=FALSE)
}

# Check files exist
if (! file.exists(opt$files)) {
  stop(paste("File does not exist:", opt$files), call.=FALSE)
}

# Prefix
if (is.null(opt$prefix)) {
  print_help(opt_parser)
  stop("Please provide result file prefix", call.=FALSE)
}

###############################################################################
#* --                                                                     -- *#
#* --                              FUNCTIONS                              -- *#
#* --                                                                     -- *#
###############################################################################

generate_output_filepath <- function(filename, outdir = '.') {
  output_filepath <- file.path(outdir, filename)
  return(output_filepath)
}

###############################################################################
#* --                                                                     -- *#
#* --                          MAIN SCRIPT                                -- *#
#* --                                                                     -- *#
###############################################################################

message("Reading list of files...")
files_to_merge <- scan(file = opt$files, what = character())
dfs_to_merge <- list()

for (fn in files_to_merge) {
  dfs_to_merge[[basename(fn)]] <- vroom(file = fn , delim = "\t", col_names = T, show_col_types = FALSE)
}

message("Start merging...")
merged_df <- data.frame()
merged_df <- rbindlist(dfs_to_merge, fill = T)
message("Done merging...")

message("Writing TSV to file...")
merged_df_filename = generate_output_filepath(paste(opt$prefix, "tsv", sep = "."), opt$out)
write.table(merged_df, file = merged_df_filename, row.names = F, sep = "\t", quote = F)
message(paste("TSV written to:", merged_df_filename))

message("Writing RDS to file...")
merged_rds_filename = generate_output_filepath(paste(opt$prefix, "rds", sep = "."), opt$rds)
saveRDS(merged_df, file = merged_rds_filename)
message(paste("RDS written to:", merged_rds_filename))