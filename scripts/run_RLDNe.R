# Usage:
# Rscript run_RLDNe.R input_file
args <- commandArgs(trailingOnly=TRUE)
if(length(args) < 1) stop("usage: Rscript run_RLDNe.R input_file")
input_file <- args[1]
out_csv <- gsub(".tsv",".out.csv",input_file)

# Loading packages
suppressPackageStartupMessages({
  library(data.table)
  library(RLDNe)
})

# helper: run RLDNe on one file
df <- fread(input_file, na.strings = c("0","NA","."), data.table = FALSE)
efgl <- readInData(df,genotypeStart = 3,pedigreeColumn = 1,nameColumn = 2)
rldne <- exportGenePop_RLDNe(EFGLdata = efgl)
rldne <- create_LDNe_params(rldne,matingsystem=1,crit_vals=c(0.01,0.02,0.05))
std_out <- run_LDNe(rldne)
Ne_out <- read_LDNeOutFile(rldne)
fwrite(Ne_out, out_csv)
cat("Writing：", out_csv, "\n")

