# ============================ gene_hits_table.R ============================ #
#' description: script for generating a table showing all of the gene hits to 
#' UHGG genomes. going into a table in the manuscript
#' 
#' NOTE: it is assumed that all necessary packages and data are loaded into the 
#' R session where this script is being run.
# =========================================================================== #

print("- reading in data file")
gene_hits <- read.csv("../../results/from-GutFunFind/Cysteine_Degradation.genes_table.txt",
                      sep = '\t', header = TRUE)

