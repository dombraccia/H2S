# ============================ gene_hits_hmm_table.R ============================ #
#' description: script for generating a table showing all of the gene hits to 
#' UHGG genomes. going into a table in the manuscript
#' 
#' NOTE: it is assumed that all necessary packages and data are loaded into the 
#' R session where this script is being run.
# =========================================================================== #

print("- reading in data files") 

## code below taken from original upset.R script
cys_genes <- read.table("results/from-GutFunFind/Cysteine_Degradation_hmm.genes_hmm.tsv",
                        header = TRUE)
colnames(cys_genes) <- c("genome", "genes")
gene_names <- c("dcyD", "yhaO", "yhaM", "mgl", "aspC", "sseA", "metC", 
                "malY", "cysK", "cysM", "mccB", "tnaA", "iscS", "mccA",
                "absent")
cys_genes <- mutate(cys_genes, dcyD = 0, yhaO = 0, yhaM = 0, mgl = 0, aspC = 0, sseA = 0, metC = 0,  
                    malY = 0, cysK = 0, cysM = 0, mccB = 0, tnaA = 0, iscS = 0, mccA = 0, absent = 0)
for (i in 1:dim(cys_genes)[1]) {
  current_genes <- strsplit(cys_genes[i, "genes"], ",")
  gene_table <- table(current_genes)
  for (j in 1:length(gene_table)) {
    if ("absent" %in% names(current_genes)) {break}
    cys_genes[i, names(gene_table)[j]] <- gene_table[j]
  }
}
# select(cys_genes, -c(absent)) # removing "absent" column

## generating dsr_overlap for future use
dsr_feature <- read.table("results/from-GutFunFind/Dissimilatory_Sulf_Reduction.feature.tsv")
colnames(dsr_feature) <- c("genome", "func")
dsr_feature %>% 
  mutate(absent = ifelse(func == "absent", 1, 0),
         dsrA = ifelse((func == "dsrA" | func == "dsrAB"), 1, 0),
         dsrB = ifelse((func == "dsrB" | func == "dsrAB"), 1, 0),
         dsrAB = ifelse(func == "dsrAB", 1, 0)) %>%
  select(genome, absent, dsrA, dsrB, dsrAB) -> dsr_overlap

