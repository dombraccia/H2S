# ================================= upset.R ================================= #
#' description: script for making the upset plot over each cys degrading enzyme
#' 
#' NOTE: it is assumed that all necessary packages and data are loaded into the 
#' R session where this script is being run.
# =========================================================================== #

## loading in cys degrading gene genome containment
cys_genes <- read.table("../results/from-GutFunFind/Cysteine_Degradation.genes.tsv")
colnames(cys_genes) <- c("genome", "genes")
gene_names <- c("sseA", "mccA", "metB", "dycD", 
                "fn0625", "fn1055", "fn1220", "fn1419", "mgl")
cys_genes <- mutate(cys_genes, sseA = 0, mccA = 0, metB = 0, dycD = 0, 
                    fn0625 = 0, fn1055 = 0, fn1220 = 0, fn1419 = 0, mgl = 0)
# cys_genes <- mutate(cys_genes, MST = 0, CBS = 0, CSE = 0, CYD = 0,
#                     MGL_0625 = 0, MGL_1055 = 0, MGL_1220 = 0, MGL_1419 = 0, MGL = 0)
for (i in 1:dim(cys_genes)[1]) {
  current_genes <- strsplit(cys_genes[i, "genes"], ",")
  for (j in 1:length(current_genes[[1]])) {
    if ("absent" %in% current_genes) {break}
    cys_genes[i, current_genes[[1]][j]] <- 1
  }
}
new_gene_names <- c("MST", "CBS", "CSE", "CYD",
                "MGL_0625", "MGL_1055", "MGL_1220", "MGL_1419", "MGL")
colnames(cys_genes) <- c("genome", "genes", new_gene_names)

## loading dsrAB gene containment
dsr_feature <- read.table("../results/from-GutFunFind/Dissimilatory_Sulf_Reduction.feature.tsv")
colnames(dsr_feature) <- c("genome", "func")
dsr_feature %>% 
  mutate(absent = ifelse(func == "absent", 1, 0),
         dsrA = ifelse((func == "dsrA" | func == "dsrAB"), 1, 0),
         dsrB = ifelse((func == "dsrB" | func == "dsrAB"), 1, 0),
         dsrAB = ifelse(func == "dsrAB", 1, 0)) %>%
  select(genome, absent, dsrA, dsrB, dsrAB) -> dsr_overlap

## TODO: combine "absent" columns from both cys met deg and diss sulf red
cys_dsr_overlap <- cbind(cys_genes[, -2], dsrAB = dsr_overlap$dsrAB) 