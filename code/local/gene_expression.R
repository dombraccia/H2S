# =========================== gene_expression.R ============================= #
#' description: processing and plotting of gene expression data for curated
#' cysteine degrading genes. Diamond (blastx mode) was used to map reads to 
#' custom database of specific genes examined in this paper.
#' 
#' NOTE: it is assumed that all necessary packages and data are loaded into the 
#' R session where this script is being run.
# =========================================================================== #



#### ======================== loading in data =========================== #####
print("- loading and processing H2S data")

## H2S count data
counts_raw_H2S <- read.csv("../../results/diamond_processed/H2S_producing/counts.tsv", 
         sep = "\t", header = TRUE)
rnames <- counts_raw_H2S$Run
counts_H2S <- t(counts_raw_H2S[, 2:21])
colnames(counts_H2S) <- rnames

## CH4 count data
counts_raw_CH4 <- read.csv("../../results/diamond_processed/CH4_producing/counts.tsv", 
                           sep = "\t", header = TRUE)
rnames_CH4 <- counts_raw_CH4$Run
counts_CH4 <- t(counts_raw_CH4[, 2:21])
colnames(counts_CH4) <- rnames_CH4

## sample metadata
metadata <- read.csv("../../data/david-2014/SraRunTable.txt", header = TRUE)
extra_metadata <- read.csv("../../data/from-xiaofang/cib_metadata.tsv", 
                           sep = "\t", header = TRUE)


#### ======================== normalizing counts ======================== #####

print("- normalizing counts to RPK -> TPM")

## getting length of genes in kb
H2S_aa_seqs <- read.fasta("../../data/from-uniprot/H2S_producing_enzymes.fa", seqtype = "AA")
CH4_aa_seqs <- read.fasta("../../data/from-uniprot/CH4_producing_enzymes.fa", seqtype = "AA")

## calculating gene lengths (*3 AA->DNA, /1000 per kilobase)
H2S_seq_len_raw <- sapply(H2S_aa_seqs, length) * 3 / 1000 
CH4_aa_seq_lengths <- sapply(CH4_aa_seqs, length) * 3 / 1000 

## averaging lengths for H2S producing genes (for some genes, multiple 
## homologues were downloaded and included in the search space)
H2S_aa_seq_lengths <- vector("numeric", length = 20)
H2S_aa_seq_lengths[1] <- mean(H2S_seq_len_raw[1:4])
H2S_aa_seq_lengths[2:16] <- H2S_seq_len_raw[5:19]
H2S_aa_seq_lengths[17] <- mean(H2S_seq_len_raw[20:21])
H2S_aa_seq_lengths[18:20] <- H2S_seq_len_raw[22:24]


## scaling raw counts to rpkm
trc_raw <- read.csv("../../data/david-2014/sample_read_counts.tsv",
                sep = "\t", header = FALSE)
colnames(trc_raw) <- c("Run", "num_reads")
trc <- trc_raw$num_reads / 1e6
names(trc) <- trc_raw$Run
norm_counts_H2S <- sweep(counts_H2S, 2, trc, FUN = '/')
norm_counts_CH4 <- sweep(counts_CH4, 2, trc, FUN = '/')

## converting read number normalized counts to RPKM
RPKM_H2S <- sweep(norm_counts_H2S, 1, H2S_aa_seq_lengths, FUN = '/')
RPKM_CH4 <- sweep(norm_counts_CH4, 1, CH4_aa_seq_lengths, FUN = '/')



#### ======================== prep for plotting ========================= #####
print("- prepping data for plotting")

## subsetting metadata based on diet and time
metadata$Time <- as.numeric(gsub("Day ", "", metadata$Time))
metadata %>% filter(Time < 0) -> baseline_metadata
metadata %>% filter(Time > 0, diet == "plant-based") -> plant_metadata
metadata %>% filter(Time > 0, diet == "animal-based") -> animal_metadata

## prepping baseline data
baseline_counts_H2S <- RPKM_H2S[, colnames(RPKM_H2S) %in% baseline_metadata$Run]
baseline_df_H2S <- melt(baseline_counts_H2S)
colnames(baseline_df_H2S) <- c("genes", "Run", "counts")
baseline_df_H2S$genes <- factor(baseline_df_H2S$genes, 
    levels =  c("mst", "cse", "cbs", "cyd", "metC", "cysM", "cysK",  "malY", "yhaO", "yhaM",
                "decR", "iscS", "fn0625", "fn1055", "fn1220", "fn1419", "mgl", "tnaA", "dsrA", "dsrB"))

plant_counts_H2S <- RPKM_H2S[, colnames(RPKM_H2S) %in% plant_metadata$Run]
plant_df_H2S <- melt(plant_counts_H2S)
colnames(plant_df_H2S) <- c("genes", "Run", "counts")
plant_df_H2S$genes <- factor(plant_df_H2S$genes, 
    levels =  c("mst", "cse", "cbs", "cyd", "metC", "cysM", "cysK",  "malY", "yhaO", "yhaM",
                "decR", "iscS", "fn0625", "fn1055", "fn1220", "fn1419", "mgl", "tnaA", "dsrA", "dsrB"))

animal_counts_H2S <- RPKM_H2S[, colnames(RPKM_H2S) %in% animal_metadata$Run]
animal_df_H2S <- melt(animal_counts_H2S)
colnames(animal_df_H2S) <- c("genes", "Run", "counts")
animal_df_H2S$genes <- factor(animal_df_H2S$genes, 
    levels =  c("mst", "cse", "cbs", "cyd", "metC", "cysM", "cysK",  "malY", "yhaO", "yhaM",
                "decR", "iscS", "fn0625", "fn1055", "fn1220", "fn1419", "mgl", "tnaA", "dsrA", "dsrB"))


#### ========================== CH4 PRODUCTION =========================== ####

print("- processing CH4 data")

## prepping baseline data
baseline_counts_CH4 <- RPKM_CH4[, colnames(RPKM_CH4) %in% baseline_metadata$Run]
baseline_df_CH4 <- melt(baseline_counts_CH4)
colnames(baseline_df_CH4) <- c("genes", "Run", "counts")
baseline_df_CH4$genes <- factor(baseline_df_CH4$genes, 
                                levels =  c("fwdA", "fwdB", "fwdD", "fwdE", "fwdF", "fwdG",
                                            "Ftr", "Hmd", "Mch", "mcrA", "mcrB", "mcrG", "Mer", 
                                            "mtrA", "mtrB", "mtrC", "mtrD", "mtrE", "mtrF", "mtrG"))

plant_counts_CH4 <- RPKM_CH4[, colnames(RPKM_CH4) %in% plant_metadata$Run]
plant_df_CH4 <- melt(plant_counts_CH4)
colnames(plant_df_CH4) <- c("genes", "Run", "counts")
plant_df_CH4$genes <- factor(plant_df_CH4$genes, 
                             levels =  c("fwdA", "fwdB", "fwdD", "fwdE", "fwdF", "fwdG",
                                         "Ftr", "Hmd", "Mch", "mcrA", "mcrB", "mcrG", "Mer", 
                                         "mtrA", "mtrB", "mtrC", "mtrD", "mtrE", "mtrF", "mtrG"))

animal_counts_CH4 <- RPKM_CH4[, colnames(RPKM_CH4) %in% animal_metadata$Run]
animal_df_CH4 <- melt(animal_counts_CH4)
colnames(animal_df_CH4) <- c("genes", "Run", "counts")
animal_df_CH4$genes <- factor(animal_df_CH4$genes, 
                              levels =  c("fwdA", "fwdB", "fwdD", "fwdE", "fwdF", "fwdG",
                                          "Ftr", "Hmd", "Mch", "mcrA", "mcrB", "mcrG", "Mer", 
                                          "mtrA", "mtrB", "mtrC", "mtrD", "mtrE", "mtrF", "mtrG"))


#### ============================ plotting ============================== #####
print("- plotting H2S plots")

## baseline boxplot
# baseline_df_H2S$counts <- baseline_df_H2S$counts + 1 
baseline_boxplot_H2S <- ggplot(data = baseline_df_H2S, aes(x = genes, fill = genes, y = counts)) +
  geom_boxplot(outlier.shape = NA) + 
  # geom_jitter(alpha = 0.8, size = 0.4, width = 0.2) + 
  theme_bw() +
  theme(legend.position = "none",
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.title.x = element_text(size = 16),
        axis.title.y = element_text(size = 16)) +
  scale_y_continuous(limits = c(0, 20), 
                     # breaks = c(1, 10, 100, 1000, 5000),
                     #trans = "log2",
                     minor_breaks = NULL) +
  xlab("baseline") + ylab("RPKM")
baseline_boxplot_H2S

## plant boxplot
plant_boxplot_H2S <- ggplot(data = plant_df_H2S, aes(x = genes, fill = genes, y = counts)) +
  geom_boxplot(outlier.shape = NA) + 
  # geom_jitter(alpha = 0.8, size = 0.4, width = 0.2) +
  theme_bw() +
  theme(legend.position = "none",
        axis.title.y = element_blank(),
        axis.text.x = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks.x = element_blank(),
        axis.ticks.y = element_blank(),
        axis.title.x = element_text(size = 16)) +
  scale_y_continuous(limits = c(0, 20), 
                     # breaks = c(1, 10, 100, 1000, 5000),
                     #trans = "log2",
                     minor_breaks = NULL) +
  xlab("plant-based")

## animal boxplot
animal_boxplot_H2S <- ggplot(data = animal_df_H2S, aes(x = genes, fill = genes, y = counts)) +
  geom_boxplot(outlier.shape = NA) + 
  # geom_jitter(alpha = 0.8, size = 0.4, width = 0.2) +
  theme_bw() +
  theme(axis.title.y = element_blank(),
        axis.text.y = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.ticks.y = element_blank(),
        axis.title.x = element_text(size = 16)) +
  scale_y_continuous(limits = c(0, 20), 
                     # breaks = c(1, 10, 100, 1000, 5000),
                     #trans = "log2",
                     minor_breaks = NULL) +
  xlab("animal-based")

baseline_boxplot_H2S | plant_boxplot_H2S | animal_boxplot_H2S
ggsave("../../figures/supplementary/david2014_RNAseq_diet_H2S_genes.svg", height = 7, width = 10)
ggsave("../../figures/supplementary/david2014_RNAseq_diet_H2S_genes.png", height = 7, width = 10)

#### ========================== plotting CH4 ============================ #####
print("- plotting CH4 data")

## baseline boxplot
baseline_boxplot_CH4 <- ggplot(data = baseline_df_CH4, aes(x = genes, fill = genes, y = counts)) +
  geom_boxplot(outlier.shape = NA) + 
  # geom_jitter(alpha = 0.8, size = 0.4, width = 0.2) + 
  theme_bw() +
  theme(legend.position = "none",
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.title.x = element_text(size = 16),
        axis.title.y = element_text(size = 16)) +
  scale_y_continuous(limits = c(0, 5), 
                     #trans = "log2",
                     minor_breaks = NULL) +
  xlab("baseline") + ylab("RPKM")

## plant boxplot
plant_boxplot_CH4 <- ggplot(data = plant_df_CH4, aes(x = genes, fill = genes, y = counts)) +
  geom_boxplot(outlier.shape = NA) + 
  # geom_jitter(alpha = 0.8, size = 0.4, width = 0.2) +
  theme_bw() +
  theme(legend.position = "none",
        axis.title.y = element_blank(),
        axis.text.x = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks.x = element_blank(),
        axis.ticks.y = element_blank(),
        axis.title.x = element_text(size = 16)) +
  scale_y_continuous(limits = c(0, 5),
                     #trans = "log2",
                     minor_breaks = NULL) +
  xlab("plant-based")

## animal boxplot
animal_boxplot_CH4 <- ggplot(data = animal_df_CH4, aes(x = genes, fill = genes, y = counts)) +
  geom_boxplot(outlier.shape = NA) + 
  # geom_jitter(alpha = 0.8, size = 0.4, width = 0.2) +
  theme_bw() +
  theme(axis.title.y = element_blank(),
        axis.text.y = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.ticks.y = element_blank(),
        axis.title.x = element_text(size = 16)) +
  scale_y_continuous(limits = c(0, 5),
                     #trans = "log2",
                     minor_breaks = NULL) +
  xlab("animal-based")

baseline_boxplot_CH4 | plant_boxplot_CH4 | animal_boxplot_CH4
ggsave("../../figures/supplementary/david2014_RNAseq_diet_CH4_genes.svg", height = 7, width = 10)
ggsave("../../figures/supplementary/david2014_RNAseq_diet_CH4_genes.png", height = 7, width = 10)