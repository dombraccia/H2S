# ================== metatranscriptomic_visualization.R ==================== #
#' description: plotting results of transcriptomic gene quantification analysis
#' 
#' NOTE: it is assumed that all necessary packages and data are loaded into the 
#' R session where this script is being run.
#' 
#' Author: Domenick J. Braccia
#' Last updated: 2021-09-20
# =========================================================================== #

# loading libraries
library(dplyr)
library(magrittr)
library(reshape2)
library(scales)
library(ggplot2)
library(patchwork)
library(cowplot)

# load in data

hpfs_cd_raw_counts <- read.csv("results/salmon_processed/hpfs/cys_counts.tsv", sep = "\t", header = TRUE)
hpfs_dsr_raw_counts <- read.csv("results/salmon_processed/hpfs/DSR_counts.tsv", sep = "\t", header = TRUE)
hpfs_ch4_raw_counts <- read.csv("results/salmon_processed/hpfs/ch4_counts.tsv",sep = "\t", header = TRUE)
david2014_cd_raw_counts <- read.csv("results/salmon_processed/david2014/cys_counts.tsv", sep = "\t", header = TRUE)
david2014_dsr_raw_counts <- read.csv("results/salmon_processed/david2014/DSR_counts.tsv", sep = "\t", header = TRUE)
david2014_ch4_raw_counts <- read.csv("results/salmon_processed/david2014/ch4_counts.tsv", sep = "\t", header = TRUE)

hpfs_info <- read.csv("data/hpfs/nreads.txt", sep = " ", header = FALSE)
david2014_info <- read.csv("data/david2014/nreads.txt", sep = " ", header = FALSE)

# process count data and metadata
colnames(hpfs_info) <- c("run", "nreads")
colnames(david2014_info) <- c("run", "nreads")
hpfs_info$run <- substr(hpfs_info$run, 1, 10)
david2014_info$run <- substr(david2014_info$run, 1, 9)

## cysteine degradation genes
cd_gene_names <- sapply(strsplit(hpfs_cd_raw_counts$Name, ';'), "[[", 4)
cd_info <- data.frame(gene_info = hpfs_cd_raw_counts$Name,
                      gene_name = cd_gene_names,
                      length = hpfs_cd_raw_counts$Length,
                      efflength = hpfs_cd_raw_counts$EffectiveLength)

hpfs_cd_raw_counts %>% select(matches("SRR*")) -> hpfs_cd_counts
david2014_cd_raw_counts %>% select(matches("*SRR")) -> david2014_cd_counts

## dissimilatory sulfate reducing genes
dsr_gene_names <- sapply(strsplit(hpfs_dsr_raw_counts$Name, ';'), "[[", 4)
dsr_info <- data.frame(gene_info = hpfs_dsr_raw_counts$Name,
                      gene_name = dsr_gene_names,
                      length = hpfs_dsr_raw_counts$Length,
                      efflength = hpfs_dsr_raw_counts$EffectiveLength)

hpfs_cd_raw_counts %>% select(matches("SRR*")) -> hpfs_dsr_counts
david2014_cd_raw_counts %>% select(matches("*SRR")) -> david2014_dsr_counts


## ch4 production genes
ch4_gene_names <- sapply(strsplit(hpfs_ch4_raw_counts$Name, ';'), "[[", 2)
ch4_info <- data.frame(gene_info = hpfs_ch4_raw_counts$Name,
                       gene_name = ch4_gene_names,
                       length = hpfs_ch4_raw_counts$Length,
                       efflength = hpfs_ch4_raw_counts$EffectiveLength)

hpfs_ch4_raw_counts %>% select(matches("SRR*")) -> hpfs_ch4_counts
david2014_ch4_raw_counts %>% select(matches("SRR*")) -> david2014_ch4_counts

## DSR genes
dsr_gene_names <- sapply(strsplit(hpfs_dsr_raw_counts$Name, ';'), "[[", 4)
dsr_info <- data.frame(gene_info = hpfs_dsr_raw_counts$Name,
                      gene_name = dsr_gene_names,
                      length = hpfs_dsr_raw_counts$Length,
                      efflength = hpfs_dsr_raw_counts$EffectiveLength)

hpfs_dsr_raw_counts %>% select(matches("SRR*")) -> hpfs_dsr_counts
david2014_dsr_raw_counts %>% select(matches("*SRR")) -> david2014_dsr_counts


# only keep info for samples which were quantified & reorder
hpfs_info <- hpfs_info[hpfs_info$run %in% colnames(hpfs_cd_counts),]
david2014_info <- david2014_info[david2014_info$run %in% colnames(david2014_cd_counts), ]
hpfs_cd_counts <- hpfs_cd_counts[hpfs_info$run]
david2014_cd_counts <- david2014_cd_counts[david2014_info$run]

# normalize counts to TPM

# 1. convert raw counts -> RPK
hpfs_cd_rpk <- sweep(x = hpfs_cd_counts, MARGIN = 1, STATS = (cd_info$length / 1e3), FUN = "/")
david2014_cd_rpk <- sweep(x = david2014_cd_counts, MARGIN = 1, STATS = (cd_info$length / 1e3), FUN = "/")
hpfs_dsr_rpk <- sweep(x = hpfs_dsr_counts, MARGIN = 1, STATS = (dsr_info$length / 1e3), FUN = "/")
david2014_dsr_rpk <- sweep(x = david2014_dsr_counts, MARGIN = 1, STATS = (dsr_info$length / 1e3), FUN = "/")
hpfs_ch4_rpk <- sweep(x = hpfs_ch4_counts, MARGIN = 1, STATS = (ch4_info$length / 1e3), FUN = "/")
david2014_ch4_rpk <- sweep(x = david2014_ch4_counts, MARGIN = 1, STATS = (ch4_info$length / 1e3), FUN = "/")

# 2. divide by scaling factor (nreads/1e6)
hpfs_cd_tpm <- sweep(x = hpfs_cd_rpk, MARGIN = 2, STATS = (hpfs_info$nreads/1e6), FUN = "/")
david2014_cd_tpm <- sweep(x = david2014_cd_rpk, MARGIN = 2, STATS = (david2014_info$nreads/1e6), FUN = "/")
hpfs_dsr_tpm <- sweep(x = hpfs_dsr_rpk, MARGIN = 2, STATS = (hpfs_info$nreads/1e6), FUN = "/")
david2014_dsr_tpm <- sweep(x = david2014_dsr_rpk, MARGIN = 2, STATS = (david2014_info$nreads/1e6), FUN = "/")
hpfs_ch4_tpm <- sweep(x = hpfs_ch4_rpk, MARGIN = 2, STATS = (hpfs_info$nreads/1e6), FUN = "/")
david2014_ch4_tpm <- sweep(x = david2014_ch4_rpk, MARGIN = 2, STATS = (david2014_info$nreads/1e6), FUN = "/")

# combine tpm of homologous genes
hpfs_cd_tpm$gene_name <- cd_info$gene_name
hpfs_cd_tpm %>% group_by(gene_name) %>% summarise_if(is.numeric, sum) -> hpfs_cd_tpm
hpfs_dsr_tpm$gene_name <- dsr_info$gene_name
hpfs_dsr_tpm %>% group_by(gene_name) %>% summarise_if(is.numeric, sum) -> hpfs_dsr_tpm
hpfs_ch4_tpm$gene_name <- ch4_info$gene_name
hpfs_ch4_tpm %>% group_by(gene_name) %>% summarise_if(is.numeric, sum) -> hpfs_ch4_tpm

david2014_cd_tpm$gene_name <- cd_info$gene_name
david2014_cd_tpm %>% group_by(gene_name) %>% summarise_if(is.numeric, sum) -> david2014_cd_tpm
david2014_dsr_tpm$gene_name <- dsr_info$gene_name
david2014_dsr_tpm %>% group_by(gene_name) %>% summarise_if(is.numeric, sum) -> david2014_dsr_tpm
david2014_ch4_tpm$gene_name <- ch4_info$gene_name
david2014_ch4_tpm %>% group_by(gene_name) %>% summarise_if(is.numeric, sum) -> david2014_ch4_tpm

# cbind-ing cysteine degradation genes and dsr genes for plotting
hpfs_h2s_tpm <- rbind(hpfs_cd_tpm, hpfs_dsr_tpm)
david2014_h2s_tpm <- rbind(david2014_cd_tpm, david2014_dsr_tpm)

#  prepping tpm dataframes 
hpfs_h2s_tpm_long <- melt(hpfs_h2s_tpm, variable.name = "run", value.name = "tpm", id = c("gene_name"))
hpfs_h2s_tpm_long$gene_name <- factor(hpfs_h2s_tpm_long$gene_name, 
                                     levels = c("dcyD", "yhaO", "yhaM", "mgl", "aspC", "sseA", "malY", 
                                                "metC", "cysK", "cysM", "mccB", "tnaA", "iscS", "mccA",
                                                "dsrA", "dsrB"))
hpfs_ch4_tpm_long <- melt(hpfs_ch4_tpm, variable.name = "run", value.name = "tpm", id = c("gene_name"))
david2014_h2s_tpm_long <- melt(david2014_h2s_tpm, variable.name = "run", value.name = "tpm", id = c("gene_name"))
david2014_h2s_tpm_long$gene_name <- factor(david2014_h2s_tpm$gene_name, 
                                          levels = c("dcyD", "yhaO", "yhaM", "mgl", "aspC", "sseA", "malY", 
                                                     "metC", "cysK", "cysM", "mccB", "tnaA", "iscS", "mccA",
                                                     "dsrA", "dsrB"))
david2014_ch4_tpm_long <- melt(david2014_ch4_tpm, variable.name = "run", value.name = "tpm", id = c("gene_name"))

# plotting normalized tpm values 
clrs1 <- c("#93c47d", "#999999", "#93c47d", "#93c47d", "#999999", "#93c47d", rep("#ffd966", 5), rep("#e06666", 3), rep("#5f89dd", 2))
h2s_gene_names <- as.character(unique(hpfs_h2s_tpm_long$gene_name))
p1 <- ggplot(hpfs_h2s_tpm_long, aes(x = gene_name, y = tpm, fill = gene_name)) + 
  geom_boxplot() +
  scale_y_continuous(limits = c(0, 400)) +
  # scale_y_discrete(breaks = c(0.5)) + 
  scale_fill_manual(values = clrs1) +
  theme_bw() +
  theme(legend.text = element_text(face = "italic"),
        axis.title.x = element_blank(),
        panel.grid.minor = element_blank(),
        panel.grid.major = element_blank(),
        axis.text.x = element_text(size = 12, angle = 45, hjust = 0.95),
        axis.text.y = element_text(size = 12),
        axis.title.y = element_text(size = 12)
        ) +
  ggtitle("Health Professionals Follow-up Study (HPFS)") + ylab("transcripts per million (TPM)")
p1

p2 <- ggplot(david2014_h2s_tpm_long, aes(x = gene_name, y = tpm, fill = gene_name)) + 
  geom_boxplot() +
  scale_y_continuous(limits = c(0, 100)) +
  scale_fill_manual(values = clrs1) +
  theme_bw() +
  theme(legend.position = "none",
        axis.title.x = element_blank(),
        panel.grid.minor = element_blank(),
        panel.grid.major = element_blank(),
        axis.text.x = element_text(size = 12, angle = 45, hjust = 0.95),
        axis.text.y = element_text(size = 12),
        axis.title.y = element_text(size = 12)
        ) +
  ggtitle("David et al. 2014")+ ylab("transcripts per million (TPM)")
p2

figure4 <- (p1 / p2) + plot_layout(guides = "collect")
figure4
ggsave("figures/figure4/cys_deg_gene_expression.svg", height = 8, width = 8, units = "in")
ggsave("figures/figure4/cys_deg_gene_expression.png", height = 8, width = 8, units = "in")


###### ================================================================== #####


# hpfs_ch4_tpm_long$tpm <- hpfs_ch4_tpm_long$tpm + 1
p3 <- ggplot(hpfs_ch4_tpm_long, aes(x = gene_name, y = tpm)) + 
  geom_boxplot() +
  scale_y_continuous(labels = comma) +
  theme_bw() +
  theme(axis.title.x = element_blank(),
        panel.grid.minor = element_blank(),
        panel.grid.major = element_blank(),
        axis.text.x = element_text(size = 12, angle = 45, hjust = 0.95),
        axis.text.y = element_text(size = 12),
        axis.title.y = element_text(size = 12)
  ) +
  ggtitle("Health Professionals Follow-up Study (HPFS)")+ ylab("transcripts per million (TPM)")
p3

p4 <- ggplot(david2014_ch4_tpm_long, aes(x = gene_name, y = tpm)) + 
  geom_boxplot() +
  theme_bw() +
  theme(axis.title.x = element_blank(),
        panel.grid.minor = element_blank(),
        panel.grid.major = element_blank(),
        axis.text.x = element_text(size = 12, angle = 45, hjust = 0.95),
        axis.text.y = element_text(size = 12),
        axis.title.y = element_text(size = 12)
  ) +
  ggtitle("David et al. 2014")+ ylab("transcripts per million (TPM)")
p4

figureS3 <- (p3 / p4) & theme(legend.justification = "center")
figureS3
ggsave("figures/supplementary/figureS3.png", height = 8, width = 8, units = "in")
ggsave("figures/supplementary/figureS3.svg", height = 8, width = 8, units = "in")

##### ================ additional statistic calculations ================ #####

## pct of samples with expression of at least one primary cysteine degrading gene.
hpfs_h2s_tpm %>% filter(gene_name == "yhaM" | 
                          gene_name == "mgl" | 
                          gene_name == "dcyD" | 
                          gene_name == "sseA") -> hpfs_primary_tpm

david2014_h2s_tpm %>% filter(gene_name == "yhaM" | 
                              gene_name == "mgl" | 
                              gene_name == "dcyD" | 
                              gene_name == "sseA") -> david2014_primary_tpm

david2014_primary_tpm_cs <- colSums(david2014_primary_tpm[,-1])
hpfs_primary_tpm_cs <- colSums(hpfs_primary_tpm[,-1])
combined_primary_tpm_cs <- c(hpfs_primary_tpm_cs, david2014_primary_tpm_cs)
num_primary_expr <- sum(combined_primary_tpm_cs > 10)
num_total <- length(combined_primary_tpm_cs)
pct_primary_expr <- round(num_primary_expr / num_total * 100, 2)

print(paste0("- RESULT: pct of samples with expression of at least one primary cysteine degrading gene: ",
             pct_primary_expr, " (", num_primary_expr, "/", num_total, ")"))

## pct of samples with expression of at least one secondary cysteine degrading gene.
hpfs_h2s_tpm %>% filter(gene_name == "malY" | 
                          gene_name == "mccB" | 
                          gene_name == "metC" | 
                          gene_name == "cysK" |
                          gene_name == "cysM"
                        ) -> hpfs_secondary_tpm
david2014_h2s_tpm %>% filter(gene_name == "malY" | 
                          gene_name == "mccB" | 
                          gene_name == "metC" | 
                          gene_name == "cysK" |
                          gene_name == "cysM") -> david2014_secondary_tpm


david2014_secondary_tpm_cs <- colSums(david2014_secondary_tpm[,-1])
hpfs_secondary_tpm_cs <- colSums(hpfs_secondary_tpm[,-1])
combined_secondary_tpm_cs <- c(hpfs_secondary_tpm_cs, david2014_secondary_tpm_cs)
num_secondary_expr <- sum(combined_secondary_tpm_cs > 10)
num_total <- length(combined_secondary_tpm_cs)
pct_secondary_expr <- round(num_secondary_expr / num_total * 100, 2)

print(paste0("- RESULT: pct of samples with expression of at least one primary cysteine degrading gene: ",
             pct_secondary_expr, " (", num_secondary_expr, "/", num_total, ")"))

## pct of samples with expression of both dsr genes and methanogenesis genes
### prepping tpm dataframes for comparison
dsr_gene_names <- hpfs_dsr_tpm$gene_name
hpfs_dsr_tpm <- as.data.frame(t(hpfs_dsr_tpm[,-1]))
colnames(hpfs_dsr_tpm) <- dsr_gene_names

ch4_gene_names <- hpfs_ch4_tpm$gene_name
hpfs_ch4_tpm <- as.data.frame(t(hpfs_ch4_tpm[,-1]))
colnames(hpfs_ch4_tpm) <- ch4_gene_names

david2014_dsr_tpm <- as.data.frame(t(david2014_dsr_tpm[,-1]))
colnames(david2014_dsr_tpm) <- dsr_gene_names

david2014_ch4_tpm <- as.data.frame(t(david2014_ch4_tpm[,-1]))
colnames(david2014_ch4_tpm) <- ch4_gene_names

### filtering results to find overlap in DSR + ch4 prod genes
hpfs_all_samples <- rownames(hpfs_dsr_tpm)
hpfs_num_samples <- length(hpfs_all_samples)

hpfs_dsr_tpm %>% 
  filter(dsrA > 0, dsrB > 0) %>% rownames() -> hpfs_dsr_expr_samples

### require that at least 80% of methanogenesis genes are expressed per sample
hpfs_chr_expr_samples <- character()
for (i in 1:nrow(hpfs_ch4_tpm)) {
  if (sum(hpfs_ch4_tpm[i,] > 0) >= (ncol(hpfs_ch4_tpm) * 0.8)) {
    hpfs_chr_expr_samples <- append(hpfs_chr_expr_samples, rownames(hpfs_ch4_tpm[i,]))
  }
}

hpfs_dsr_ch4_expr_num <- length(intersect(hpfs_dsr_expr_samples, hpfs_chr_expr_samples))
hpfs_dsr_ch4_expr_pct <- round(hpfs_dsr_ch4_expr_num / hpfs_num_samples * 100, 2)

print(paste0("- RESULT: the pct of hpfs samples which show expression of both dsr and ch4 producing genes: ",
             hpfs_dsr_ch4_expr_pct, " (", hpfs_dsr_ch4_expr_num, "/", hpfs_num_samples, ")"))

david2014_all_samples <- rownames(david2014_dsr_tpm)
david2014_num_samples <- length(david2014_all_samples)

david2014_dsr_tpm %>% 
  filter(dsrA > 0, dsrB > 0) %>% rownames() -> david2014_dsr_expr_samples

### require that at least 80% of methanogenesis genes are expressed per sample
david2014_chr_expr_samples <- character()
for (i in 1:nrow(david2014_ch4_tpm)) {
  if (sum(david2014_ch4_tpm[i,] > 0) >= (ncol(david2014_ch4_tpm) * 0.8)) {
    david2014_chr_expr_samples <- append(david2014_chr_expr_samples, rownames(david2014_ch4_tpm[i,]))
  }
}

david2014_dsr_ch4_expr_num <- length(intersect(david2014_dsr_expr_samples, david2014_chr_expr_samples))
david2014_dsr_ch4_expr_pct <- round(david2014_dsr_ch4_expr_num / length(david2014_all_samples) * 100, 2)

print(paste0("- RESULT: the pct of david2014 samples which show expression of both dsr and ch4 producing genes: ",
             david2014_dsr_ch4_expr_pct, " (", david2014_dsr_ch4_expr_num, "/", david2014_num_samples, ")"))

## getting just pct of samples expressing DSR genes (for Results section)
hpfs_dsr_expr_pct <- round(length(hpfs_dsr_expr_samples) / length(hpfs_all_samples) * 100, 2)
print(paste0("- RESULT: the pct of hpfs samples which show expression of DSR genes: ",
             hpfs_dsr_expr_pct, " (", length(hpfs_dsr_expr_samples), "/", length(hpfs_all_samples), ")"))

david_dsr_expr_pct <- round(length(david2014_dsr_expr_samples) / length(david2014_all_samples) * 100, 2)
print(paste0("- RESULT: the pct of david2014 samples which show expression of DSR genes: ",
             hpfs_dsr_expr_pct, " (", length(david2014_dsr_expr_samples), "/", length(david2014_all_samples), ")"))

total_num_dsr_expr_samples <- length(hpfs_dsr_expr_samples) + length(david2014_dsr_expr_samples)
total_num_samples <- length(hpfs_all_samples) + length(david2014_all_samples)
dsr_expr_pct <- round(total_num_dsr_expr_samples / total_num_samples * 100, 2)
print(paste0("- RESULT: the pct of hpfs AND david2014 samples which show expression of DSR genes: ",
             dsr_expr_pct, " (", total_num_dsr_expr_samples, "/", total_num_samples, ")"))


##### ================= SUPPLEMENTARY FIG 4 20210910 ==================== #####

## prepping david2014_ch4_tpm data for plotting
david2014_ch4_tpm_temp <- as.data.frame(t(david2014_ch4_tpm))
david2014_ch4_tpm_temp$gene_name <- rownames(david2014_ch4_tpm_temp)
david2014_ch4_tpm_temp <- david2014_ch4_tpm_temp[, c(ncol(david2014_ch4_tpm_temp), 1:(ncol(david2014_ch4_tpm_temp) - 1))]
david2014_ch4_tpm_temp <- david2014_ch4_tpm_temp[, match(colnames(david2014_h2s_tpm), colnames(david2014_ch4_tpm_temp))]## reorder using match

david2014_tpm <- rbind(david2014_h2s_tpm, david2014_ch4_tpm_temp)
david2014_tpm$pathway <- c("primary", rep("secondary", 2), "primary", "erroneous", "secondary", "erroneous", "secondary", "secondary", rep("primary", 2), "erroneous", "primary", "upstream process", rep("DSR", 2), rep("methanogenesis", 16))
david2014_tpm$pathway <- factor(david2014_tpm$pathway, levels = c("primary", "secondary", "erroneous", "DSR", "methanogenesis", "upstream process"))
david2014_tpm_long <- melt(david2014_tpm, variable.name = "sample", value.name = "tpm", id = c("gene_name", "pathway"))
david2014_tpm_long %>% group_by(sample) %>% summarise(tpm_sum = sum(tpm)) %>% arrange(desc(tpm_sum)) -> david2014_tpm_sum
david2014_relevel <- as.character(david2014_tpm_sum$sample)

david2014_tpm_long$sample <- as.character(david2014_tpm_long$sample)
david2014_tpm_long$sample <- factor(david2014_tpm_long$sample, levels = c(david2014_relevel))

## prepping hpfs data for plotting
hpfs_ch4_tpm_temp <- as.data.frame(t(hpfs_ch4_tpm))
hpfs_ch4_tpm_temp$gene_name <- rownames(hpfs_ch4_tpm_temp)
hpfs_ch4_tpm_temp <- hpfs_ch4_tpm_temp[,c(ncol(hpfs_ch4_tpm_temp),1:(ncol(hpfs_ch4_tpm_temp) - 1))]
hpfs_ch4_tpm_temp <- hpfs_ch4_tpm_temp[, match(colnames(hpfs_h2s_tpm), colnames(hpfs_ch4_tpm_temp))]## reorder using match

hpfs_tpm <- rbind(hpfs_h2s_tpm, hpfs_ch4_tpm_temp)
pathways <- c("primary", rep("secondary", 2), "primary", "erroneous", "secondary", "erroneous", "secondary", "secondary", rep("primary", 2), "erroneous", "primary", "upstream process", rep("DSR", 2), rep("methanogenesis", 16))
hpfs_tpm$pathway <- pathways
hpfs_tpm$pathway <- factor(hpfs_tpm$pathway, levels = c("primary", "secondary", "erroneous", "DSR", "methanogenesis", "upstream process"))
hpfs_tpm_long <- melt(hpfs_tpm, variable.name = "sample", value.name = "tpm", id = c("gene_name", "pathway"))
hpfs_tpm_long %>% 
  group_by(sample) %>% 
  summarise(tpm_sum = sum(tpm)) %>% 
  arrange(desc(tpm_sum)) -> hpfs_tpm_sum

## releveling the 'sample' column to plot stacked bar plot in descending order
hpfs_relevel <- as.character(hpfs_tpm_sum$sample)
hpfs_tpm_long$sample <- as.character(hpfs_tpm_long$sample)
hpfs_tpm_long$sample <- factor(hpfs_tpm_long$sample, levels = c(hpfs_relevel))

hpfs_tpm_long %>% 
  group_by(sample) %>% 
  summarise(tpm_sum = sum(tpm)) %>% 
  filter(tpm_sum < 1e3) %>% 
  arrange(desc(tpm_sum)) %>% 
  select(sample)  -> hpfs_tpm_low_sum
hpfs_low_samples <- hpfs_tpm_low_sum$sample
hpfs_tpm_long %>% 
  filter(sample %in% hpfs_low_samples) -> hpfs_tpm_long_low

##### ============================= PLOTTING ============================ #####

## plotting stacked barplot (david 2014)
clrs <- c("#93c47d", "#ffd966", "#e06666", "#5f89dd", "#c17a1b", "#999999")
ggplot(david2014_tpm_long, aes(fill = pathway, y = tpm, x = sample)) + 
  geom_bar(position = "stack", stat = "identity") +
  scale_y_continuous(expand = c(0,0)) +
  theme_bw() + theme(panel.grid.major = element_blank(),
                     panel.grid.minor = element_blank(),
                     axis.text.x = element_blank(),
                     axis.ticks.x = element_blank(),
                     text = element_text(size = 16)) +
  scale_fill_manual(values = c(clrs)) +
  ggtitle("David et al. 2014") + ylab("transcripts per million (TPM)")
ggsave("figures/supplementary/TPM_vs_sample_david2014.png", height = 5, width = 10, units = "in")
ggsave("figures/supplementary/TPM_vs_sample_david2014.svg", height = 5, width = 10, units = "in")

##### ALL SAMPLES #####

hpfs_expr_plot <- ggplot(hpfs_tpm_long, aes(fill = pathway, y = tpm, x = sample)) + 
  geom_bar(position = "stack", stat = "identity", width = 1) +
  # geom_hline(yintercept = 1e3, linetype = "dashed") +
  scale_y_continuous(expand = c(0,0), labels = comma) +
  theme_bw() + theme(panel.grid.major = element_blank(),
                     panel.grid.minor = element_blank(),
                     axis.text.x = element_blank(),
                     axis.ticks.x = element_blank(),
                     text = element_text(size = 16)) +
  scale_fill_manual(values = c(clrs)) +
  ggtitle("Health Professionals Follow-up Study (HPFS)") + ylab("transcripts per million (TPM)")
hpfs_expr_plot
ggsave("figures/supplementary/TPM_vs_sample_hpfs.png", height = 5, width = 10, units = "in")

inset_plot <- ggplot(hpfs_tpm_long_low, aes(fill = pathway, y = tpm, x = sample)) + 
  geom_bar(position = "stack", stat = "identity", width = 1) +
  scale_y_continuous(expand = c(0,0), limits = c(0,1e3),
                     labels = comma) +
  theme_bw() + theme(panel.grid.major = element_blank(),
                     panel.grid.minor = element_blank(),
                     axis.title.x = element_blank(),
                     axis.title.y = element_blank(),
                     axis.text.x = element_blank(),
                     axis.ticks.x = element_blank(),
                     legend.position = "none",
                     text = element_text(size = 16)) +
  scale_fill_manual(values = c(clrs)) 
inset_plot

hpfs_expr_plot + inset_element(
  inset_plot, 
  left = 0.06, 
  bottom = 0.1, 
  right = unit(1, 'npc') - unit(1, 'cm'), 
  top = unit(1, 'npc') - unit(1, 'cm')
)
ggsave("figures/supplementary/TPM_vs_sample_hpfs.png", height = 5, width = 10, units = "in")
ggsave("figures/supplementary/TPM_vs_sample_hpfs.svg", height = 5, width = 10, units = "in")
