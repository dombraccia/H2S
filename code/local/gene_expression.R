# =========================== gene_expression.R ============================= #
#' description: processing and plotting of gene expression data for curated
#' cysteine degrading genes. Diamand (blastx mode) was used to map reads to 
#' custom database of specific genes examined in this paper.
#' 
#' NOTE: it is assumed that all necessary packages and data are loaded into the 
#' R session where this script is being run.
# =========================================================================== #



#### ======================== loading in data =========================== #####
print("- loading and processing data")

## count data
counts_raw <- read.csv("../../results/diamond_processed/counts.tsv", 
         sep = "\t", header = TRUE)
rnames <- counts_raw$Run
counts <- t(counts_raw[, 2:8])
colnames(counts) <- rnames

## metadata
metadata <- read.csv("../../data/david-2014/metadata.txt", header = TRUE)

#### ==================================================================== #####

#### ======================== prep for plotting ========================= #####
print("- prepping data for plotting")

## subsetting metadata based on diet and time
metadata$Time <- as.numeric(gsub("Day ", "", metadata$Time))
metadata %>% filter(Time < 0) -> baseline_metadata
metadata %>% filter(Time > 0, diet == "plant-based") -> plant_metadata
metadata %>% filter(Time > 0, diet == "animal-based") -> animal_metadata

## prepping baseline data
baseline_counts <- counts[, colnames(counts) %in% baseline_metadata$Run]
baseline_df <- melt(baseline_counts)
colnames(baseline_df) <- c("genes", "Run", "counts")
# baseline_df$counts <- log2(as.numeric(baseline_df$counts))
# baseline_df$counts[!is.finite(baseline_df$counts)] <- 0
baseline_df$genes <- factor(baseline_df$genes, 
                            levels =  c("MST", "CYD", "CBS", "CSE", 
                                        "MGL", "dsrA", "dsrB"))

plant_counts <- counts[, colnames(counts) %in% plant_metadata$Run]
plant_df <- melt(plant_counts)
colnames(plant_df) <- c("genes", "Run", "counts")
# plant_df$counts <- log2(as.numeric(plant_df$counts))
# plant_df$counts[!is.finite(plant_df$counts)] <- 0
plant_df$genes <- factor(plant_df$genes, 
                            levels =  c("MST", "CYD", "CBS", "CSE", 
                                        "MGL", "dsrA", "dsrB"))

animal_counts <- counts[, colnames(counts) %in% animal_metadata$Run]
animal_df <- melt(animal_counts)
colnames(animal_df) <- c("genes", "Run", "counts")
# animal_df$counts <- log2(as.numeric(animal_df$counts))
# animal_df$counts[!is.finite(animal_df$counts)] <- 0
animal_df$genes <- factor(animal_df$genes, 
                            levels =  c("MST", "CYD", "CBS", "CSE", 
                                        "MGL", "dsrA", "dsrB"))

#### ==================================================================== #####

#### ============================ plotting ============================== #####
print("- plotting data")

## baseline boxplot
baseline_boxplot <- ggplot(data = baseline_df, aes(x = genes, fill = genes, y = counts)) +
  geom_boxplot(outlier.shape = NA) + geom_jitter(alpha = 0.8, size = 0.4, width = 0.2) + 
  theme_bw() +
  theme(legend.position = "none",
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.title.x = element_text(size = 16),
        axis.title.y = element_text(size = 16)) +
  scale_y_continuous(limits = c(1, 6000), 
                     breaks = c(1, 10, 100, 1000, 5000),
                     trans = "log2",
                     minor_breaks = NULL) +
  scale_fill_manual(values = c("#d3411a", "#d85b1d", "#dd9e2f", "#e2c137", 
                               "#ebf268", "#1f97c1", "#1e5bc6")) +
  xlab("baseline") + ylab("Reads mapped per gene (raw counts)")

## plant boxplot
plant_boxplot <- ggplot(data = plant_df, aes(x = genes, fill = genes, y = counts)) +
  geom_boxplot(outlier.shape = NA) + geom_jitter(alpha = 0.8, size = 0.4, width = 0.2) +
  theme_bw() +
  theme(legend.position = "none",
        axis.title.y = element_blank(),
        axis.text.x = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks.x = element_blank(),
        axis.ticks.y = element_blank(),
        axis.title.x = element_text(size = 16)) +
  scale_y_continuous(limits = c(1, 6000), 
                     breaks = c(1, 10, 100, 1000, 5000),
                     trans = "log2",
                     minor_breaks = NULL) +
  scale_fill_manual(values = c("#d3411a", "#d85b1d", "#dd9e2f", "#e2c137", 
                               "#ebf268", "#1f97c1", "#1e5bc6")) +
  xlab("plant-based")

## animal boxplot
animal_boxplot <- ggplot(data = animal_df, aes(x = genes, fill = genes, y = counts)) +
  geom_boxplot(outlier.shape = NA) + geom_jitter(alpha = 0.8, size = 0.4, width = 0.2) +
  theme_bw() +
  theme(axis.title.y = element_blank(),
        axis.text.y = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.ticks.y = element_blank(),
        axis.title.x = element_text(size = 16)) +
  scale_y_continuous(limits = c(1, 6000), 
                     breaks = c(1, 10, 100, 1000, 5000),
                     trans = "log2",
                     minor_breaks = NULL) +
  scale_fill_manual(values = c("#d3411a", "#d85b1d", "#dd9e2f", "#e2c137", 
                               "#ebf268", "#1f97c1", "#1e5bc6")) +
  xlab("animal-based")

baseline_boxplot | plant_boxplot | animal_boxplot
ggsave("../../figures/supplementary/david2014_RNAseq_diet.svg", height = 7, width = 10)
ggsave("../../figures/supplementary/david2014_RNAseq_diet.png", height = 7, width = 10)
#### ==================================================================== #####

