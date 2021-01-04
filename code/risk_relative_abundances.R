# ======================== risk_relative_abundances.R ======================= #
#' description: script for comparing the relative abundances of cysteine 
#' degrading and sulfite reducing bacteria in stool samples from cMD.
#' 
#' NOTE: it is assumed that all necessary packages and data are loaded into the 
#' R session where this script is being run.
# =========================================================================== #

print("- calculating full column sums")
### cMD
cys_k2_RA <- k2_RA[cys_taxa, ]
cys_k2_cs <- colSums(cys_k2_RA)
dsrAB_k2_RA <- k2_RA[dsrAB_taxa, ]
dsrAB_k2_cs <- colSums(dsrAB_k2_RA)
### hmp2
cys_k2_hmp2_RA <- k2_hmp2_RA[cys_taxa, ]
cys_k2_hmp2_cs <- colSums(cys_k2_hmp2_RA)
dsrAB_k2_hmp2_RA <- k2_hmp2_RA[dsrAB_taxa, ]
dsrAB_k2_hmp2_cs <- colSums(dsrAB_k2_hmp2_RA)

print("- filtering pData to only include populations for plotting")
k2_pData %>%
  filter(study_condition == c("IBD", "CRC")) -> ibd_crc
unique(ibd_crc$dataset_name) -> ibd_crc_studies
kraken_risk_pData <- k2_pData[k2_pData$dataset_name %in% ibd_crc_studies, ]
### cMD
cys_k2_cs_risk <- cys_k2_cs[names(cys_k2_cs) %in% rownames(kraken_risk_pData)]
dsrAB_k2_cs_risk <- dsrAB_k2_cs[names(dsrAB_k2_cs) %in% rownames(kraken_risk_pData)]
### HMP2
cys_k2_hmp2_cs_risk <- cys_k2_hmp2_cs[names(cys_k2_hmp2_cs) %in% k2_hmp2_pData$External.ID]
# tmp <- k2_hmp2_pData[k2_hmp2_pData$External.ID.mod %in% names(cys_k2_hmp2_cs),]
dsrAB_k2_hmp2_cs_risk <- dsrAB_k2_hmp2_cs[names(dsrAB_k2_hmp2_cs) %in% k2_hmp2_pData$External.ID]

print("- prepping dataframes for plotting")
### cMD
hmp2_study_cond <- as.character(kraken_risk_pData$study_condition)
names(hmp2_study_cond) <- rownames(kraken_risk_pData)

cys_k2_df <- as.data.frame(cbind(cys_k2_cs_risk, hmp2_study_cond = hmp2_study_cond[names(cys_k2_cs_risk)]))
colnames(cys_k2_df) <- c("RA", "population")
cys_k2_df$RA <- as.numeric(cys_k2_df$RA)
cys_k2_df$func <- rep("cysteine", length(cys_k2_cs_risk))
cys_k2_df$population <- factor(cys_k2_df$population, levels = c("control", "IBD", "adenoma", "CRC"))

dsrAB_k2_df <- as.data.frame(cbind(dsrAB_k2_cs_risk, hmp2_study_cond = hmp2_study_cond[names(dsrAB_k2_cs_risk)]))
colnames(dsrAB_k2_df) <- c("RA", "population")
dsrAB_k2_df$RA <- as.numeric(dsrAB_k2_df$RA)
dsrAB_k2_df$func <- rep("dsrAB", length(dsrAB_k2_cs_risk))
dsrAB_k2_df$population <- factor(dsrAB_k2_df$population, levels = c("control", "IBD", "adenoma", "CRC"))

### HMP2
subset_k2_pData <- select(k2_hmp2_pData, External.ID, diagnosis)
k2_risk_pData <- subset_k2_pData$diagnosis
names(k2_risk_pData) <- subset_k2_pData$External.ID

cys_k2_hmp2_df <- as.data.frame(cbind(cys_k2_hmp2_cs_risk, k2_risk_pData = k2_risk_pData[names(cys_k2_hmp2_cs_risk)]))
colnames(cys_k2_hmp2_df) <- c("RA", "diagnosis")
cys_k2_hmp2_df$RA <- as.numeric(cys_k2_hmp2_df$RA)
cys_k2_hmp2_df$func <- rep("cysteine", length(cys_k2_hmp2_cs_risk))
cys_k2_hmp2_df$diagnosis <- factor(cys_k2_hmp2_df$diagnosis, levels = c("nonIBD", "UC", "CD"))
# colnames(cys_k2_hmp2_df) <- c("RA", "diagnosis", "func")

dsrAB_k2_hmp2_df <- as.data.frame(cbind(dsrAB_k2_hmp2_cs_risk, k2_risk_pData = k2_risk_pData[names(dsrAB_k2_hmp2_cs_risk)]))
colnames(dsrAB_k2_hmp2_df) <- c("RA", "diagnosis")
dsrAB_k2_hmp2_df$RA <- as.numeric(dsrAB_k2_hmp2_df$RA)
dsrAB_k2_hmp2_df$func <- rep("dsrAB", length(dsrAB_k2_hmp2_cs_risk))
dsrAB_k2_hmp2_df$diagnosis <- factor(dsrAB_k2_hmp2_df$diagnosis, levels = c("nonIBD", "UC", "CD"))
# colnames(dsrAB_k2_hmp2_df) <- c("RA", "diagnosis", "func")

print("- plotting risk cMD sample data")
clrs <- c("#999999", "#E69F00", "#56B4E9", "#0072B2")
cys_k2_boxplot <- ggplot(cys_k2_df) + 
  geom_boxplot(aes(x = func, y = RA, fill = population)) +
  theme_bw() +
  theme(legend.position = "none",
        axis.text.x = element_blank()) +
  scale_y_continuous(limits = c(0.001, 100),
                     breaks = c(0.01, 0.1, 1, 10, 100),
                     trans = "log2",
                     minor_breaks = NULL) +
  scale_fill_manual(values = clrs) +
  ylab("relative bacterial abundance, log2[%]") + xlab("cysteine_degrading")
# cys_k2_boxplot

dsrAB_k2_boxplot <- ggplot(dsrAB_k2_df) + 
  geom_boxplot(aes(x = func, y = RA, fill = population)) +
  theme_bw() +
  theme(axis.title.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.text.y = element_blank(),
        axis.text.x = element_blank()) +
  scale_y_continuous(limits = c(0.001, 100),
                     breaks = c(0.01, 0.1, 1, 10, 100),
                     trans = "log2", 
                     minor_breaks = NULL) +
  scale_fill_manual(values = clrs) +
  xlab("sulfite_reducing")
# dsrAB_k2_boxplot

print("- arranging cMD plots")
cys_k2_boxplot | dsrAB_k2_boxplot
ggsave("../figures/figure3/k2_cMD_RA_risk_populations.svg", height = 7, width = 7)
ggsave("../figures/figure3/k2_cMD_RA_risk_populations.png", height = 7, width = 7)

# =========================================================================== #

print("- plotting hmp2 data")
clrs <- c("#999999", "#F0E442", "#D55E00")
cys_k2_hmp2_boxplot <- ggplot(cys_k2_hmp2_df) + 
  geom_boxplot(aes(x = func, y = RA, fill = diagnosis)) +
  theme_bw() +
  theme(legend.position = "none",
        axis.text.x = element_blank()) +
  scale_y_continuous(limits = c(0.001, 100),
                     breaks = c(0.01, 0.1, 1, 10, 100),
                     trans = "log2") +
  scale_fill_manual(values = clrs) +
  ylab("relative bacterial abundance, log2[%]") + xlab("cysteine_degrading")
# cys_k2_hmp2_boxplot

dsrAB_k2_hmp2_boxplot <- ggplot(dsrAB_k2_hmp2_df) + 
  geom_boxplot(aes(x = func, y = RA, fill = diagnosis)) +
  theme_bw() +
  theme(axis.title.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.text.y = element_blank(),
        axis.text.x = element_blank()) +
  scale_y_continuous(limits = c(0.001, 100),
                     breaks = c(0.01, 0.1, 1, 10, 100),
                     trans = "log2", 
                     minor_breaks = NULL) +
  scale_fill_manual(values = clrs) +
  xlab("sulfite_reducing")
# dsrAB_k2_hmp2_boxplot

# =========================================================================== #

print("- arranging hmp2 plots")
cys_k2_hmp2_boxplot | dsrAB_k2_hmp2_boxplot
ggsave("../figures/figure3/k2_hmp2_RA_risk_populations.svg", height = 7, width = 7)
ggsave("../figures/figure3/k2_hmp2_RA_risk_populations.png", height = 7, width = 7)
