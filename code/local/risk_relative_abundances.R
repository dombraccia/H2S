# ======================== risk_relative_abundances.R ======================= #
#' description: script for comparing the relative abundances of cysteine 
#' degrading and sulfite reducing bacteria in stool samples from cMD.
#' 
#' NOTE: it is assumed that all necessary packages and data are loaded into the 
#' R session where this script is being run.
# =========================================================================== #
##### =================================================================== #####
print("- calculating full column sums")
##### =================================================================== #####

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
### prism
cys_k2_prism_RA <- k2_prism_RA[cys_taxa, ]
cys_k2_prism_cs <- colSums(cys_k2_prism_RA)
dsrAB_k2_prism_RA <- k2_prism_RA[dsrAB_taxa, ]
dsrAB_k2_prism_cs <- colSums(dsrAB_k2_prism_RA)
### CIB
cys_k2_cib_RA <- k2_cib_RA[cys_taxa, ]
cys_k2_cib_cs <- colSums(cys_k2_cib_RA)
dsrAB_k2_cib_RA <- k2_cib_RA[dsrAB_taxa, ]
dsrAB_k2_cib_cs <- colSums(dsrAB_k2_cib_RA)

##### =================================================================== #####
print("- filtering pData to only include populations for plotting")
##### =================================================================== #####

### cMD

# filter risk samples
k2_pData %>%
  filter(study_condition == c("IBD", "CRC")) -> ibd_crc
unique(ibd_crc$dataset_name) -> ibd_crc_studies
kraken_risk_pData <- k2_pData[k2_pData$dataset_name %in% ibd_crc_studies, ]

#get colSums
cys_k2_cs_risk <- cys_k2_cs[names(cys_k2_cs) %in% rownames(kraken_risk_pData)]
dsrAB_k2_cs_risk <- dsrAB_k2_cs[names(dsrAB_k2_cs) %in% rownames(kraken_risk_pData)]

### hmp2
cys_k2_hmp2_cs_risk <- cys_k2_hmp2_cs[names(cys_k2_hmp2_cs) %in% k2_hmp2_pData$External.ID]
dsrAB_k2_hmp2_cs_risk <- dsrAB_k2_hmp2_cs[names(dsrAB_k2_hmp2_cs) %in% k2_hmp2_pData$External.ID]

### prism
cys_k2_prism_cs_risk <- cys_k2_prism_cs[names(cys_k2_prism_cs) %in% k2_prism_pData$Sample]
dsrAB_k2_prism_cs_risk <- dsrAB_k2_prism_cs[names(dsrAB_k2_prism_cs) %in% k2_prism_pData$Sample]

### cib
## filter risk samples
subset_k2_cib_metadata <- cbind(select(k2_cib_pData, Sample_SRA),
                                select(k2_cib_metadata, disease, Time, Sample.Name))
subset_k2_cib_metadata$Time[is.na(subset_k2_cib_metadata$Time)] <- 1
subset_k2_cib_metadata <- subset_k2_cib_metadata %>%
  filter(Time == c(1)) 
k2_cib_risk_metadata <- subset_k2_cib_metadata$disease
names(k2_cib_risk_metadata) <- subset_k2_cib_metadata$Sample_SRA 

## get colSums
cys_k2_cib_cs_risk <- cys_k2_cib_cs[names(cys_k2_cib_cs) %in% names(k2_cib_risk_metadata)]
dsrAB_k2_cib_cs_risk <- dsrAB_k2_cib_cs[names(dsrAB_k2_cib_cs) %in% names(k2_cib_risk_metadata)]


##### =================================================================== #####
print("- prepping dataframes for plotting")
##### =================================================================== #####

### cMD
cMD_study_cond <- as.character(kraken_risk_pData$study_condition)
names(cMD_study_cond) <- rownames(kraken_risk_pData)

cys_k2_df <- as.data.frame(cbind(cys_k2_cs_risk, cMD_study_cond = cMD_study_cond[names(cys_k2_cs_risk)]))
colnames(cys_k2_df) <- c("RA", "population")
cys_k2_df$RA <- as.numeric(cys_k2_df$RA)
cys_k2_df$func <- rep("cysteine", length(cys_k2_cs_risk))
cys_k2_df$population <- factor(cys_k2_df$population, levels = c("control", "IBD", "adenoma", "CRC"))

dsrAB_k2_df <- as.data.frame(cbind(dsrAB_k2_cs_risk, cMD_study_cond = cMD_study_cond[names(dsrAB_k2_cs_risk)]))
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

### prism
subset_k2_prism_pData <- select(k2_prism_pData, Sample, Diagnosis)
k2_prism_risk_pData <- subset_k2_prism_pData$Diagnosis
names(k2_prism_risk_pData) <- subset_k2_prism_pData$Sample

cys_k2_prism_df <- as.data.frame(cbind(cys_k2_prism_cs_risk, k2_prism_risk_pData = k2_prism_risk_pData[names(cys_k2_prism_cs_risk)]))
colnames(cys_k2_prism_df) <- c("RA", "diagnosis")
cys_k2_prism_df$RA <- as.numeric(cys_k2_prism_df$RA)
cys_k2_prism_df$func <- rep("cysteine", length(cys_k2_prism_cs_risk))
cys_k2_prism_df$diagnosis <- factor(cys_k2_prism_df$diagnosis, levels = c("HC", "UC", "CD"))

dsrAB_k2_prism_df <- as.data.frame(cbind(dsrAB_k2_prism_cs_risk, k2_prism_risk_pData = k2_prism_risk_pData[names(dsrAB_k2_prism_cs_risk)]))
colnames(dsrAB_k2_prism_df) <- c("RA", "diagnosis")
dsrAB_k2_prism_df$RA <- as.numeric(dsrAB_k2_prism_df$RA)
dsrAB_k2_prism_df$func <- rep("dsrAB", length(dsrAB_k2_prism_cs_risk))
dsrAB_k2_prism_df$diagnosis <- factor(dsrAB_k2_prism_df$diagnosis, levels = c("HC", "UC", "CD"))

### cib
cys_k2_cib_df <- as.data.frame(cbind(cys_k2_cib_cs_risk, k2_cib_risk_metadata ))
colnames(cys_k2_cib_df) <- c("RA", "diagnosis")
cys_k2_cib_df$RA <- as.numeric(cys_k2_cib_df$RA)
cys_k2_cib_df$func <- rep("cysteine", length(cys_k2_cib_cs_risk))
cys_k2_cib_df$diagnosis <- factor(cys_k2_cib_df$diagnosis, levels = c("Control", "Crohn"))

dsrAB_k2_cib_df <- as.data.frame(cbind(dsrAB_k2_cib_cs_risk, k2_cib_risk_metadata ))
colnames(dsrAB_k2_cib_df) <- c("RA", "diagnosis")
dsrAB_k2_cib_df$RA <- as.numeric(dsrAB_k2_cib_df$RA)
dsrAB_k2_cib_df$func <- rep("dsrAB", length(dsrAB_k2_cib_cs_risk))
dsrAB_k2_cib_df$diagnosis <- factor(dsrAB_k2_cib_df$diagnosis, levels = c("Control", "Crohn"))

##### =================================================================== #####
print("- calculating p-values between selected populations")
##### =================================================================== #####

#### cMD ####
cys_cMD_CRC_control <- wilcox.test(filter(cys_k2_df, population == "CRC")$RA, filter(cys_k2_df, population == "control")$RA)
cys_cMD_adenoma_control <- wilcox.test(filter(cys_k2_df, population == "control")$RA, filter(cys_k2_df, population == "adenoma")$RA)
cys_cMD_IBD_control <- wilcox.test(filter(cys_k2_df, population == "IBD")$RA,filter(cys_k2_df, population == "control")$RA)
dsrAB_cMD_CRC_control <- wilcox.test(filter(dsrAB_k2_df, population == "CRC")$RA,filter(dsrAB_k2_df, population == "control")$RA)
dsrAB_cMD_adenoma_control <- wilcox.test(filter(dsrAB_k2_df, population == "adenoma")$RA,filter(dsrAB_k2_df, population == "control")$RA)
#### ####

#### HMP2 ####
cys_hmp2_CD_nonIBD <- wilcox.test(filter(cys_k2_hmp2_df, diagnosis == "CD")$RA,filter(cys_k2_hmp2_df, diagnosis == "nonIBD")$RA)
cys_hmp2_UC_nonIBD <- wilcox.test(filter(cys_k2_hmp2_df, diagnosis == "UC")$RA,filter(cys_k2_hmp2_df, diagnosis == "nonIBD")$RA)
dsrAB_hmp2_CD_nonIBD <- wilcox.test(filter(dsrAB_k2_hmp2_df, diagnosis == "CD")$RA,filter(dsrAB_k2_hmp2_df, diagnosis == "nonIBD")$RA)
dsrAB_hmp2_UC_nonIBD <- wilcox.test(filter(dsrAB_k2_hmp2_df, diagnosis == "UC")$RA,filter(dsrAB_k2_hmp2_df, diagnosis == "nonIBD")$RA)
#### ####

#### prism ####
cys_prism_CD_HC <- wilcox.test(filter(cys_k2_prism_df, diagnosis == "CD")$RA,filter(cys_k2_prism_df, diagnosis == "HC")$RA)
cys_prism_UC_HC <- wilcox.test(filter(cys_k2_prism_df, diagnosis == "UC")$RA,filter(cys_k2_prism_df, diagnosis == "HC")$RA)
dsrAB_prism_CD_HC <- wilcox.test(filter(dsrAB_k2_prism_df, diagnosis == "CD")$RA,filter(dsrAB_k2_prism_df, diagnosis == "HC")$RA)
dsrAB_prism_UC_HC <- wilcox.test(filter(dsrAB_k2_prism_df, diagnosis == "UC")$RA,filter(dsrAB_k2_prism_df, diagnosis == "HC")$RA)
#### ####

#### cib ####
cys_cib_Crohn_control <- wilcox.test(filter(cys_k2_cib_df, diagnosis == "Crohn")$RA,filter(cys_k2_cib_df, diagnosis == "Control")$RA)
dsrAB_cib_Crohn_control <- wilcox.test(filter(dsrAB_k2_cib_df, diagnosis == "Crohn")$RA,filter(dsrAB_k2_cib_df, diagnosis == "Control")$RA)
#### ####

##### =================================================================== #####
print("- plotting risk cMD sample data")
##### =================================================================== #####

clrs <- c("#999999", "#E69F00", "#56B4E9", "#0072B2")
cys_k2_boxplot <- ggplot(cys_k2_df) + 
  geom_boxplot(aes(x = func, y = RA, fill = population)) +
  theme_bw() +
  theme(legend.position = "none",
        axis.text.x = element_blank(),
        axis.title.x = element_text(size = 8),
        axis.title.y = element_text(size = 8)) +
  scale_y_continuous(limits = c(0.001, 100),
                     breaks = c(0.01, 0.1, 1, 10, 100),
                     trans = "log2",
                     minor_breaks = NULL) +
  scale_fill_manual(values = clrs) +
  ylab("relative bacterial abundance, log2[%]") + xlab("PCDB")
# cys_k2_boxplot

dsrAB_k2_boxplot <- ggplot(dsrAB_k2_df) + 
  geom_boxplot(aes(x = func, y = RA, fill = population)) +
  theme_bw() +
  theme(legend.title = element_text(size = 8),
        legend.text = element_text(size = 6),
        legend.key.size = unit(0.3, "cm"),
        axis.title.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.text.y = element_blank(),
        axis.text.x = element_blank(),
        axis.title.x = element_text(size = 8)) +
  scale_y_continuous(limits = c(0.001, 100),
                     breaks = c(0.01, 0.1, 1, 10, 100),
                     trans = "log2", 
                     minor_breaks = NULL) +
  scale_fill_manual(values = clrs) +
  xlab("SRB")
dsrAB_k2_boxplot

##### =================================================================== #####
print("- arranging cMD plots")
##### =================================================================== #####
(cys_k2_boxplot | dsrAB_k2_boxplot) +
  plot_annotation(title = "curatedMetagneomicData",
                  theme = theme(plot.title = element_text(size = 22)))
# ggsave("../../figures/figure4/k2_cMD_RA_risk_populations.svg", height = 7, width = 7)
# ggsave("../../figures/figure4/k2_cMD_RA_risk_populations.png", height = 7, width = 7)
##### =================================================================== #####
print("- plotting hmp2 data")
##### =================================================================== #####

clrs <- c("#999999", "#F0E442", "#D55E00")
cys_k2_hmp2_boxplot <- ggplot(cys_k2_hmp2_df) + 
  geom_boxplot(aes(x = func, y = RA, fill = diagnosis)) +
  theme_bw() +
  theme(legend.position = "none",
        axis.text.x = element_blank(),
        panel.grid.minor.y = element_blank(),
        axis.title.x = element_text(size = 8),
        axis.title.y = element_text(size = 8)) +
  scale_y_continuous(limits = c(0.001, 100),
                     breaks = c(0.01, 0.1, 1, 10, 100),
                     trans = "log2") +
  scale_fill_manual(values = clrs) +
  ylab("relative bacterial abundance, log2[%]") + xlab("PCDB")
# cys_k2_hmp2_boxplot

dsrAB_k2_hmp2_boxplot <- ggplot(dsrAB_k2_hmp2_df) + 
  geom_boxplot(aes(x = func, y = RA, fill = diagnosis)) +
  theme_bw() +
  theme(legend.title = element_text(size = 8),
        legend.text = element_text(size = 6),
        axis.title.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.text.y = element_blank(),
        axis.text.x = element_blank(),
        panel.grid.minor.y = element_blank(),
        axis.title.x = element_text(size = 8)) +
  scale_y_continuous(limits = c(0.001, 100),
                     breaks = c(0.01, 0.1, 1, 10, 100),
                     trans = "log2", 
                     minor_breaks = NULL) +
  scale_fill_manual(values = clrs) +
  xlab("SRB")
# dsrAB_k2_hmp2_boxplot
##### =================================================================== #####
print("- arranging hmp2 plots")
##### =================================================================== #####

(cys_k2_hmp2_boxplot | dsrAB_k2_hmp2_boxplot) +
  plot_annotation(title = "HMP2",
                  theme = theme(plot.title = element_text(size = 22)))
# ggsave("../../figures/figure4/k2_hmp2_RA_risk_populations.svg", height = 7, width = 7)
# ggsave("../../figures/figure4/k2_hmp2_RA_risk_populations.png", height = 7, width = 7)
##### =================================================================== #####
print("- plotting prism data")
##### =================================================================== #####

clrs <- c("#999999", "#F0E442", "#D55E00")
cys_k2_prism_boxplot <- ggplot(cys_k2_prism_df) + 
  geom_boxplot(aes(x = func, y = RA, fill = diagnosis)) +
  theme_bw() +
  theme(legend.position = "none",
        axis.text.x = element_blank(),
        panel.grid.minor.y = element_blank(),
        axis.title.x = element_text(size = 8),
        axis.title.y = element_text(size = 8)) +
  scale_y_continuous(limits = c(0.001, 100),
                     breaks = c(0.01, 0.1, 1, 10, 100),
                     trans = "log2") +
  scale_fill_manual(values = clrs) +
  ylab("relative bacterial abundance, log2[%]") + xlab("PCDB")
# cys_k2_prism_boxplot

dsrAB_k2_prism_boxplot <- ggplot(dsrAB_k2_prism_df) + 
  geom_boxplot(aes(x = func, y = RA, fill = diagnosis)) +
  theme_bw() +
  theme(legend.title = element_text(size = 8),
        legend.text = element_text(size = 6),
        axis.title.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.text.y = element_blank(),
        axis.text.x = element_blank(),
        axis.title.x = element_text(size = 8)) +
  scale_y_continuous(limits = c(0.001, 100),
                     breaks = c(0.01, 0.1, 1, 10, 100),
                     trans = "log2", 
                     minor_breaks = NULL) +
  scale_fill_manual(values = clrs) +
  xlab("SRB")
# dsrAB_k2_prism_boxplot
##### =================================================================== #####
print("- arranging prism plots")
##### =================================================================== #####

(cys_k2_prism_boxplot | dsrAB_k2_prism_boxplot) +
  plot_annotation(title = "PRISM",
                  theme = theme(plot.title = element_text(size = 22)))
# ggsave("../../figures/figure4/k2_prism_RA_risk_populations.svg", height = 7, width = 7)
# ggsave("../../figures/figure4/k2_prism_RA_risk_populations.png", height = 7, width = 7)

##### =================================================================== #####
print("- plotting cib data")
##### =================================================================== #####

clrs <- c("#999999", "#D55E00")
cys_k2_cib_boxplot <- ggplot(cys_k2_cib_df) + 
  geom_boxplot(aes(x = func, y = RA, fill = diagnosis)) +
  theme_bw() +
  theme(legend.position = "none",
        axis.text.x = element_blank(),
        panel.grid.minor.y = element_blank(),
        axis.title.x = element_text(size = 8),
        axis.title.y = element_text(size = 8)) +
  scale_y_continuous(limits = c(0.001, 100),
                     breaks = c(0.01, 0.1, 1, 10, 100),
                     trans = "log2") +
  scale_fill_manual(values = clrs) +
  ylab("relative bacterial abundance, log2[%]") + xlab("PCDB")
# cys_k2_cib_boxplot

dsrAB_k2_cib_boxplot <- ggplot(dsrAB_k2_cib_df) + 
  geom_boxplot(aes(x = func, y = RA, fill = diagnosis)) +
  theme_bw() +
  theme(legend.title = element_text(size = 8),
        legend.text = element_text(size = 6),
        axis.title.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.text.y = element_blank(),
        axis.text.x = element_blank(),
        axis.title.x = element_text(size = 8)) +
  scale_y_continuous(limits = c(0.001, 100),
                     breaks = c(0.01, 0.1, 1, 10, 100),
                     trans = "log2", 
                     minor_breaks = NULL) +
  scale_fill_manual(values = clrs) +
  xlab("SRB")
# dsrAB_k2_cib_boxplot
##### =================================================================== #####
print("- arranging cib plots")
##### =================================================================== #####

(cys_k2_cib_boxplot | dsrAB_k2_cib_boxplot) +
  plot_annotation(title = "Pediatric Crohn's (Lewis et. al 2015)",
                  theme = theme(plot.title = element_text(size = 22)))
# ggsave("../../figures/figure4/k2_cib_RA_risk_populations.svg", height = 7, width = 7)
# ggsave("../../figures/figure4/k2_cib_RA_risk_populations.png", height = 7, width = 7)

