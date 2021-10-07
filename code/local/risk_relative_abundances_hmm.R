# ======================== risk_relative_abundances.R ======================= #
#' description: script for comparing the relative abundances of cysteine 
#' degrading and sulfite reducing bacteria in stool samples from cMD.
#' 
#' NOTE: it is assumed that all necessary packages and data are loaded into the 
#' R session where this script is being run.
# =========================================================================== #

##### =================================================================== #####
print("- calculating full column sums")

### cMD
primary_k2_RA <- k2_RA[primary_taxa$taxa_ID, ]
primary_k2_cs <- colSums(primary_k2_RA)
secondary_k2_RA <- k2_RA[secondary_taxa$taxa_ID, ]
secondary_k2_cs <- colSums(secondary_k2_RA)
erroneous_k2_RA <- k2_RA[erroneous_taxa$taxa_ID, ]
erroneous_k2_cs <- colSums(erroneous_k2_RA)
dsrAB_k2_RA <- k2_RA[dsrAB_taxa, ]
dsrAB_k2_cs <- colSums(dsrAB_k2_RA)

### hmp2
primary_k2_hmp2_RA <- k2_hmp2_RA[primary_taxa$taxa_ID, ]
primary_k2_hmp2_cs <- colSums(primary_k2_hmp2_RA)
secondary_k2_hmp2_RA <- k2_hmp2_RA[secondary_taxa$taxa_ID, ]
secondary_k2_hmp2_cs <- colSums(secondary_k2_hmp2_RA)
erroneous_k2_hmp2_RA <- k2_hmp2_RA[erroneous_taxa$taxa_ID, ]
erroneous_k2_hmp2_cs <- colSums(erroneous_k2_hmp2_RA)
dsrAB_k2_hmp2_RA <- k2_hmp2_RA[dsrAB_taxa, ]
dsrAB_k2_hmp2_cs <- colSums(dsrAB_k2_hmp2_RA)

### prism
primary_k2_prism_RA <- k2_prism_RA[primary_taxa$taxa_ID, ]
primary_k2_prism_cs <- colSums(primary_k2_prism_RA)
secondary_k2_prism_RA <- k2_prism_RA[secondary_taxa$taxa_ID, ]
secondary_k2_prism_cs <- colSums(secondary_k2_prism_RA)
erroneous_k2_prism_RA <- k2_prism_RA[erroneous_taxa$taxa_ID, ]
erroneous_k2_prism_cs <- colSums(erroneous_k2_prism_RA)
dsrAB_k2_prism_RA <- k2_prism_RA[dsrAB_taxa, ]
dsrAB_k2_prism_cs <- colSums(dsrAB_k2_prism_RA)

### CIB
primary_k2_cib_RA <- k2_cib_RA[primary_taxa$taxa_ID, ]
primary_k2_cib_cs <- colSums(primary_k2_cib_RA)
secondary_k2_cib_RA <- k2_cib_RA[secondary_taxa$taxa_ID, ]
secondary_k2_cib_cs <- colSums(secondary_k2_cib_RA)
erroneous_k2_cib_RA <- k2_cib_RA[erroneous_taxa$taxa_ID, ]
erroneous_k2_cib_cs <- colSums(erroneous_k2_cib_RA)
dsrAB_k2_cib_RA <- k2_cib_RA[dsrAB_taxa, ]
dsrAB_k2_cib_cs <- colSums(dsrAB_k2_cib_RA)

##### =================================================================== #####
print("- filtering pData to only include populations for plotting")

### cMD

# filter risk samples
k2_pData %>%
  filter(study_condition == c("IBD", "CRC")) -> ibd_crc
unique(ibd_crc$dataset_name) -> ibd_crc_studies
kraken_risk_pData <- k2_pData[k2_pData$dataset_name %in% ibd_crc_studies, ]

#get colSums
primary_k2_cs_risk <- primary_k2_cs[names(primary_k2_cs) %in% rownames(kraken_risk_pData)]
secondary_k2_cs_risk <- secondary_k2_cs[names(secondary_k2_cs) %in% rownames(kraken_risk_pData)]
erroneous_k2_cs_risk <- erroneous_k2_cs[names(erroneous_k2_cs) %in% rownames(kraken_risk_pData)]
dsrAB_k2_cs_risk <- dsrAB_k2_cs[names(dsrAB_k2_cs) %in% rownames(kraken_risk_pData)]

### hmp2
primary_k2_hmp2_cs_risk <- primary_k2_hmp2_cs[names(primary_k2_hmp2_cs) %in% k2_hmp2_pData$External.ID]
secondary_k2_hmp2_cs_risk <- secondary_k2_hmp2_cs[names(secondary_k2_hmp2_cs) %in% k2_hmp2_pData$External.ID]
erroneous_k2_hmp2_cs_risk <- erroneous_k2_hmp2_cs[names(erroneous_k2_hmp2_cs) %in% k2_hmp2_pData$External.ID]
dsrAB_k2_hmp2_cs_risk <- dsrAB_k2_hmp2_cs[names(dsrAB_k2_hmp2_cs) %in% k2_hmp2_pData$External.ID]

### prism
primary_k2_prism_cs_risk <- primary_k2_prism_cs[names(primary_k2_prism_cs) %in% k2_prism_pData$Sample]
secondary_k2_prism_cs_risk <- secondary_k2_prism_cs[names(secondary_k2_prism_cs) %in% k2_prism_pData$Sample]
erroneous_k2_prism_cs_risk <- erroneous_k2_prism_cs[names(erroneous_k2_prism_cs) %in% k2_prism_pData$Sample]
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
primary_k2_cib_cs_risk <- primary_k2_cib_cs[names(primary_k2_cib_cs) %in% k2_cib_pData$Sample_SRA]
secondary_k2_cib_cs_risk <- secondary_k2_cib_cs[names(secondary_k2_cib_cs) %in% k2_cib_pData$Sample_SRA]
erroneous_k2_cib_cs_risk <- erroneous_k2_cib_cs[names(erroneous_k2_cib_cs) %in% k2_cib_pData$Sample_SRA]
dsrAB_k2_cib_cs_risk <- dsrAB_k2_cib_cs[names(dsrAB_k2_cib_cs) %in% names(k2_cib_risk_metadata)]

##### =================================================================== #####
print("- prepping dataframes for plotting")
##### =================================================================== #####

### cMD
cMD_study_cond <- as.character(kraken_risk_pData$study_condition)
names(cMD_study_cond) <- rownames(kraken_risk_pData)

primary_k2_df <- as.data.frame(cbind(primary_k2_cs_risk, cMD_study_cond = cMD_study_cond[names(primary_k2_cs_risk)]))
colnames(primary_k2_df) <- c("RA", "population")
primary_k2_df$RA <- as.numeric(primary_k2_df$RA)
primary_k2_df$func <- rep("cysteine", length(primary_k2_cs_risk))
primary_k2_df$population <- factor(primary_k2_df$population, levels = c("control", "IBD", "adenoma", "CRC"))

secondary_k2_df <- as.data.frame(cbind(secondary_k2_cs_risk, cMD_study_cond = cMD_study_cond[names(secondary_k2_cs_risk)]))
colnames(secondary_k2_df) <- c("RA", "population")
secondary_k2_df$RA <- as.numeric(secondary_k2_df$RA)
secondary_k2_df$func <- rep("cysteine", length(secondary_k2_cs_risk))
secondary_k2_df$population <- factor(secondary_k2_df$population, levels = c("control", "IBD", "adenoma", "CRC"))

erroneous_k2_df <- as.data.frame(cbind(erroneous_k2_cs_risk, cMD_study_cond = cMD_study_cond[names(erroneous_k2_cs_risk)]))
colnames(erroneous_k2_df) <- c("RA", "population")
erroneous_k2_df$RA <- as.numeric(erroneous_k2_df$RA)
erroneous_k2_df$func <- rep("cysteine", length(erroneous_k2_cs_risk))
erroneous_k2_df$population <- factor(erroneous_k2_df$population, levels = c("control", "IBD", "adenoma", "CRC"))

dsrAB_k2_df <- as.data.frame(cbind(dsrAB_k2_cs_risk, cMD_study_cond = cMD_study_cond[names(dsrAB_k2_cs_risk)]))
colnames(dsrAB_k2_df) <- c("RA", "population")
dsrAB_k2_df$RA <- as.numeric(dsrAB_k2_df$RA)
dsrAB_k2_df$func <- rep("dsrAB", length(dsrAB_k2_cs_risk))
dsrAB_k2_df$population <- factor(dsrAB_k2_df$population, levels = c("control", "IBD", "adenoma", "CRC"))

#### HMP2 ####
subset_k2_pData <- select(k2_hmp2_pData, External.ID, diagnosis)
k2_risk_pData <- subset_k2_pData$diagnosis
names(k2_risk_pData) <- subset_k2_pData$External.ID

primary_k2_hmp2_df <- as.data.frame(cbind(primary_k2_hmp2_cs_risk, k2_risk_pData = k2_risk_pData[names(primary_k2_hmp2_cs_risk)]))
colnames(primary_k2_hmp2_df) <- c("RA", "diagnosis")
primary_k2_hmp2_df$RA <- as.numeric(primary_k2_hmp2_df$RA)
primary_k2_hmp2_df$func <- rep("cysteine", length(primary_k2_hmp2_cs_risk))
primary_k2_hmp2_df$diagnosis <- factor(primary_k2_hmp2_df$diagnosis, levels = c("nonIBD", "UC", "CD"))
# colnames(primary_k2_hmp2_df) <- c("RA", "diagnosis", "func")

secondary_k2_hmp2_df <- as.data.frame(cbind(secondary_k2_hmp2_cs_risk, k2_risk_pData = k2_risk_pData[names(secondary_k2_hmp2_cs_risk)]))
colnames(secondary_k2_hmp2_df) <- c("RA", "diagnosis")
secondary_k2_hmp2_df$RA <- as.numeric(secondary_k2_hmp2_df$RA)
secondary_k2_hmp2_df$func <- rep("cysteine", length(secondary_k2_hmp2_cs_risk))
secondary_k2_hmp2_df$diagnosis <- factor(secondary_k2_hmp2_df$diagnosis, levels = c("nonIBD", "UC", "CD"))
# colnames(secondary_k2_hmp2_df) <- c("RA", "diagnosis", "func")

erroneous_k2_hmp2_df <- as.data.frame(cbind(erroneous_k2_hmp2_cs_risk, k2_risk_pData = k2_risk_pData[names(erroneous_k2_hmp2_cs_risk)]))
colnames(erroneous_k2_hmp2_df) <- c("RA", "diagnosis")
erroneous_k2_hmp2_df$RA <- as.numeric(erroneous_k2_hmp2_df$RA)
erroneous_k2_hmp2_df$func <- rep("cysteine", length(erroneous_k2_hmp2_cs_risk))
erroneous_k2_hmp2_df$diagnosis <- factor(erroneous_k2_hmp2_df$diagnosis, levels = c("nonIBD", "UC", "CD"))
# colnames(erroneous_k2_hmp2_df) <- c("RA", "diagnosis", "func")

dsrAB_k2_hmp2_df <- as.data.frame(cbind(dsrAB_k2_hmp2_cs_risk, k2_risk_pData = k2_risk_pData[names(dsrAB_k2_hmp2_cs_risk)]))
colnames(dsrAB_k2_hmp2_df) <- c("RA", "diagnosis")
dsrAB_k2_hmp2_df$RA <- as.numeric(dsrAB_k2_hmp2_df$RA)
dsrAB_k2_hmp2_df$func <- rep("dsrAB", length(dsrAB_k2_hmp2_cs_risk))
dsrAB_k2_hmp2_df$diagnosis <- factor(dsrAB_k2_hmp2_df$diagnosis, levels = c("nonIBD", "UC", "CD"))
# colnames(dsrAB_k2_hmp2_df) <- c("RA", "diagnosis", "func")

#### PRISM ####
subset_k2_prism_pData <- select(k2_prism_pData, Sample, Diagnosis)
k2_prism_risk_pData <- subset_k2_prism_pData$Diagnosis
names(k2_prism_risk_pData) <- subset_k2_prism_pData$Sample

primary_k2_prism_df <- as.data.frame(cbind(primary_k2_prism_cs_risk, k2_prism_risk_pData = k2_prism_risk_pData[names(primary_k2_prism_cs_risk)]))
colnames(primary_k2_prism_df) <- c("RA", "diagnosis")
primary_k2_prism_df$RA <- as.numeric(primary_k2_prism_df$RA)
primary_k2_prism_df$func <- rep("cysteine", length(primary_k2_prism_cs_risk))
primary_k2_prism_df$diagnosis <- factor(primary_k2_prism_df$diagnosis, levels = c("HC", "UC", "CD"))

secondary_k2_prism_df <- as.data.frame(cbind(secondary_k2_prism_cs_risk, k2_prism_risk_pData = k2_prism_risk_pData[names(secondary_k2_prism_cs_risk)]))
colnames(secondary_k2_prism_df) <- c("RA", "diagnosis")
secondary_k2_prism_df$RA <- as.numeric(secondary_k2_prism_df$RA)
secondary_k2_prism_df$func <- rep("cysteine", length(secondary_k2_prism_cs_risk))
secondary_k2_prism_df$diagnosis <- factor(secondary_k2_prism_df$diagnosis, levels = c("HC", "UC", "CD"))

erroneous_k2_prism_df <- as.data.frame(cbind(erroneous_k2_prism_cs_risk, k2_prism_risk_pData = k2_prism_risk_pData[names(erroneous_k2_prism_cs_risk)]))
colnames(erroneous_k2_prism_df) <- c("RA", "diagnosis")
erroneous_k2_prism_df$RA <- as.numeric(erroneous_k2_prism_df$RA)
erroneous_k2_prism_df$func <- rep("cysteine", length(erroneous_k2_prism_cs_risk))
erroneous_k2_prism_df$diagnosis <- factor(erroneous_k2_prism_df$diagnosis, levels = c("HC", "UC", "CD"))

dsrAB_k2_prism_df <- as.data.frame(cbind(dsrAB_k2_prism_cs_risk, k2_prism_risk_pData = k2_prism_risk_pData[names(dsrAB_k2_prism_cs_risk)]))
colnames(dsrAB_k2_prism_df) <- c("RA", "diagnosis")
dsrAB_k2_prism_df$RA <- as.numeric(dsrAB_k2_prism_df$RA)
dsrAB_k2_prism_df$func <- rep("dsrAB", length(dsrAB_k2_prism_cs_risk))
dsrAB_k2_prism_df$diagnosis <- factor(dsrAB_k2_prism_df$diagnosis, levels = c("HC", "UC", "CD"))

#### CIB ####
primary_k2_cib_df <- as.data.frame(cbind(primary_k2_cib_cs_risk, k2_cib_risk_metadata ))
colnames(primary_k2_cib_df) <- c("RA", "diagnosis")
primary_k2_cib_df$RA <- as.numeric(primary_k2_cib_df$RA)
primary_k2_cib_df$func <- rep("primary", length(primary_k2_cib_cs_risk))
primary_k2_cib_df$diagnosis <- factor(primary_k2_cib_df$diagnosis, levels = c("Control", "Crohn"))

secondary_k2_cib_df <- as.data.frame(cbind(secondary_k2_cib_cs_risk, k2_cib_risk_metadata ))
colnames(secondary_k2_cib_df) <- c("RA", "diagnosis")
secondary_k2_cib_df$RA <- as.numeric(secondary_k2_cib_df$RA)
secondary_k2_cib_df$func <- rep("secondary", length(secondary_k2_cib_cs_risk))
secondary_k2_cib_df$diagnosis <- factor(secondary_k2_cib_df$diagnosis, levels = c("Control", "Crohn"))

erroneous_k2_cib_df <- as.data.frame(cbind(erroneous_k2_cib_cs_risk, k2_cib_risk_metadata ))
colnames(erroneous_k2_cib_df) <- c("RA", "diagnosis")
erroneous_k2_cib_df$RA <- as.numeric(erroneous_k2_cib_df$RA)
erroneous_k2_cib_df$func <- rep("erroneous", length(erroneous_k2_cib_cs_risk))
erroneous_k2_cib_df$diagnosis <- factor(erroneous_k2_cib_df$diagnosis, levels = c("Control", "Crohn"))

dsrAB_k2_cib_df <- as.data.frame(cbind(dsrAB_k2_cib_cs_risk, k2_cib_risk_metadata ))
colnames(dsrAB_k2_cib_df) <- c("RA", "diagnosis")
dsrAB_k2_cib_df$RA <- as.numeric(dsrAB_k2_cib_df$RA)
dsrAB_k2_cib_df$func <- rep("dsrAB", length(dsrAB_k2_cib_cs_risk))
dsrAB_k2_cib_df$diagnosis <- factor(dsrAB_k2_cib_df$diagnosis, levels = c("Control", "Crohn"))

##### =================================================================== #####
print("- [NOT RUN] calculating p-values between selected populations")

#### cMD ####
primary_cMD_CRC_control <- wilcox.test(filter(primary_k2_df, population == "CRC")$RA, 
                                       filter(primary_k2_df, population == "control")$RA)
primary_cMD_adenoma_control <- wilcox.test(filter(primary_k2_df, population == "adenoma")$RA, 
                                           filter(primary_k2_df, population == "control")$RA)
primary_cMD_adenoma_CRC <- wilcox.test(filter(primary_k2_df, population == "adenoma")$RA, 
                                       filter(primary_k2_df, population == "CRC")$RA)
primary_cMD_IBD_control <- wilcox.test(filter(primary_k2_df, population == "IBD")$RA, 
                                       filter(primary_k2_df, population == "control")$RA)
secondary_cMD_CRC_control <- wilcox.test(filter(secondary_k2_df, population == "CRC")$RA, 
                                         filter(secondary_k2_df, population == "control")$RA)
secondary_cMD_CRC_adenoma <- wilcox.test(filter(secondary_k2_df, population == "CRC")$RA, 
                                         filter(secondary_k2_df, population == "adenoma")$RA)
secondary_cMD_IBD_control <- wilcox.test(filter(secondary_k2_df, population == "IBD")$RA, 
                                         filter(secondary_k2_df, population == "control")$RA)
erroneous_cMD_CRC_control <- wilcox.test(filter(erroneous_k2_df, population == "CRC")$RA, 
                                         filter(erroneous_k2_df, population == "control")$RA)
erroneous_cMD_CRC_adenoma <- wilcox.test(filter(erroneous_k2_df, population == "CRC")$RA, 
                                         filter(erroneous_k2_df, population == "adenoma")$RA)
erroneous_cMD_IBD_adenoma <- wilcox.test(filter(erroneous_k2_df, population == "IBD")$RA, 
                                         filter(erroneous_k2_df, population == "control")$RA)
dsrAB_cMD_CRC_control <- wilcox.test(filter(dsrAB_k2_df, population == "CRC")$RA,
                                     filter(dsrAB_k2_df, population == "control")$RA)
dsrAB_cMD_adenoma_control <- wilcox.test(filter(dsrAB_k2_df, population == "adenoma")$RA,
                                         filter(dsrAB_k2_df, population == "control")$RA)
dsrAB_cMD_CRC_adenoma <- wilcox.test(filter(dsrAB_k2_df, population == "CRC")$RA,
                                     filter(dsrAB_k2_df, population == "adenoma")$RA)
dsrAB_cMD_IBD_control <- wilcox.test(filter(dsrAB_k2_df, population == "IBD")$RA,
                                     filter(dsrAB_k2_df, population == "control")$RA)


#### HMP2 ####
primary_hmp2_CD_nonIBD <- wilcox.test(filter(primary_k2_hmp2_df, diagnosis == "CD")$RA,filter(primary_k2_hmp2_df, diagnosis == "nonIBD")$RA)
primary_hmp2_UC_nonIBD <- wilcox.test(filter(primary_k2_hmp2_df, diagnosis == "UC")$RA,filter(primary_k2_hmp2_df, diagnosis == "nonIBD")$RA)
secondary_hmp2_CD_nonIBD <- wilcox.test(filter(secondary_k2_hmp2_df, diagnosis == "CD")$RA,filter(secondary_k2_hmp2_df, diagnosis == "nonIBD")$RA)
secondary_hmp2_UC_nonIBD <- wilcox.test(filter(secondary_k2_hmp2_df, diagnosis == "UC")$RA,filter(secondary_k2_hmp2_df, diagnosis == "nonIBD")$RA)
dsrAB_hmp2_CD_nonIBD <- wilcox.test(filter(dsrAB_k2_hmp2_df, diagnosis == "CD")$RA,filter(dsrAB_k2_hmp2_df, diagnosis == "nonIBD")$RA)
dsrAB_hmp2_UC_nonIBD <- wilcox.test(filter(dsrAB_k2_hmp2_df, diagnosis == "UC")$RA,filter(dsrAB_k2_hmp2_df, diagnosis == "nonIBD")$RA)

#### prism ####
primary_prism_CD_HC <- wilcox.test(filter(primary_k2_prism_df, diagnosis == "CD")$RA,filter(primary_k2_prism_df, diagnosis == "HC")$RA)
primary_prism_UC_HC <- wilcox.test(filter(primary_k2_prism_df, diagnosis == "UC")$RA,filter(primary_k2_prism_df, diagnosis == "HC")$RA)
secondary_prism_CD_HC <- wilcox.test(filter(secondary_k2_prism_df, diagnosis == "CD")$RA,filter(secondary_k2_prism_df, diagnosis == "HC")$RA)
secondary_prism_UC_HC <- wilcox.test(filter(secondary_k2_prism_df, diagnosis == "UC")$RA,filter(secondary_k2_prism_df, diagnosis == "HC")$RA)
dsrAB_prism_CD_HC <- wilcox.test(filter(dsrAB_k2_prism_df, diagnosis == "CD")$RA,filter(dsrAB_k2_prism_df, diagnosis == "HC")$RA)
dsrAB_prism_UC_HC <- wilcox.test(filter(dsrAB_k2_prism_df, diagnosis == "UC")$RA,filter(dsrAB_k2_prism_df, diagnosis == "HC")$RA)
#### ####

#### cib ####
primary_cib_Crohn_control <- wilcox.test(filter(primary_k2_cib_df, diagnosis == "Crohn")$RA,filter(primary_k2_cib_df, diagnosis == "Control")$RA)
secondary_cib_Crohn_control <- wilcox.test(filter(secondary_k2_cib_df, diagnosis == "Crohn")$RA,filter(secondary_k2_cib_df, diagnosis == "Control")$RA)
dsrAB_cib_Crohn_control <- wilcox.test(filter(dsrAB_k2_cib_df, diagnosis == "Crohn")$RA,filter(dsrAB_k2_cib_df, diagnosis == "Control")$RA)
#### ####


# all_pvals <- data.frame(primary_cMD_CRC_control$p.value,
# primary_cMD_adenoma_control$p.value,
# primary_cMD_adenoma_CRC$p.value,
# primary_cMD_IBD_control$p.value,
# secondary_cMD_CRC_control$p.value,
# secondary_cMD_CRC_adenoma$p.value,
# secondary_cMD_IBD_control$p.value,
# dsrAB_cMD_CRC_control$p.value,
# dsrAB_cMD_adenoma_control$p.value,
# primary_hmp2_CD_nonIBD$p.value,
# primary_hmp2_UC_nonIBD$p.value,
# secondary_hmp2_CD_nonIBD$p.value,
# secondary_hmp2_UC_nonIBD$p.value,
# dsrAB_hmp2_CD_nonIBD$p.value,
# dsrAB_hmp2_UC_nonIBD$p.value,
# primary_prism_CD_HC$p.value,
# primary_prism_UC_HC$p.value,
# secondary_prism_CD_HC$p.value,
# secondary_prism_UC_HC$p.value,
# dsrAB_prism_CD_HC$p.value,
# dsrAB_prism_UC_HC$p.value,
# primary_cib_Crohn_control$p.value,
# secondary_cib_Crohn_control$p.value,
# dsrAB_cib_Crohn_control$p.value)

# hist(-log10(all_pvals))
# axis(side=1, at=seq(0,1, 15),labels = seq(0,3,15))


##### ==================== plotting risk cMD sample data ================ #####
print("- plotting risk cMD sample data")

clrs1 <- c("#999999", "#E69F00", "#56B4E9", "#0072B2")

primary_k2_boxplot <- ggplot(primary_k2_df) + 
  geom_boxplot(aes(x = func, y = RA, fill = population)) +
  theme_bw() +
  theme(legend.position = "none",
        axis.text.x = element_blank(),
        axis.title.x = element_text(size = 8),
        # axis.title.y = element_blank(),
        panel.grid = element_blank()) +
  scale_y_continuous(limits = c(0, 100),
                     minor_breaks = NULL) +
  scale_fill_manual(values = clrs1) +
  ylab("relative abundance, (%)") + xlab("primary")
primary_k2_boxplot

secondary_k2_boxplot <- ggplot(secondary_k2_df) + 
  geom_boxplot(aes(x = func, y = RA, fill = population)) +
  theme_bw() +
  theme(legend.position = "none",
        axis.text.x = element_blank(),
        axis.title.x = element_text(size = 8),
        axis.title.y = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        panel.grid = element_blank()) +
  scale_y_continuous(limits = c(0, 100),
                     # breaks = c(0.01, 0.1, 1, 10, 100),
                     # trans = "log2",
                     minor_breaks = NULL) +
  scale_fill_manual(values = clrs1) +
  xlab("secondary")
secondary_k2_boxplot

erroneous_k2_boxplot <- ggplot(erroneous_k2_df) + 
  geom_boxplot(aes(x = func, y = RA, fill = population)) +
  theme_bw() +
  theme(legend.position = "none",
        axis.text.x = element_blank(),
        axis.title.x = element_text(size = 8),
        # axis.title.y = element_blank(),
        # axis.text.y = element_blank(),
        # axis.ticks.y = element_blank(),
        panel.grid = element_blank()) +
  scale_y_continuous(limits = c(0.001, 100),
                     breaks = c(0.01, 0.1, 1, 10, 100),
                     trans = "log2",
                     minor_breaks = NULL) +
  scale_fill_manual(values = clrs1) +
  ylab("relative abundance, log2(%)") + xlab("erroneous")
erroneous_k2_boxplot

dsrAB_k2_boxplot <- ggplot(dsrAB_k2_df) + 
  geom_boxplot(aes(x = func, y = RA, fill = population)) +
  theme_bw() +
  theme(legend.position = "none",
        axis.title.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.text.y = element_blank(),
        axis.text.x = element_blank(),
        axis.title.x = element_text(size = 8),
        panel.grid = element_blank()) +
  scale_y_continuous(limits = c(0.001, 100),
                     breaks = c(0.01, 0.1, 1, 10, 100),
                     trans = "log2",
                     minor_breaks = NULL) +
  scale_fill_manual(values = clrs1) +
  ylab("relative abundance, log2(%)") + xlab("SRB")
dsrAB_k2_boxplot

## just for getting the legend for plotting in illustrator
p1_legend_plt <- ggplot(dsrAB_k2_df) + 
  geom_boxplot(aes(x = func, y = RA, fill = population)) +
  theme_bw() +
  scale_fill_manual(values = clrs1)
grid.newpage()
grid.draw(get_legend(p1_legend_plt))
ggsave("../../figures/figure3/p1_legend_raw.svg", height = 2, width = 2)
ggsave("../../figures/figure3/p1_legend_raw.png", height = 2, width = 2)

print("- arranging cMD plots")

# (cys_k2_boxplot | dsrAB_k2_boxplot) +
#   plot_annotation(title = "curatedMetagneomicData",
#                   theme = theme(plot.title = element_text(size = 22)))
# ggsave("../../figures/figure4/k2_cMD_RA_risk_populations.svg", height = 7, width = 7)
# ggsave("../../figures/figure4/k2_cMD_RA_risk_populations.png", height = 7, width = 7)
##### ==================== plotting hmp2 data =========================== #####
print("- plotting hmp2 data")

clrs2 <- c("#999999", "#F0E442", "#D55E00")

primary_k2_hmp2_boxplot <- ggplot(primary_k2_hmp2_df) + 
  geom_boxplot(aes(x = func, y = RA, fill = diagnosis)) +
  theme_bw() +
  theme(legend.position = "none",
        axis.text.x = element_blank(),
        panel.grid.minor.y = element_blank(),
        axis.title.x = element_text(size = 8),
        # axis.title.y = element_blank(),
        panel.grid = element_blank()) +
  scale_y_continuous(limits = c(0, 100)) +
  scale_fill_manual(values = clrs2) +
  ylab("relative abundance (%)") + xlab("primary")
primary_k2_hmp2_boxplot

secondary_k2_hmp2_boxplot <- ggplot(secondary_k2_hmp2_df) + 
  geom_boxplot(aes(x = func, y = RA, fill = diagnosis)) +
  theme_bw() +
  theme(legend.position = "none",
        axis.text.x = element_blank(),
        panel.grid.minor.y = element_blank(),
        axis.title.x = element_text(size = 8),
        axis.title.y = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        panel.grid = element_blank()) +
  scale_y_continuous(limits = c(0, 100)) +
  scale_fill_manual(values = clrs2) +
  ylab("relative bacterial abundance (%)") + xlab("secondary")
secondary_k2_hmp2_boxplot

erroneous_k2_hmp2_boxplot <- ggplot(erroneous_k2_hmp2_df) + 
  geom_boxplot(aes(x = func, y = RA, fill = diagnosis)) +
  theme_bw() +
  theme(legend.position = "none",
        axis.text.x = element_blank(),
        panel.grid.minor.y = element_blank(),
        axis.title.x = element_text(size = 8),
        # axis.title.y = element_blank(),
        # axis.text.y = element_blank(),
        # axis.ticks.y = element_blank(),
        panel.grid = element_blank()) +
  scale_y_continuous(limits = c(0.001, 100),
                     breaks = c(0.01, 0.1, 1, 10, 100),
                     trans = "log2",
                     minor_breaks = NULL) +
  scale_fill_manual(values = clrs2) +
  ylab("relative abundance, log2(%)") + xlab("erroneous")
erroneous_k2_hmp2_boxplot

dsrAB_k2_hmp2_boxplot <- ggplot(dsrAB_k2_hmp2_df) + 
  geom_boxplot(aes(x = func, y = RA, fill = diagnosis)) +
  theme_bw() +
  theme(legend.position = "none",
        axis.title.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.text.y = element_blank(),
        axis.text.x = element_blank(),
        axis.title.x = element_text(size = 8),
        panel.grid = element_blank()) +
  scale_y_continuous(limits = c(0.001, 100),
                     breaks = c(0.01, 0.1, 1, 10, 100),
                     trans = "log2",
                     minor_breaks = NULL) +
  scale_fill_manual(values = clrs2) +
  xlab("SRB")
dsrAB_k2_hmp2_boxplot

## just for getting the legend for plotting in illustrator
p2_legend_plt <- ggplot(dsrAB_k2_hmp2_df) + 
  geom_boxplot(aes(x = func, y = RA, fill = diagnosis)) +
  theme_bw() +
  scale_fill_manual(values = clrs2)
grid.newpage()
grid.draw(get_legend(p2_legend_plt))
ggsave("../../figures/figure3/p2_legend_raw.svg", height = 2, width = 2)
ggsave("../../figures/figure3/p2_legend_raw.png", height = 2, width = 2)

##### ==================== plotting prism data ========================== #####
print("- plotting prism data")

clrs3 <- c("#999999", "#F0E442", "#D55E00")

primary_k2_prism_boxplot <- ggplot(primary_k2_prism_df) + 
  geom_boxplot(aes(x = func, y = RA, fill = diagnosis)) +
  theme_bw() +
  theme(legend.position = "none",
        axis.text.x = element_blank(),
        panel.grid.minor.y = element_blank(),
        axis.title.x = element_text(size = 8),
        # axis.title.y = element_blank(),
        panel.grid = element_blank()) +
  scale_y_continuous(limits = c(0, 100)) +
  scale_fill_manual(values = clrs3) +
  ylab("relative abundance, (%)") + xlab("primary")
primary_k2_prism_boxplot

secondary_k2_prism_boxplot <- ggplot(secondary_k2_prism_df) + 
  geom_boxplot(aes(x = func, y = RA, fill = diagnosis)) +
  theme_bw() +
  theme(legend.position = "none",
        axis.text.x = element_blank(),
        panel.grid.minor.y = element_blank(),
        axis.title.x = element_text(size = 8),
        axis.title.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.text.y = element_blank(),
        panel.grid = element_blank()) +
  scale_y_continuous(limits = c(0, 100)) +
  scale_fill_manual(values = clrs3) +
  ylab("relative bacterial abundance, (%)") + xlab("secondary")
secondary_k2_prism_boxplot

erroneous_k2_prism_boxplot <- ggplot(erroneous_k2_prism_df) + 
  geom_boxplot(aes(x = func, y = RA, fill = diagnosis)) +
  theme_bw() +
  theme(legend.position = "none",
        axis.text.x = element_blank(),
        panel.grid.minor.y = element_blank(),
        axis.title.x = element_text(size = 8),
        # axis.title.y = element_blank(),
        # axis.ticks.y = element_blank(),
        # axis.text.y = element_blank(),
        panel.grid = element_blank()) +
  scale_y_continuous(limits = c(0.001, 100),
                     breaks = c(0.01, 0.1, 1, 10, 100),
                     trans = "log2",
                     minor_breaks = NULL) +
  scale_fill_manual(values = clrs3) +
  ylab("relative abundance, log2(%)") + xlab("erroneous")
erroneous_k2_prism_boxplot

dsrAB_k2_prism_boxplot <- ggplot(dsrAB_k2_prism_df) + 
  geom_boxplot(aes(x = func, y = RA, fill = diagnosis)) +
  theme_bw() +
  theme(legend.position = "none",
        axis.title.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.text.y = element_blank(),
        axis.text.x = element_blank(),
        axis.title.x = element_text(size = 8),
        panel.grid = element_blank()) +
  scale_y_continuous(limits = c(0.001, 100),
                     breaks = c(0.01, 0.1, 1, 10, 100),
                     trans = "log2",
                     minor_breaks = NULL) +
  scale_fill_manual(values = clrs3) +
  xlab("SRB")
# dsrAB_k2_prism_boxplot

## just for getting the legend for plotting in illustrator
p3_legend_plt <- ggplot(dsrAB_k2_prism_df) + 
  geom_boxplot(aes(x = func, y = RA, fill = diagnosis)) +
  theme_bw() +
  scale_fill_manual(values = clrs3)
grid.newpage()
grid.draw(get_legend(p3_legend_plt))
ggsave("../../figures/figure3/p3_legend_raw.svg", height = 2, width = 2)
ggsave("../../figures/figure3/p3_legend_raw.png", height = 2, width = 2)

##### ==================== plotting cib data ============================ #####
print("- plotting cib data")

clrs4 <- c("#999999", "#D55E00")

primary_k2_cib_boxplot <- ggplot(primary_k2_cib_df) + 
  geom_boxplot(aes(x = func, y = RA, fill = diagnosis)) +
  theme_bw() +
  theme(legend.position = "none",
        axis.text.x = element_blank(),
        panel.grid.minor.y = element_blank(),
        axis.title.x = element_text(size = 8),
        # axis.title.y = element_blank(),
        panel.grid = element_blank()) +
  scale_y_continuous(limits = c(0, 100)) +
  scale_fill_manual(values = clrs4) +
  ylab("relative abundance, (%)") + xlab("primary")
primary_k2_cib_boxplot

secondary_k2_cib_boxplot <- ggplot(secondary_k2_cib_df) + 
  geom_boxplot(aes(x = func, y = RA, fill = diagnosis)) +
  theme_bw() +
  theme(legend.position = "none",
        axis.text.x = element_blank(),
        panel.grid.minor.y = element_blank(),
        axis.title.x = element_text(size = 8),
        axis.title.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.text.y = element_blank(),
        panel.grid = element_blank()) +
  scale_y_continuous(limits = c(0, 100)) +
  scale_fill_manual(values = clrs4) +
  ylab("relative bacterial abundance, log2[%]") + xlab("secondary")
secondary_k2_cib_boxplot

erroneous_k2_cib_boxplot <- ggplot(erroneous_k2_cib_df) + 
  geom_boxplot(aes(x = func, y = RA, fill = diagnosis)) +
  theme_bw() +
  theme(legend.position = "none",
        axis.text.x = element_blank(),
        panel.grid.minor.y = element_blank(),
        axis.title.x = element_text(size = 8),
        # axis.title.y = element_blank(),
        # axis.ticks.y = element_blank(),
        # axis.text.y = element_blank(),
        panel.grid = element_blank()) +
  scale_y_continuous(limits = c(0.001, 100),
                     breaks = c(0.01, 0.1, 1, 10, 100),
                     trans = "log2",
                     minor_breaks = NULL) +
  scale_fill_manual(values = clrs4) +
  ylab("relative abundance, log2(%)") + xlab("erroneous")
erroneous_k2_cib_boxplot

dsrAB_k2_cib_boxplot <- ggplot(dsrAB_k2_cib_df) + 
  geom_boxplot(aes(x = func, y = RA, fill = diagnosis)) +
  theme_bw() +
  theme(legend.position = "none",
        axis.title.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.text.y = element_blank(),
        axis.text.x = element_blank(),
        axis.title.x = element_text(size = 8),
        panel.grid = element_blank()) +
  scale_y_continuous(limits = c(0.001, 100),
                     breaks = c(0.01, 0.1, 1, 10, 100),
                     trans = "log2",
                     minor_breaks = NULL) +
  scale_fill_manual(values = clrs4) +
  xlab("SRB")
dsrAB_k2_cib_boxplot

## just for getting the legend for plotting in illustrator
p4_legend_plt <- ggplot(dsrAB_k2_cib_df) + 
  geom_boxplot(aes(x = func, y = RA, fill = diagnosis)) +
  theme_bw() +
  scale_fill_manual(values = clrs4)
grid.newpage()
grid.draw(get_legend(p4_legend_plt))
ggsave("../../figures/figure3/p4_legend_raw.svg", height = 2, width = 2)
ggsave("../../figures/figure3/p4_legend_raw.png", height = 2, width = 2)

