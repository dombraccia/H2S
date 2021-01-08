# ============================= process_data.R ============================== #
#' description: initial processing and subetting step for imported data
#' 
#' NOTE: it is assumed that all necessary packages and data are loaded into the 
#' R session where this script is being run.
# =========================================================================== #

# ========================= cMD DATA (Metaphlan) ============================ #

print("- loading in saved cMD data (metaphlan RA)")
if(file.exists("../data/all_data.RDS")) {
  m2_data <- readRDS("../data/all_data.RDS")
} else {
  tic()
  m2_data <- curatedMetagenomicData("*metaphlan_bugs_list.stool*", dryrun = FALSE)
  saveRDS(all_data, file = "../data/all_data.RDS")
  toc()
}

# ====================== kraken2 data from Xiaofang ========================= #
print("- loading kraken2 RA from Xiaofang")
k2_counts <- readRDS("../data/from-xiaofang/kraken2_output.rds")
k2_pData <- readRDS("../data/from-xiaofang/metadata.rds")
k2_hmp2_counts <- read.csv("../data/from-xiaofang/HMP2-IBD.mpa.txt", 
                           sep = "\t", header = TRUE, row.names = 1)
k2_hmp2_pData <- read.csv("../data/from-xiaofang/hmp2_metadata.MGX.csv", 
                          header = TRUE) # NOTE: `diagnosis` column contains study condition information.
k2_hmppilot_counts <- read.csv("../data/from-xiaofang/HMP-IBD-pilot.mpa.txt", 
                               sep = "\t", header = TRUE)
k2_prism_counts <- read.csv("../data/from-xiaofang/merged.mpa.txt", 
                            sep = "\t", header = TRUE, row.names = 1) %>%
                          as.matrix()
k2_prism_pData <- read.csv("../data/prism_combined_metadata_20201209.txt", 
                           sep = "\t", header = TRUE)
k2_cib_counts <- read.csv("../data/from-xiaofang/cib.mpa.txt", 
                          sep = "\t", row.names = 1) %>% 
                          as.matrix()
k2_cib_pData <- read.csv("../data/from-xiaofang/cib_metadata.tsv",
                         sep = "\t")

# ============================ METAPHLAN2 DATA ============================== #
print("- processing metaphlan RA data")
## merging all datasets
m2_data_merged <- mergeData(m2_data)
m2_pData <- pData(m2_data_merged)
study_cond <- m2_pData$study_condition
study_cond_NA <- m2_data_merged[,is.na(m2_pData$study_condition) ]

## removing NAs from pData
m2_pData <- m2_pData[!is.na(m2_pData$study_condition), ]

## subsetting pData to study_condition == "controls"
m2_pData_controls <- m2_pData %>% filter(study_condition == "control")

## getting rid of NAs in relative abundance (RA) data
m2_RA <- exprs(m2_data_merged)
m2_RA <- m2_RA[, !is.na(m2_pData$study_condition)]

## subsetting RA down to just present in control samples
m2_RA_controls <- m2_RA[, colnames(m2_RA) %in% rownames(m2_pData_controls)]

# =============================== KRAKEN DATA =============================== #
print("- processing kraken RA data (hmp2)")
## computing RA values from count info
k2_counts_spp <- k2_counts[grep("s__", rownames(k2_counts)), ]
k2_RA <- prop.table(k2_counts_spp, 2) * 100 ## counts -> RA (%)
k2_hmp2_counts_spp <-  k2_hmp2_counts[grep("s__", rownames(k2_hmp2_counts)), ]
k2_hmp2_RA <- prop.table(as.matrix(k2_hmp2_counts_spp), 2) * 100 ## counts -> RA (%)

## getting rid of RA data whose sample has study_condition == NA
study_cond_NA <- k2_pData %>% filter(is.na(study_condition)) %>% rownames()
k2_RA <- select(as.data.frame(k2_RA), -study_cond_NA)
k2_pData <- filter(k2_pData, !is.na(study_condition))
k2_pData %>% 
  filter(study_condition == "control") -> k2_pData_controls
k2_RA_controls <- k2_RA[, colnames(k2_RA) %in% rownames(k2_pData_controls)]

## removing duplicate metadata sample entries
k2_hmp2_pData <- k2_hmp2_pData[-grep("_", k2_hmp2_pData$External.ID), ]

## modifying k2_hmp2_pData$External.ID column for compatability with RA
# k2_hmp2_pData$External.ID.mod <- sapply(strsplit(k2_hmp2_pData$External.ID, "_"), "[", 1)

# =========================================================================== #

print("- processing kraken RA data (hmp2)")
## computing RA values from count info
k2_prism_counts_spp <- k2_prism_counts[grep("s__", rownames(k2_prism_counts)), ]
k2_prism_RA <- prop.table(k2_prism_counts_spp, 2) * 100 ## counts -> RA (%)
k2_prism_counts_spp <-  k2_prism_counts_spp[grep("s__", rownames(k2_prism_counts_spp)), ]
k2_prism_RA <- prop.table(as.matrix(k2_prism_counts_spp), 2) * 100 ## counts -> RA (%)


#### ===================================================================== ####
print("- processing kraken RA data (cib)")
#### ===================================================================== ####

## computing RA values from count info
k2_cib_counts_spp <- k2_cib_counts[grep("s__", rownames(k2_cib_counts)), ]
k2_cib_RA <- prop.table(k2_cib_counts_spp, 2) * 100 ## counts -> RA (%)
k2_cib_counts_spp <-  k2_cib_counts_spp[grep("s__", rownames(k2_cib_counts_spp)), ]
k2_cib_RA <- prop.table(as.matrix(k2_cib_counts_spp), 2) * 100 ## counts -> RA (%)

## remove pData samples not present in counts
k2_cib_pData <- k2_cib_pData[k2_cib_pData$Sample_SRA %in% colnames(k2_cib_RA), ]

