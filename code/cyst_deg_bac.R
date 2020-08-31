## running summary plots from curatedMetagenomicData on the prevalence 
## cysteine-degrading bacteria

##### ============================ LOAD LIBS ============================ #####

library(curatedMetagenomicData)
library(tictoc)
library(dplyr)

##### ========================== DOWNLOAD DATA ========================== #####

if(file.exists("data/all_data.RDS")) {
  all_data <- readRDS("data/all_data.RDS")
} else {
  tic()
  all_data <- curatedMetagenomicData("*metaphlan_bugs_list.stool*", dryrun = FALSE)
  saveRDS(all_data, file = "data/all_data.RDS")
  toc()
}

##### ======================== EXPLORING DATA =========================== #####

## merging all datasets
all_data_merged <- mergeData(all_data)

all_pData <- pData(all_data_merged)
disease_type <- all_pData %>% 
  select("study_condition", "disease", "disease_subtype", "disease_stage")
sum(is.na(all_pData$disease))
sum(all_pData$disease == "healthy")
study_cond <- all_pData$study_condition
study_cond_NA <- all_data_merged[,is.na(all_pData$study_condition) ]

##### ============================ SUBSET DATA ============================ #####

## removing NAs from pData
all_pData <- all_pData[!is.na(all_pData$study_condition), ]

## subsetting pData to study_condition == "controls"
all_pData_controls <- all_pData %>% filter(study_condition == "control")

## getting rid of NAs in relative abundance (RA) data
all_RA <- exprs(all_data_merged)
all_RA <- all_RA[, !is.na(all_pData$study_condition)]

## subsetting RA down to just present in control samples
control_RA <- all_RA[, colnames(all_RA) %in% rownames(all_pData_controls)]


##### ============================ <> ============================ #####
