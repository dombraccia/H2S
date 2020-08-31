## running summary plots from curatedMetagenomicData on cysteine-degrading bacteria

##### ============================ LOAD LIBS ============================ #####

library(curatedMetagenomicData)
library(tictoc)

##### ========================== DOWNLOAD DATA ========================== #####

tic()
tmp <- curatedMetagenomicData("*metaphlan_bugs_list.stool*", dryrun = FALSE)
toc()

##### ============================ <> ============================ #####
