# ======================= H2S_CH4_presence_absence.R ======================== #
#' description: []
#' 
#' NOTE: it is assumed that all necessary packages and data are loaded into the 
#' R session where this script is being run.
# =========================================================================== #

dim(RPKM_H2S)
dim(RPKM_CH4)

df <- data.frame(matrix(nrow = ncol(RPKM_H2S), ncol = 3))
rownames(df) <- colnames(RPKM_H2S)
colnames(df) <- c("cysteine_degrading", "sulfate_reducing", "methane_producing")

for (i in 1:nrow(df)) {
  # print("-- ON ", i)
  if (any(RPKM_H2S[c(1:10, 12:18), i] > 1)) {
    df[i, 1] <- 1
  } else {
    df[i, 1] <- 0
  }
  
  if (RPKM_H2S[19, i] > 1 && RPKM_H2S[20, i] > 1) {
    df[i, 2] <- 1
  } else {
    df[i, 2] <- 0
  }
  
  if (sum(RPKM_CH4[, i] > 0) / length(RPKM_CH4[, 1] > 0) > 0.9) {
    df[i, 3] <- 1
  } else {
    df[i, 3] <- 0
  }
}
