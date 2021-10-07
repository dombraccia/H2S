# ==================== healthy_relative_abundances_hmm.R ==================== #
#' description: script for comparing the relative abundances of putative 
#' cysteine degrading and sulfite reducing bacteria in stool samples from cMD.
#' 
#' NOTE: it is assumed that all necessary packages and data are loaded into the 
#' R session where this script is being run.
# =========================================================================== #

print("- importing taxa hit information and feature information")
dsr_taxa_hits_raw <- read.csv(
  file = "results/from-GutFunFind/Dissimilatory_Sulf_Reduction.taxa_hits.txt", 
  sep = "\t", header = FALSE)
cys_taxa_hits_raw <- read.csv(
  file = "results/from-GutFunFind/Cysteine_Degradation_hmm.taxa_hits_hmm.tsv", 
  sep = "\t", header = TRUE)
cys_feature <- read.csv(
  file = "results/from-GutFunFind/Cysteine_Degradation_hmm.feature_hmm.tsv", 
  sep = "\t", header = TRUE)

uhgg_genid2taxa <- read.csv("data/from-xiaofang/spec.txt", sep = '\t', header = FALSE)
uhgg_genid2taxa$V2 <- gsub(";", "|", uhgg_genid2taxa$V2)
colnames(uhgg_genid2taxa) <- c("genome_ID", "taxa_ID")

dsr_taxa_hits <- as.vector(t(dsr_taxa_hits_raw), mode = "character") 
dsr_taxa_hits <- gsub(";", "|", dsr_taxa_hits)
cys_taxa_hits <- as.vector(t(cys_taxa_hits_raw), mode = "character")
cys_taxa_hits <- gsub(";", "|", cys_taxa_hits)

print("- selecting SRB that contain dsrAB")
dsr_overlap %>% 
  filter(dsrAB == 1) %>% 
  select(genome) -> dsrAB_genids

dsrAB_taxa <- uhgg_genid2taxa$taxa_ID[uhgg_genid2taxa$genome_ID %in% as.vector(t(dsrAB_genids), 
                                                                   mode = "character")]

print("- subsetting kraken2 RA controls data for SRB")
dsr_k2_RA_controls <- k2_RA_controls[dsrAB_taxa, ]
dsr_k2_cs_controls <- colSums(dsr_k2_RA_controls)

print("- subsetting kraken2 RA data for cys degrading bacteria")

taxa_feature <- tibble(taxa_ID = uhgg_genid2taxa$taxa_ID, 
                       output_type = cys_feature$output_type)

taxa_feature %>% 
  filter(output_type == "primary") %>% 
  select(taxa_ID) -> primary_taxa
primary_k2_RA_controls <- k2_RA_controls[primary_taxa$taxa_ID, ]
  primary_k2_cs_controls <- colSums(primary_k2_RA_controls)

taxa_feature %>% 
  filter(output_type == "secondary") %>% 
  select(taxa_ID) -> secondary_taxa
secondary_k2_RA_controls <- k2_RA_controls[secondary_taxa$taxa_ID, ]
secondary_k2_cs_controls <- colSums(secondary_k2_RA_controls)

taxa_feature %>% 
  filter(output_type == "erroneous") %>% 
  select(taxa_ID) -> erroneous_taxa
erroneous_k2_RA_controls <- k2_RA_controls[erroneous_taxa$taxa_ID, ]
erroneous_k2_cs_controls <- colSums(erroneous_k2_RA_controls)

## prepping data for plotting
df <- data.frame(primary = primary_k2_cs_controls,
                 secondary = secondary_k2_cs_controls,
                 erroneous = erroneous_k2_cs_controls,
                 SRB = dsr_k2_cs_controls)
dim(df)
df <- filter(df, primary <= 100 & secondary <= 100 & erroneous <= 100)

## calculating p-values comparing relative abundance between PCDB and SRB
cMD_healthy_k2_primary_SRB <- 
  wilcox.test(df$primary, 
              df$SRB)
cMD_healthy_k2_primary_SRB

print("- plotting abundances of sulfite red bac and cys deg bac")
cys_met_dsr_healthy_RA <- ggplot(melt(df), aes(x = factor(variable), y = value, fill = factor(variable))) +
  # geom_violin(scale = "width") +
  geom_boxplot(outlier.alpha = 0.12) + #geom_jitter(aes(alpha = 0.0001), width = 0.25) +
  # stat_summary(fun = median, geom = "point", shape = 1, size = 4) +
  # stat_summary(fun = mean, geom = "point", shape = 4, size = 4) +
  scale_shape_manual("", values=c("mean" = "X")) +
  # scale_fill_manual(values = c("#dd9e2f", "#1c6bc9")) +
  # scale_y_continuous(limits = c(0.001, 100),
                     # breaks = c(0.01, 0.1, 1, 10, 100),  
                     # trans = "log2", 
                     # labels = label_comma(accuracy = 0.01)) +
  ggtitle("Relative microbial abundances across H2S production pathways") +
  ylab("relative microbial abundances, log2(%)") + labs(fill = "Function") + 
  theme_bw() +
  theme(axis.title.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.text.x = element_text(size = 14),
        axis.text.y = element_text(size = 14),
        axis.title.y = element_text(size = 15),
        panel.grid.minor.x = element_blank(),
        panel.grid.minor.y = element_blank(),
        plot.title = element_text(size = 14),
        # plot.title = element_blank(),
        legend.position = "none") +
  ggtitle("Relative Abundances of various H2S producing bacteria")
cys_met_dsr_healthy_RA
ggsave("figures/figure3/k2_RAs_healthycontrols.svg", 
       plot = cys_met_dsr_healthy_RA, width = 6, height = 7)
ggsave("figures/figure3/k2_RAs_healthycontrols.png", 
       plot = cys_met_dsr_healthy_RA, width = 6, height = 7)
