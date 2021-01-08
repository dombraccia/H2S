# ====================== healthy_relative_abundances.R ======================= #
#' description: script for comparing the relative abundances of cysteine 
#' degrading and sulfite reducing bacteria in stool samples from cMD.
#' 
#' NOTE: it is assumed that all necessary packages and data are loaded into the 
#' R session where this script is being run.
# =========================================================================== #

print("- importing taxa hit information and feature information")
dsr_taxa_hits_raw <- read.csv(
  file = "../results/from-GutFunFind/Dissimilatory_Sulf_Reduction.taxa_hits.txt", 
  sep = "\t", header = FALSE)
cys_taxa_hits_raw <- read.csv(
  file = "../results/from-GutFunFind/Cysteine_Degradation.taxa_hits.txt", 
  sep = "\t", header = FALSE)
uhgg_genid2taxa <- read.csv("../data/from-xiaofang/spec.txt", sep = '\t', header = FALSE)
uhgg_genid2taxa$V2 <- gsub(";", "|", uhgg_genid2taxa$V2)

dsr_taxa_hits <- as.vector(t(dsr_taxa_hits_raw), mode = "character") 
dsr_taxa_hits <- gsub(";", "|", dsr_taxa_hits)
cys_taxa_hits <- as.vector(t(cys_taxa_hits_raw), mode = "character")
cys_taxa_hits <- gsub(";", "|", cys_taxa_hits)

print("- selecting SRB that contain dsrAB")
dsr_overlap %>% 
  filter(dsrAB == 1) %>% 
  select(genome) -> dsrAB_genids
dsrAB_taxa <- uhgg_genid2taxa$V2[uhgg_genid2taxa$V1 %in% as.vector(t(dsrAB_genids), 
                                                                   mode = "character")]

## subsetting kraken2 RA controls data for SRB
dsr_k2_RA_controls <- k2_RA_controls[dsrAB_taxa, ]
dsr_k2_cs_controls <- colSums(dsr_k2_RA_controls)

print("- subsetting kraken2 RA data for cys degrading bacteria")
cys_genes %>%
  filter(genes != "absent") %>%
  select(genome) -> cys_genids
cys_taxa <- uhgg_genid2taxa$V2[uhgg_genid2taxa$V1 %in% as.vector(t(cys_genids), mode = "character")]
cys_k2_RA_controls <- k2_RA_controls[cys_taxa, ]
cys_k2_cs_controls <- colSums(cys_k2_RA_controls)

## prepping data for plotting
df <- data.frame(cysteine_degrading = cys_k2_cs_controls,
                 sulfite_reducing = dsr_k2_cs_controls)

print("- plotting abundances of sulfite red bac and cys deg bac")
cys_met_dsr_healthy_RA <- ggplot(melt(df), aes(x = factor(variable), y = value, fill = factor(variable))) +
  geom_violin() +
  # geom_boxplot(outlier.alpha = 0.12) + #geom_jitter(aes(alpha = 0.001)) +
  stat_summary(fun = median, geom = "point", shape = 1, size = 4) +
  stat_summary(fun = mean, geom = "point", shape = 4, size = 4) +
  scale_shape_manual("", values=c("mean" = "X")) +
  scale_fill_manual(values = c("#dd9e2f", "#1c6bc9")) +
  scale_y_continuous(limits = c(0.001, 100),
                     breaks = c(0.01, 0.1, 1, 10, 100),  
                     trans = "log2", 
                     labels = label_comma(accuracy = 0.01)) +
  ggtitle("Relative microbial abundances across H2S production pathways") +
  ylab("relative microbial abundances, log2(%)") + labs(fill = "Function") + 
  theme_bw() +
  theme(axis.title.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.text.x = element_text(size = 13),
        axis.text.y = element_text(size = 13),
        axis.title.y = element_text(size = 13),
        # plot.title = element_text(size = 14),
        plot.title = element_blank(),
        legend.position = "none") 
cys_met_dsr_healthy_RA
ggsave("../figures/figure2/k2_RAs_healthycontrols.svg", 
       plot = cys_met_dsr_healthy_RA, width = 5, height = 7)
ggsave("../figures/figure2/k2_RAs_healthycontrols.png", 
       plot = cys_met_dsr_healthy_RA, width = 5, height = 7)
