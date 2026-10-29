# Figure 4A: AREG expression and glycolysis in baseline OPMDs (GSE26549)
# Requires GEOquery (Bioconductor) and ggplot2 (CRAN).
# Run from the desired project directory; GEO download requires internet access.

suppressPackageStartupMessages({
  library(GEOquery)
  library(ggplot2)
})

dir.create("data_raw", showWarnings = FALSE, recursive = TRUE)
dir.create("results", showWarnings = FALSE, recursive = TRUE)

gse <- GEOquery::getGEO("GSE26549", GSEMatrix = TRUE, getGPL = TRUE,
                        destdir = "data_raw")
eset <- gse[[1]]
meta <- Biobase::pData(eset)
clinical <- data.frame(
  sample = as.character(meta$geo_accession),
  biopsy_time = as.character(meta$`time of biopsy:ch1`),
  stringsAsFactors = FALSE
)
clinical <- clinical[clinical$biopsy_time == "biopsy at baseline", , drop = FALSE]
stopifnot(nrow(clinical) == 62L, !anyDuplicated(clinical$sample))

mat <- Biobase::exprs(eset)
areg_probes <- c("8095736", "8095744")
gly_probes <- c(SLC2A1 = "7915472", HK2 = "8042942", PFKP = "7925876",
                ALDOA = "7994737", ENO1 = "7912198", PKM = "7990151",
                LDHA = "7938777", PDK1 = "8046408", SLC16A1 = "7918622",
                SLC16A3 = "8010770")
stopifnot(all(clinical$sample %in% colnames(mat)),
          all(c(areg_probes, unname(gly_probes)) %in% rownames(mat)))

samples <- clinical$sample
clinical$AREG <- colMeans(mat[areg_probes, samples, drop = FALSE])
gly_expr <- t(mat[unname(gly_probes), samples, drop = FALSE])
colnames(gly_expr) <- names(gly_probes)
# Each gene is Z-standardized across the 62 baseline samples.
clinical$Glycolysis_score <- rowMeans(scale(gly_expr))
clinical$AREG_group <- factor(
  ifelse(clinical$AREG > median(clinical$AREG), "AREG-high", "AREG-low"),
  levels = c("AREG-low", "AREG-high")
)
stopifnot(!anyNA(clinical[, c("AREG", "Glycolysis_score")]))

spearman <- cor.test(clinical$AREG, clinical$Glycolysis_score,
                     method = "spearman", exact = FALSE)
wilcox <- wilcox.test(Glycolysis_score ~ AREG_group,
                      data = clinical, exact = FALSE)
print(spearman)
print(wilcox)
print(table(clinical$AREG_group))
write.csv(
  data.frame(test = c("Spearman AREG vs glycolysis", "Wilcoxon AREG groups"),
             estimate = c(unname(spearman$estimate), NA_real_),
             p_value = c(spearman$p.value, wilcox$p.value)),
  "results/Figure_4A_statistics.csv", row.names = FALSE
)

p <- ggplot(clinical, aes(AREG_group, Glycolysis_score,
                          fill = AREG_group, color = AREG_group)) +
  geom_violin(trim = FALSE, alpha = 0.35, width = 0.85, linewidth = 0.65) +
  geom_boxplot(aes(group = AREG_group), width = 0.20, outlier.shape = NA,
               fill = "white", color = "black", alpha = 0.65,
               linewidth = 0.65) +
  geom_jitter(width = 0.08, size = 2.1, alpha = 0.95) +
  scale_fill_manual(values = c("AREG-low" = "#A9D9ED", "AREG-high" = "#F5C5AC")) +
  scale_color_manual(values = c("AREG-low" = "#0876B9", "AREG-high" = "#E65B12")) +
  labs(x = NULL, y = "Glycolysis score",
       title = "Baseline oral potentially malignant disorders",
       subtitle = sprintf("GSE26549 | Wilcoxon P = %.4f", wilcox$p.value)) +
  theme_classic(base_size = 13) +
  theme(legend.position = "none", plot.title = element_text(face = "bold"))

print(p)
# Export manually in RStudio: Plots > Export > Save as PDF.
writeLines(capture.output(sessionInfo()), "results/Figure_4A_sessionInfo.txt")
