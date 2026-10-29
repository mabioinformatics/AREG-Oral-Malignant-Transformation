# Figure 4B: 10-year oral cancer-free survival in baseline OPMDs (GSE26549)
# Requires GEOquery (Bioconductor), ggplot2, survival, and survminer (CRAN).
# Run from the desired project directory; GEO download requires internet access.

suppressPackageStartupMessages({
  library(GEOquery)
  library(ggplot2)
  library(survival)
  library(survminer)
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
  outcome = as.character(meta$`outcome:ch1`),
  survival_years = as.numeric(as.character(
    meta$`oral cancer-free survival time (years):ch1`)),
  stringsAsFactors = FALSE
)
clinical <- clinical[clinical$biopsy_time == "biopsy at baseline", , drop = FALSE]
clinical$event <- as.integer(clinical$outcome == "oral cancer development")
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
stopifnot(!anyNA(clinical[, c("AREG", "Glycolysis_score", "event",
                              "survival_years")]),
          all(clinical$survival_years >= 0))

# High groups are strictly above the baseline median; ties belong to low.
clinical$AREG_group <- factor(
  ifelse(clinical$AREG > median(clinical$AREG), "AREG-high", "AREG-low"),
  levels = c("AREG-low", "AREG-high")
)
clinical$Glycolysis_group <- ifelse(
  clinical$Glycolysis_score > median(clinical$Glycolysis_score),
  "Glycolysis-high", "Glycolysis-low"
)
clinical$Coupled_group <- factor(
  ifelse(clinical$AREG_group == "AREG-high" &
           clinical$Glycolysis_group == "Glycolysis-high",
         "AREG-high / Glycolysis-high", "Other OPMDs"),
  levels = c("Other OPMDs", "AREG-high / Glycolysis-high")
)
stopifnot(sum(clinical$Coupled_group == "Other OPMDs") == 45L,
          sum(clinical$Coupled_group == "AREG-high / Glycolysis-high") == 17L)

# Separate display factor keeps risk-table labels aligned with group counts.
clinical$Plot_group <- factor(
  clinical$Coupled_group,
  levels = levels(clinical$Coupled_group),
  labels = c("Other OPMDs", "AREG-high /\nGlycolysis-high")
)
stopifnot(identical(as.integer(clinical$Plot_group),
                    as.integer(clinical$Coupled_group)))

# Administratively censor follow-up at 10 years.
clinical$time_10yr <- pmin(clinical$survival_years, 10)
clinical$event_10yr <- as.integer(
  clinical$event == 1 & clinical$survival_years <= 10
)
print(table(clinical$Coupled_group, clinical$event_10yr))

fit <- survfit(Surv(time_10yr, event_10yr) ~ Plot_group, data = clinical)
logrank <- survdiff(Surv(time_10yr, event_10yr) ~ Coupled_group,
                    data = clinical)
p_logrank <- pchisq(logrank$chisq, df = 1, lower.tail = FALSE)
cox <- coxph(Surv(time_10yr, event_10yr) ~ Coupled_group, data = clinical)
print(logrank)
print(p_logrank)
print(summary(cox))
write.csv(
  data.frame(logrank_chisq = logrank$chisq, logrank_p = p_logrank,
             hazard_ratio = unname(exp(coef(cox)))),
  "results/Figure_4B_statistics.csv", row.names = FALSE
)

plot_obj <- ggsurvplot(
  fit, data = clinical,
  risk.table = TRUE, risk.table.height = 0.22,
  risk.table.y.text = TRUE, risk.table.y.text.col = FALSE,
  pval = TRUE, conf.int = FALSE,
  censor = TRUE, censor.shape = "|", censor.size = 3,
  palette = c("#0876B9", "#E65B12"),
  xlim = c(0, 10), break.time.by = 2,
  xlab = "Follow-up (years)", ylab = "Malignancy-free survival",
  title = NULL, legend.title = "",
  legend.labs = c("Other OPMDs", "AREG-high /\nGlycolysis-high"),
  ggtheme = theme_classic(base_size = 14)
)
plot_obj$plot <- plot_obj$plot +
  theme(
    axis.title.y = element_text(size = 13, margin = margin(r = -12)),
    axis.title.x = element_text(size = 13),
    axis.text = element_text(size = 11),
    legend.title = element_blank(),
    legend.text = element_text(size = 10, lineheight = 0.9)
  )

print(plot_obj)
# Export manually in RStudio: Plots > Export > Save as PDF.
# Verify that the PDF includes the risk table and correct baseline counts (45, 17).
writeLines(capture.output(sessionInfo()), "results/Figure_4B_sessionInfo.txt")
