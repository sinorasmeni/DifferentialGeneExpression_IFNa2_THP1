library(edgeR)
library(ggplot2)
library(ggrepel)



#list of paths to each sample
sample_folders <- c("path/to/salmon_outputs/IFNa_rep1/",
                    "path/to/salmon_outputs/IFNa_rep2/",
                    "path/to/salmon_outputs/IFNa_rep3/",
                    "path/to/salmon_outputs/cnt_rep1/",
                    "path/to/salmon_outputs/cnt_rep2/",
                    "path/to/salmon_outputs/cnt_rep3/")

files <- file.path(sample_folders, "quant.sf")

#Import Salmon quants using tximport
txi_det <- tximport(files, type = "salmon", txOut = TRUE, ignoreTxVersion = TRUE)
dte <- DGEList(counts = txi_det$counts)


#metadata
sample_metadata <- data.frame(
  sample_id = c("trim2transcripts_quant2","trim2transcripts_quant3", "trim2transcripts_quant4", "trim2transcripts_quant5", "trim2transcripts_quant6", "trim2transcripts_quant7"),
  group = c("IFNa", "IFNa", "IFNa", "cnt", "cnt", "cnt")
)

dte$samples$group <- sample_metadata$group



#filter lowly expressed transcripts
keep <- filterByExpr(dte, group = dte$samples$group)
dte <- dte[keep, , keep.lib.sizes = FALSE]



#normalise for library sizes
dte <- calcNormFactors(dte)

#design matrix based on condition
design <- model.matrix(~ group, data = sample_metadata)

#estimate dispersion and fit model
dte <- estimateDisp(dte, design)
fit_det <- glmFit(dte, design)

#Differential transcript expression testing 
lrt_det <- glmLRT(fit_det, coef = 2)


topTags(lrt_det)

lrt_det$table$FDR <- p.adjust(lrt_det$table$PValue, method = "BH")

result_det <- topTags(lrt_det, n = Inf)

head(result_det$table)

det_results <- result_det$table

# Assuming lrt$table contains the results with logFC and p-value columns
volcano_data_det <- data.frame(
  logFC = lrt_det$table$logFC,
  pValue = -log10(lrt_det$table$PValue),
  Gene = rownames(lrt_det$table)
)



# Define significance thresholds
significance_threshold <- 0.05
log2fc_threshold <- 1  # Corresponds to a fold change of 2

volcano_data_det$Significant <- ifelse(
  volcano_data_det$pValue > -log10(significance_threshold) & abs(volcano_data_det$logFC) > log2fc_threshold,
  ifelse(volcano_data_det$logFC > 0, 'Upregulated', 'Downregulated'),
  'Not Significant'
)

# Plot
ggplot(volcano_data_det, aes(x = logFC, y = pValue, color = Significant)) +
  geom_point(alpha = 0.7) +
  scale_color_manual(values = c('Upregulated' = 'red', 'Downregulated' = 'blue', 'Not Significant' = 'grey')) +
  geom_vline(xintercept = c(-log2fc_threshold, log2fc_threshold), linetype = "dashed", color = "black") +  # Vertical lines for log2 fold change thresholds
  geom_hline(yintercept = -log10(significance_threshold), linetype = "dashed", color = "black") +  # Horizontal line for p-value threshold
  xlab("Log2 Fold Change") +
  ylab("-Log10(P-Value)") +
  ggtitle("Volcano Plot") +
  theme_minimal() +
  theme(legend.title = element_blank())  # Remove legend title
