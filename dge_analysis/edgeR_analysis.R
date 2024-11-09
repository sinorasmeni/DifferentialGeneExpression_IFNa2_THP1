library(edgeR)
library(ggplot2)
library(ggrepel)
library(tximport)
library(org.Hs.eg.db)
library(tidyverse)

###Differential transcript expression

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


### Differential gene expression

#Path to quant files
files <- file.path(sample_folders, "quant.sf")

#Path to GTF file
gtf_file <- "path/to/gtf_file"

#Create a TxDb object from the GTF file
txdb <- makeTxDbFromGFF(gtf_file, format = "gtf")

# Extract transcript-to-gene mapping (tx2gene) using transcripts() and gene mapping from the TxDb object
tx2gene <- transcripts(txdb, columns = c("TXNAME", "GENEID"))

# Convert to data frame for easier manipulation
tx2gene_df <- as.data.frame(tx2gene) %>%
  select(TXNAME, GENEID)  # Select only relevant columns

# Check the tx2gene mapping
head(tx2gene_df)

# Convert columns to character vectors if necessary
tx2gene_df$TXNAME <- as.character(tx2gene_df$TXNAME)
tx2gene_df$GENEID <- as.character(tx2gene_df$GENEID)

tx2gene_df$TXNAME <- gsub("\\..*$", "", tx2gene_df$TXNAME)
# Use tximport to read Salmon quantifications and summarize to the gene level
txi <- tximport(files, type = "salmon", tx2gene = tx2gene_df, ignoreTxVersion = TRUE)

# Create a DGEList object from the txi counts
dge <- DGEList(counts = txi$counts)

#Metadata
sample_metadata <- data.frame(sample_id = c("IFNa_rep1", "IFNa_rep2", "IFNa_rep3", 
                                            "cnt_rep1", "cnt_rep2", "cnt_rep3"),
                              group = c("IFNa", "IFNa", "IFNa", "cnt", "cnt", "cnt")
)

# Ensure the sample names in dge match those in sample_metadata
rownames(dge$samples) <- sample_metadata$sample_id
dge$samples$group <- sample_metadata$group

keep <- filterByExpr(dge, group = dge$samples$group)
dge <- dge[keep, , keep.lib.sizes = FALSE]

# Normalize library sizes
dge <- calcNormFactors(dge)

# Design matrix
design <- model.matrix(~ group, data = dge$samples)

# Estimate dispersion
dge <- estimateDisp(dge, design)

# Fit the model and perform differential expression testing
fit <- glmFit(dge, design)
lrt <- glmLRT(fit, coef = 2) 

# View top tags
top_results <- topTags(lrt, n = Inf)
length(top_results)

# Create volcano plot data
volcano_data <- data.frame(
  logFC = top_results$table$logFC,
  pValue = -log10(top_results$table$PValue),
  Gene = rownames(top_results$table)
)

genes <-  rownames(top_results$table)
# Define significance thresholds
significance_threshold <- 0.05
log2_fold_change_threshold <- 1  # log2 fold change threshold of 1

# Determine significance
volcano_data$Significant <- ifelse(
  volcano_data$pValue > -log10(significance_threshold) & abs(volcano_data$logFC) > log2_fold_change_threshold,
  ifelse(volcano_data$logFC > 0, 'Upregulated', 'Downregulated'),
  'Not Significant'
)

# Count the number of genes in each category
significance_counts <- table(volcano_data$Significant)

# Print the counts
print(significance_counts)

# Plot with labels for top genes
ggplot(volcano_data, aes(x = logFC, y = pValue, color = Significant)) +
  geom_point(alpha = 0.7) +
  scale_color_manual(values = c('Upregulated' = 'purple', 'Downregulated' = 'skyblue', 'Not Significant' = 'grey')) +
  xlab("Log2 Fold Change") +
  ylab("-Log10(P-Value)") +
  ggtitle("Volcano Plot") +
  theme_minimal() +
  # Add labels for top 10 genes
  geom_label_repel(data = top_genes, aes(label = gene_symbol),
                   size = 3, box.padding = 0.1, point.padding = 0.1,
                   max.overlaps = Inf)
