library(TCGAbiolinks)
library(SummarizedExperiment)
library(edgeR)
library(ggplot2)
library(dplyr)
library(stringr)
library(magrittr)
library(tidyr)
library(readr)
library(tidyverse)


# build a query to retrieve DNA expression data
query <- GDCquery(project = 'TCGA-KIRP', #Kidney renal papillary cell carcinoma, use getGDCprojects()$project_id for full option list
                  data.category = 'Transcriptome Profiling',
                  experimental.strategy = 'RNA-Seq',
                  workflow.type = 'STAR - Counts',
                  access = 'open',
                  data.type = 'Gene Expression Quantification')

# View sample metadata before downloading (Optional)
output_query <- getResults(query)
dim(output_query)

# download data
GDCdownload(query, files.per.chunk = 165)

# prepare data
## run to save prepared data as TCGA_KIRP.rda
kirp_data <- GDCprepare(query, 
                        summarizedExperiment = TRUE, #Default
                        save=TRUE, 
                        save.filename = 'TCGA_KIRP.rda')
View(kirp_data)

#Preprocessing gives the expression data matrix.
#This isn't required if DESeq2 or limma are being used for the DE analysis as those require raw counts
preprocessed_kirp <- TCGAanalyze_Preprocessing(
  kirp_data,
  cor.cut = 0.6, #Samples with correlation coefficients of less than 0.6 are considered otliers and removed
  datatype = "unstranded",   # Make sure this matches your RNA-seq data type
  filename = "KIRP_preprocessing_result2.png"
)

# normalising by gene length and GC content
preprocessed_kirp <- TCGAanalyze_Normalization(preprocessed_kirp, 
                                               geneInfo = geneInfoHT, 
                                               method='geneLength')

preprocessed_kirp <- TCGAanalyze_Normalization(preprocessed_kirp, 
                                               geneInfo = geneInfoHT, 
                                               method='gcContent')

#Extract the sample_info and gene_info
#Gene_info
gene_info <- rowData(kirp_data) %>% as.data.frame()
View(gene_info)

#Sample_info
sample_info <- colData(kirp_data) |>
  as.data.frame() |>
  select(-treatments) |>
  apply(2, as.character) |>
  as.data.frame()
View(sample_info)

#Subsetting sample_info and removing rows with NA
sample_info <- sample_info %>% 
  select(barcode,
         patient,
         sample_type,
         tissue_type,
         tumor_descriptor,
         ajcc_pathologic_stage,
         vital_status,
         gender,
         race,
         ethnicity,
         age_at_index,
         days_to_death,
         days_to_last_follow_up) %>% 
  filter(complete.cases(.)) %>% #Remove NA rows
  mutate( #Change pathological stage to NA if tissue type is normal
    ajcc_pathologic_stage = ifelse(
      tissue_type == "Normal",
      "NA",
      ajcc_pathologic_stage
    )
  )

#Checking number of tumor and normal samples
sum(sample_info$tissue_type == "Normal")
sum(sample_info$tissue_type == "Tumor")

#Using a 2:1 tumor to sample ratio to reduce bias while keeping statistical power
normal_tissue <- sample_info[sample_info$tissue_type == "Normal", ]
tumor_tissue <- sample_info[sample_info$tissue_type == "Tumor", ]

#Selecting 48 random tumor samples with sample() function
tumor_tissue <- tumor_tissue[sample(nrow(tumor_tissue), 48), ]

#Binding tumor and normal samples
working_sample_info <- rbind(normal_tissue, tumor_tissue)

#Making the samples (columns) of expression data align with samples (rows) of sample_info
working_sample_info <- working_sample_info %>%
  rownames_to_column(var = "number_id") %>% #Converts numbered rownames to column
  mutate(barcode_copy = barcode) %>%       # Create a new column to keep barcodes
  column_to_rownames(var = "barcode") %>%  # Move original 'barcode' to rownames
  select(-number_id) #Remove the newly formed numbered column

#Select the desired samples based on the selected samples for working_sample_info
preprocessed_kirp <- preprocessed_kirp[, rownames(working_sample_info)]
dim(preprocessed_kirp)

#Checking if matching was done correctly and in the right order
table(rownames(working_sample_info) %in% colnames(preprocessed_kirp))
table(rownames(working_sample_info) == colnames(preprocessed_kirp))

#Filter the lowly expressed genes using rowMeans
preprocessed_kirp <- preprocessed_kirp[rowMeans(preprocessed_kirp) >= 10,]
dim(preprocessed_kirp)

#Matched gene names in gene_info, with expression data
gene_info <- gene_info[rownames(gene_info) == rownames(preprocessed_kirp), ]

#Save your data to files
write.csv(gene_info, 'gene_info.csv')
write.csv(working_sample_info, 'working_sample_info.csv')
write.csv(preprocessed_kirp, 'preprocessed_kirp')

#Read back your data files as variables in R
gene_info <- read.csv('gene_info.csv', row.names = 1)
working_sample_info <- read.csv('working_sample_info.csv', row.names = 1)
preprocessed_kirp <- read.csv('preprocessed_kirp', row.names = 1)

##########################################################################
################# DIFFERENTIAL EXPRESSION ANALYSIS #######################
##########################################################################
#Set my group to store each sample's tissue type
group <- working_sample_info$tissue_type

#Carry out my DEA based on my set groups
de_results <- TCGAanalyze_DEA(
  mat1 = preprocessed_kirp[, group == "Normal"],
  mat2 = preprocessed_kirp[, group == "Tumor"],
  Cond1type = "Normal",
  Cond2type = "Tumor"
)

#Viewing result statistics and visualizing
head(de_results)
summary(de_results)

p <- TCGAVisualize_volcano(de_results$logFC,
                      de_results$FDR,
                      names(de_results$logFC),
                      x.cut = 1, 
                      y.cut = 0.01,
                      xlab = "log2 Fold Change",
                      ylab = "-log10 adjusted p value",
                      title = "Volcano plot for DE genes")
ggsave("tcga_volcano_plot.png", plot = p, width = 15, height = 10, dpi = 300)

#Subsetting my significantly DEGs based on logFC and FDR
sig_de_results <- de_results[abs(de_results$logFC) > 2 & de_results$FDR < 0.05,]
write.csv(sig_de_results, "sig_de_results.csv")

#Get your top upregulated and downregulated genes using dplyr
top_25_upregulated <- sig_de_results %>% 
  arrange(desc(logFC)) %>% 
  head(., 25)
  
top_25_downregulated <- sig_de_results %>% 
  arrange(logFC) %>% 
  head(., 25)

#Prepare data for counts and values of upregulated genes
#Select the top 25 upregulated genes and change to long format so you can left_join with sample_info
upregulated_plot_data <- preprocessed_kirp[rownames(top_25_upregulated), ] %>% 
  as.data.frame() %>% 
  rownames_to_column("genes") %>% 
  pivot_longer(cols = -genes,
               names_to = "samples",
               values_to = "count") %>% 
  left_join(., working_sample_info, by = c("samples" = "barcode_copy"))

#Add log2FC and FDR values from DE results
upregulated_plot_data <-  as.data.frame(upregulated_plot_data) %>% 
  left_join(., sig_de_results %>% 
              rownames_to_column(var = "genes") %>% 
              select(genes, logFC, FDR),
            by = "genes")

# Create the plot
upregulated_plot <-  ggplot(upregulated_plot_data, aes(x = tissue_type, y = log2(count + 1), fill = tissue_type)) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(width = 0.2, size = 0.5, alpha = 0.5) +
  facet_wrap(~ genes, scales = "free_y", ncol = 5) +
  scale_fill_manual(values = c("Tumor" = "aquamarine3", "Normal" = "yellow")) +
  labs(y = "log2CPM", x = NULL) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        strip.background = element_blank(),
        strip.text = element_text(face = "bold"),
        legend.position = "none") +
  geom_text(data = upregulated_plot_data %>% group_by(genes) %>% slice(1),
            aes(x = tissue_type[1], y = Inf, label = sprintf("q = %.2e", FDR)),
            vjust = 1.25, size = 3, inherit.aes = FALSE)


#Prepare data for counts and values of downregulated genes, manipulate and create the plot just as for upregulated genes
downregulated_plot_data <- preprocessed_kirp[rownames(top_25_downregulated), ] %>% 
  as.data.frame() %>% 
  rownames_to_column("genes") %>% 
  pivot_longer(cols = -genes,
               names_to = "samples",
               values_to = "count") %>% 
  left_join(., working_sample_info, by = c("samples" = "barcode_copy"))

#Add log2FC and FDR values
downregulated_plot_data <-  as.data.frame(downregulated_plot_data) %>% 
  left_join(., sig_de_results %>% 
              rownames_to_column(var = "genes") %>% 
              select(genes, logFC, FDR),
            by = "genes")

# Create the plot
downregulated_plot <- ggplot(downregulated_plot_data, aes(x = tissue_type, y = log2(count + 1), fill = tissue_type)) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(width = 0.2, size = 0.5, alpha = 0.5) +
  facet_wrap(~ genes, scales = "free_y", ncol = 5) +
  scale_fill_manual(values = c("Tumor" = "lightblue2", "Normal" = "lightpink")) +
  labs(y = "log2CPM", x = NULL) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        strip.background = element_blank(),
        strip.text = element_text(face = "bold"),
        legend.position = "none") +
  geom_text(data = downregulated_plot_data %>% group_by(genes) %>% slice(1),
            aes(x = tissue_type[1], y = Inf, label = sprintf("q = %.2e", FDR)),
            vjust = 1.25, size = 3, inherit.aes = FALSE)

#Save your plots and data
write.csv(upregulated_plot_data, "upregulated_plot_data.csv")
write.csv(downregulated_plot_data, "downregulated_plot_data.csv")
ggsave("top_25_upregulated_genes.png", plot = upregulated_plot, width = 15, height = 10, dpi = 300)
ggsave("top_25_downregulated_genes.png", plot = downregulated_plot, width = 15, height = 10, dpi = 300)

###################################################################
#################### ENRICHMENT ANALYSIS ##########################
###################################################################
library(clusterProfiler)
library(org.Hs.eg.db) # Human gene annotation database

# Extract gene IDs
genes <- rownames(sig_de_results)

# Convert Ensembl IDs to Entrez IDs as enrichKEGG only accepts Entrez gene IDs
gene_conversion <- bitr(genes, fromType = "ENSEMBL", toType = "ENTREZID", OrgDb = org.Hs.eg.db)

# Create gene list for enrichment
entrez_gene_list <- gene_conversion$ENTREZID

# Run KEGG enrichment
kegg_enrichment <- enrichKEGG(
  gene = entrez_gene_list,
  organism = 'hsa', # hsa = Homo sapiens
  pvalueCutoff = 0.05
)

# Check the top enriched pathways
head(kegg_enrichment)

# Save the result to CSV
write.csv(as.data.frame(kegg_enrichment), "KEGG_enrichment_results.csv")

# Barplot
kegg_barplot <- barplot(kegg_enrichment, showCategory = 10)
ggsave("KEGG_Enrichment_Barplot.png", plot = kegg_barplot, width = 15, height = 10, dpi = 300)

# Dotplot
kegg_dotplot <-  dotplot(kegg_enrichment, showCategory = 10)
ggsave("KEGG_Enrichment_Dotplot.png", plot = kegg_dotplot, width = 15, height = 10, dpi = 300)





















































































