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
query <- GDCquery(project = 'TCGA-KIRP', #Lung squamous cell carcinoma, getGDCprojects()$project_id for full option list
                  data.category = 'Transcriptome Profiling',
                  experimental.strategy = 'RNA-Seq',
                  workflow.type = 'STAR - Counts',
                  access = 'open',
                  data.type = 'Gene Expression Quantification')

# View sample metadata before downloading
output_query <- getResults(query)
dim(output_query)

# download data
GDCdownload(query, files.per.chunk = 165)

# prepare data
## run to save prepared data as TCGA_LGG.rda
kirp_data <- GDCprepare(query, 
                        summarizedExperiment = TRUE, #Default
                        save=TRUE, 
                        save.filename = 'TCGA_KIRP.rda')
View(kirp_data)

preprocessed_kirp <- TCGAanalyze_Preprocessing(
  kirp_data,
  cor.cut = 0.6,
  datatype = "unstranded",   # Make sure this matches your RNA-seq data type
  filename = "KIRP_preprocessing_result.png"
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

#Selecting 48 random tumor samples
tumor_tissue <- tumor_tissue[sample(nrow(tumor_tissue), 48), ]

#Binding tumor and normal samples
working_sample_info <- rbind(normal_tissue, tumor_tissue)

#Making the samples (columns) of expression data align with samples (rows) of gene data
working_sample_info <- working_sample_info %>%
  rownames_to_column(var = "number_id") %>% #Converts numbered rownames to column
  mutate(barcode_copy = barcode) %>%       # Create a new column to keep barcodes
  column_to_rownames(var = "barcode") %>%  # Move original 'barcode' to rownames
  select(-number_id) #Remove the newly formed numbered column

#Select the desired samples based on the selected samples for working_sample_info
preprocessed_kirp <- preprocessed_kirp[, rownames(working_sample_info)]
dim(preprocessed_kirp)

#Checking if matching was done correctly
table(rownames(working_sample_info) %in% colnames(preprocessed_kirp))
table(rownames(working_sample_info) == colnames(preprocessed_kirp))

#Filter the lowly expressed genes using rowMeans
preprocessed_kirp <- preprocessed_kirp[rowMeans(preprocessed_kirp) >= 10,]
dim(preprocessed_kirp)

#Matched gene names in gene_info, with expression data
gene_info <- gene_info[rownames(gene_info) %in% rownames(preprocessed_kirp), ]

#Checking if matching was done correctly
table(rownames(gene_info) %in% rownames(preprocessed_kirp))
table(rownames(gene_info) == rownames(preprocessed_kirp))

#Save your data to files
write.csv(gene_info, 'gene_info.csv')
write.csv(working_sample_info, 'working_sample_info.csv')
write.csv(preprocessed_kirp, 'preprocessed_kirp')


















































































