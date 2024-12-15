library(survminer)
library(survival)
library(maxstat)
library(TCGAbiolinks)
library(DESeq2)
library(TCGAbiolinksGUI.data)
library(SummarizedExperiment)
library(maftools)
library(tidyverse)
project_name <- "TCGA-LUAD"
Gene_name <- "GPS1"


lung_query <- GDCquery(project = project_name,
                       data.category = "Transcriptome Profiling",
                       experimental.strategy = "RNA-Seq",
                       workflow.type = "STAR - Counts",
                       data.type = "Gene Expression Quantification",
                       access = "open")

query2 <- GDCquery_clinic(project = project_name)
query2$deceased <- ifelse(query2$vital_status == "Alive", FALSE, TRUE)

####################################
query2$days_to_last_follow_up2 <- query2$days_to_last_follow_up/30
query2$days_to_death2 <- query2$days_to_death/30

query2$overall_survival <- ifelse(query2$vital_status == "Alive",
                                  query2$days_to_last_follow_up2,
                                  query2$days_to_death2)
################################
query2$days_to_last_follow_up3 <- query2$days_to_last_follow_up/365.25
query2$days_to_death3 <- query2$days_to_death/365.25

query2$overall_survival <- ifelse(query2$vital_status == "Alive",
                                  query2$days_to_last_follow_up3,
                                  query2$days_to_death3)

########################################
query2$overall_survival <- ifelse(query2$vital_status == "Alive",
                                  query2$days_to_last_follow_up,
                                  query2$days_to_death)
###########
GDCdownload(lung_query)
tcga_lungs_data1 <- GDCprepare(lung_query, summarizedExperiment = TRUE)

#####################
gene_id <- "ENSG00000169727.12"
filtered_tcga_lungs_data1 <- tcga_lungs_data1[rownames(tcga_lungs_data1) == gene_id, ]
print(filtered_tcga_lungs_data1)
###################

lung_matrix1 <- assay(tcga_lungs_data1,"unstranded")

gene_metadata <- as.data.frame(rowData(tcga_lungs_data1))
coldata <- as.data.frame(colData(tcga_lungs_data1))

dds1 <- DESeqDataSetFromMatrix(countData = lung_matrix1,
                               colData = coldata,
                               design = ~ 1)


keep <- rowSums(counts(dds1)) >= 10
dds1 <- dds1[keep,]

vsd1 <- vst(dds1, blind=FALSE)
lungs_matrix_vst1 <- assay(vsd1)


lung_gps1a <- lungs_matrix_vst1 %>% 
  as.data.frame() %>% 
  rownames_to_column(var = 'gene_id') %>% 
  gather(key = 'case_id', value = 'counts', -gene_id) %>% 
  left_join(., gene_metadata, by = "gene_id") %>% 
  filter(gene_name == Gene_name)


#median_value <- median(lung_gps1a$counts)
lung_gps1a$case_id2 <- gsub('-01.*', '', lung_gps1a$case_id)
lung_gps1a <- merge(lung_gps1a, query2, by.x = 'case_id2', by.y = 'submitter_id')
####
#
###


surv_data <- lung_gps1a[c("overall_survival", "deceased", "counts")]
colnames(surv_data) <- c("time", "event", "GPS1")


cutpoint_result <- surv_cutpoint(
  data = surv_data,
  time = "time",
  event = "event",
  variables = "GPS1",
  minprop = 0.1
)


optimal_cutoff <- cutpoint_result$cutpoint[1]
print(optimal_cutoff)


##
##

lung_gps1a$gene_expression_level <- ifelse(lung_gps1a$counts >= 12.13568, "HIGH", "LOW")


high <- lung_gps1a[lung_gps1a$gene_expression_level=="HIGH",]
low <- lung_gps1a[lung_gps1a$gene_expression_level=="LOW",]
fit_exp <- survfit(Surv(overall_survival,deceased)~gene_expression_level, data = lung_gps1a)
ggsurvplot(fit_exp,
           data = lung_gps1a,
           pval = T,
           risk.table = T,
           title="Yearly")

high_ids <- high$case_id
low_ids <- low$case_id

lung_gps1a$age_at_diagnosis_years <- lung_gps1a$age_at_diagnosis / 365.25



samplesHIGH <- TCGAquery_SampleTypes(
  barcode = colnames(lung_gps1a),
  typesample = c("HIGH")
)