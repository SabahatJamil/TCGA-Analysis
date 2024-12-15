library(survminer)
library(survival)
library(maxstat)
library(TCGAbiolinks)
library(DESeq2)
library(TCGAbiolinksGUI.data)
library(SummarizedExperiment)
library(maftools)
library(tidyverse)



survival_anlaysis_func<- function(project_name,Gene_name){
  
  lung_query <- GDCquery(project = project_name,
                         data.category = "Transcriptome Profiling",
                         experimental.strategy = "RNA-Seq",
                         workflow.type = "STAR - Counts",
                         data.type = "Gene Expression Quantification",
                         access = "open")
  
  query2 <- GDCquery_clinic(project = project_name)
  query2$deceased <- ifelse(query2$vital_status == "Alive", FALSE, TRUE)
  
  query2$overall_survival <- ifelse(query2$vital_status == "Alive",
                                    query2$days_to_last_follow_up,
                                    query2$days_to_death)
  
  query2$days_to_last_follow_up_montly <- query2$days_to_last_follow_up/30
  query2$days_to_death_montly <- query2$days_to_death/30
  
  query2$overall_survival_montly <- ifelse(query2$vital_status == "Alive",
                                           query2$days_to_last_follow_up_montly,
                                           query2$days_to_death_montly)
  ################################
  query2$days_to_last_follow_up_yearly <- query2$days_to_last_follow_up/365.25
  query2$days_to_death_yearly <- query2$days_to_death/365.25
  
  query2$overall_survival_yearly <- ifelse(query2$vital_status == "Alive",
                                           query2$days_to_last_follow_up_yearly,
                                           query2$days_to_death_yearly)
  
  
  GDCdownload(lung_query)
  tcga_lungs_data1 <- GDCprepare(lung_query, summarizedExperiment = TRUE)
  
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
  
  surv_data <- lung_gps1a[c("overall_survival", "deceased", "counts")]
  colnames(surv_data) <- c("time", "event", Gene_name)
  
  
  cutpoint_result <- surv_cutpoint(
    data = surv_data,
    time = "time",
    event = "event",
    variables = Gene_name,
    minprop = 0.1
  )
  
  
  optimal_cutoff <- cutpoint_result$cutpoint[1]
  optimal_cutoff <- as.numeric(optimal_cutoff)
  print(optimal_cutoff)
  
  
  lung_gps1a$gene_expression_level <- ifelse(lung_gps1a$counts >= optimal_cutoff, "HIGH", "LOW")
  
  fit_exp <- survfit(Surv(overall_survival,deceased)~gene_expression_level, data = lung_gps1a)
  p1 <- ggsurvplot(fit_exp,
                   data = lung_gps1a,
                   pval = T,
                   risk.table = T,
                   title="Daily")
  print(p1)
  
  #@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
  surv_data<- lung_gps1a[c("overall_survival_montly", "deceased", "counts")]
  colnames(surv_data) <- c("time", "event", Gene_name)
  
  
  cutpoint_result <- surv_cutpoint(
    data = surv_data,
    time = "time",
    event = "event",
    variables = Gene_name,
    minprop = 0.1
  )
  
  
  optimal_cutoff <- cutpoint_result$cutpoint[1]
  optimal_cutoff <- as.numeric(optimal_cutoff)
  
  
  lung_gps1a$gene_expression_level <- ifelse(lung_gps1a$counts >= optimal_cutoff, "HIGH", "LOW")
  
  fit_exp <- survfit(Surv(overall_survival_montly,deceased)~gene_expression_level, data = lung_gps1a)
  pm <- ggsurvplot(fit_exp,
                   data = lung_gps1a,
                   pval = T,
                   risk.table = T,
                   title="Monthly")
  print(pm)
  #@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
  surv_data<- lung_gps1a[c("overall_survival_yearly", "deceased", "counts")]
  colnames(surv_data) <- c("time", "event", Gene_name)
  
  
  cutpoint_result <- surv_cutpoint(
    data = surv_data,
    time = "time",
    event = "event",
    variables = Gene_name,
    minprop = 0.1
  )
  
  
  optimal_cutoff <- cutpoint_result$cutpoint[1]
  optimal_cutoff <- as.numeric(optimal_cutoff)
  
  
  lung_gps1a$gene_expression_level <- ifelse(lung_gps1a$counts >= optimal_cutoff, "HIGH", "LOW")
  
  fit_exp <- survfit(Surv(overall_survival_yearly,deceased)~gene_expression_level, data = lung_gps1a)
  py <- ggsurvplot(fit_exp,
                   data = lung_gps1a,
                   pval = T,
                   risk.table = T,
                   title="Yearly")
  print(py)
}

survival_anlaysis_func(project_name = "TCGA-COAD", Gene_name = "GPS1")
