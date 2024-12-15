# Load required packages
library(TCGAbiolinks)
library(limma)
library(edgeR)
library(dplyr)
library(DESeq2)


BiocManager::install(c("locfit", "Rcpp"))
BiocManager::install("biomaRt")
library(biomaRt)

install.packages("BiocManager") 
BiocManager::install("edgeR")

query <- GDCquery(project = "TCGA-LUAD", 
                  data.category = "Transcriptome Profiling",
                  data.type = "Gene Expression Quantification",
                  workflow.type = "STAR - Counts")

GDCdownload(query)

brca.exp <- GDCprepare(query)

write.table(assay(brca.exp), file="tcga_luad_gene_expression.txt", quote=FALSE, sep="\t")
counts <- read.delim("tcga_luad_gene_expression.txt", row.names=1)
counts <- counts[!duplicated(rownames(counts)),]


####################################################

filtered_counts <- counts[grepl("^ENSG00000169727", rownames(counts)), ]

# Check the filtered data
print(filtered_counts)

###########################
ensembl <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")
gene_info <- getBM(attributes = c("ensembl_gene_id", "hgnc_symbol"), 
                   filters = "ensembl_gene_id", 
                   values = rownames(counts), 
                   mart = ensembl)
gene_info <- gene_info[!duplicated(gene_info$ensembl_gene_id), ]

counts <- counts[rownames(counts) %in% gene_info$ensembl_gene_id, ]
rownames(counts) <- gene_info$hgnc_symbol[match(rownames(counts), gene_info$ensembl_gene_id)]
counts <- counts[!is.na(rownames(counts)),]

###############################
d <- DGEList(counts=filtered_counts)
d <- DGEList(counts)
keep <- filterByExpr(d)
d <- d[keep,,keep.lib.sizes=FALSE] 
d <- calcNormFactors(d)

# Divide samples into GPS1 high and low groups
gps1_expr <- d$samples
gps1_high <- gps1_expr[gps1_expr$lib.size > median(gps1_expr$lib.size),]
gps1_low <- gps1_expr[gps1_expr$lib.size <= median(gps1_expr$lib.size),]
group <- factor(c(rep("high", length(gps1_high)), rep("low", length(gps1_low))))


d$samples$gene_expression_level <- ifelse(d$samples$lib.size >= median(d$samples$lib.size), "HIGH", "LOW")

# Divide samples into GPS1 high and low groups on the basis of cutoff
gps1_expr <- d$samples$lib.size
#lung_gps$counts >= 12.13568
gps1_high <- names(gps1_expr[gps1_expr > 12.13568 ])
gps1_low <- names(gps1_expr[gps1_expr <= 12.13568 ])
group <- factor(c(rep("high", length(gps1_high)), rep("low", length(gps1_low))))

design <- model.matrix(~0+group)
colnames(design) <- levels(group)

d <- estimateDisp(d, design)
fit <- glmFit(d, design)
lrt <- glmLRT(fit, coef=2)
topTags <- topTags(lrt, n=Inf)$table

deg_list <- topTags %>%
  as.data.frame() %>%
  rownames_to_column(var="gene") %>%
  filter(abs(logFC) > 1, FDR < 0.05)

volcanoplot(deg_list$logFC, -log10(deg_list$FDR), main="Volcano Plot", pch=20, highlight=rownames(deg_list))

top_degs <- rownames(deg_list)[1:50]
heatmap_data <- d$pseudo[top_degs,]
heatmap_data <- heatmap_data[,c(gps1_high, gps1_low)]
heatmap(heatmap_data, scale="row", labRow=top_degs)
library(clusterProfiler)
ego <- enrichGO(gene=rownames(deg_list),
                OrgDb=org.Hs.eg.db,
                keyType="SYMBOL",
                ont="BP",
                pAdjustMethod="BH",
                qvalueCutoff=0.05)