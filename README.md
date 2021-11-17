---
title: "Ex3"
author: "Maor Berkovich"
date: "11/11/2021"
output: html_document
---

```{r setup, include=FALSE}
knitr::opts_chunk$set(echo = TRUE)
```

```{r}
library("compGenomRData")
counts_file <- system.file("extdata/rna-seq/SRP029880.raw_counts.tsv",
package = "compGenomRData")
coldata_file <- system.file("extdata/rna-seq/SRP029880.colData.tsv",
package = "compGenomRData")

counts <- read.table(counts_file, header = T, sep = '\t')
#remove the 'width' column
countData <- as.matrix(subset(counts, select = c(-width)))
colData <- read.table(coldata_file, header = T, sep = '\t', stringsAsFactors = TRUE) 

# Set up a DESeqDataSet object.
library(DESeq2)
#create a DESeq dataset object from the count matrix and the colData
dds <- DESeqDataSetFromMatrix(countData = countData, colData = colData, design = ~ group)
print(dds)

#Filter out genes with low counts.
filteration <- rowSums(DESeq2::counts(dds)) >= 10
dds <- dds[filteration, ]
dds <- DESeq(dds)

#Run DESeq2 contrasting the CASE sample with CONTROL samples
DEresults = results(dds, contrast = c("group", 'CASE', 'CTRL'))
DEresultsDf <-as.data.frame(DEresults)
```
##### 1. 
```{r}
library(ggplot2)

ggplot()+
  geom_point(data = DEresultsDf, 
            mapping = aes(x = log2FoldChange, y = -log10(pvalue)), 
                   col = "grey80", size =1) +
  geom_point(data = subset(DEresultsDf, log2FoldChange > 1 & pvalue < 0.05), 
             aes(log2FoldChange, -log10(pvalue)), color = "red", size =1)+
  geom_point(data = subset(DEresultsDf, log2FoldChange < 1 & pvalue < 0.05), 
             aes(log2FoldChange, -log10(pvalue)), color = "steelblue2", size =1)+
  theme_bw()+
  theme(legend.title =element_blank())+
  labs(x = "log2 (Fold change of case vs. control)", y= "-log10(P-value)")

# The up-regulated genes in case samples are in red
#The down-regulated genes in control samples are in blue
```

##### 2.
```{r}
DESeq2::plotDispEsts(dds)
```

##### 3.
```{r}
# The default value of "lfcThreshold" argument is zero. 
DESeq2::results(dds, lfcThreshold=0)

#Changing the value to one.
# The result is that the pvalue and padj are bigger than one. 
DESeq2::results(dds, lfcThreshold=1)
```

##### 4.
```{r}
# Independent filtering is what the "results" function in DESeq2 preforms. The function use the mean of normalized counts as a filter, and use a threshold to found which optimizes the number of adjusted p values lower than a significance level alpha. The function filter the counts regardless of the source of the count i.e. independent of the group (for example CASE\CTRL). 
# The goal of independent filtering is to filter out those tests from the procedure that have no, or little chance of showing significant evidence, without even looking at their test statistic.
# If we don't use the independent filtering we cant filter out genes with little chance of showing evidence for sufficient differential expression (genes with very low counts). 
```

##### 5. 
```{r}
library(edgeR)
# Putting the data into a DGEList object:
d <- DGEList(counts=countData)
design <- model.matrix(~ group, data = colData)

# estimation of  the NB dispersion for the dataset. 
# The square root of the common dispersion gives the coefficient of variation of biological variation. 

y <- estimateDisp(d, design, robust=TRUE)
y$common.dispersion

# The dispersion estimates can be viewed in a BCV plot:
plotBCV(y)

# Determine differentiated expressed genes.
fit <- glmFit(d, design, dispersion=y$common.dispersion)
lrt <- glmLRT(fit)
topTags

# Plot log-fold change against log-counts per million.
plotMD(lrt)
abline(h=c(-1, 1), col="blue")

# compare the results between DESeq2 vs edgeR
summary(decideTests(lrt))
summary(DEresults)
```

##### 6. 
```{r}
library(compcodeR)
```