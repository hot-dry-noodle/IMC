library(ReactomePA)
library(tidyverse)
library(data.table)
library(org.Mm.eg.db)
library(clusterProfiler)
library(biomaRt)
library(enrichplot)
library(openxlsx)
library(DESeq2)
#readtable
base <- read.xlsx("IMC/rnaseq/gene_count_exp.annotation.xlsx")
base = base[!duplicated(base$gene_name),]
rownames(base) <- base$gene_name
base <- base[,-1]
base <- base[,-c(10:16)]
colnames(base) <- c("C1","C2","C3","P1","P2","P3","D1","D2","D3")
database <- base[,c(1:6)]
condition <- factor(c(rep("C",3),rep("P",3)))
coldata <- data.frame(row.names = colnames(database), condition)
coldata    #显示coldata值,看看分组信息与真实数据是否一致，很关键！！！！
dds <- DESeqDataSetFromMatrix(countData=database, colData=coldata, design=~condition)
keep <- rowSums(counts(dds) >= 10) >= 3  #过滤低表达基因，至少有3个样品都满足10个以上的reads数  
dds <- dds[keep, ] 
dds1 <- DESeq(dds)    # 将数据标准化，必要步骤！！！
resultsNames(dds1)    # 查看结果的名称。
dds1$condition        #默认后者的处理组比前面的对照组。
res <- results(dds1)  # 必要步骤！！！
summary(res)          #看一下结果的概要信息，p值默认小于0.1。
table(res$padj < 0.05)        #padj 即矫正后的P值。看看有多少差异基因满足所设的P值要求。TRUE的数值为满足要求的基因个数。
res <- res[order(res$padj),]  #按照padj 进行升序排列
resdata <- merge(as.data.frame(res), as.data.frame(counts(dds1, normalized=TRUE)),by="row.names",sort=FALSE)
resdata$group=ifelse(resdata$pvalue>0.01,'stable',
                     ifelse( resdata$log2FoldChange >0.8,'up',
                     ifelse( resdata$log2FoldChange < -0.8,'down','stable') ))
resdata <- resdata[which(resdata$group != "stable"),]

resdata <- resdata[order(resdata$log2FoldChange,decreasing = TRUE),]
geneList <- resdata[,3]
names(geneList) <- resdata[,1]
Go_gseresult <- gseGO(geneList, 'org.Mm.eg.db', keyType = "SYMBOL", ont="BP", nPerm = 1000, minGSSize = 10, maxGSSize = 1000, pvalueCutoff=1)
View(Go_gseresult@result)








