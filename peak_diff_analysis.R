library(GenomicRanges)
library(magrittr)
mPeak = GRanges()
histL = c("H3K9ac_Pi", "H3K9ac_D")
repL = c("1","2")
for(hist in histL){
  for(rep in repL){
    peakRes = read.table(paste0("/home/wyydata/wangyuyao/CUT_TAG_Exp/peakcalling/",hist , rep,"/", hist , rep, "_seacr4_top0.01.peaks.stringent.bed"), header = FALSE, fill = TRUE)
    mPeak = GRanges(seqnames = peakRes$V1, IRanges(start = peakRes$V2, end = peakRes$V3), strand = "*") %>% append(mPeak, .)
  }
}
masterPeak = reduce(mPeak)
library(DESeq2)
for(hist in histL){
  for(rep in repL){
    countMat = matrix(NA,length(masterPeak),length(histL)*length(repL))
  }
}
library(chromVAR)
i = 1
for(hist in histL){
  for(rep in repL){
    bamFile = paste0("/home/wyydata/wangyuyao/CUT_TAG_Exp/alignment/bam4/",hist,rep,"/bowtie2.redup.bam")
    fragment_counts <- getCounts(bamFile, masterPeak, paired = TRUE, by_rg = FALSE, format = "bam")
    countMat[, i] = counts(fragment_counts)[,1]
    i = i + 1
  }
}
colnames(countMat) = c("H3K9ac_Pi1","H3K9ac_Pi2","H3K9ac_D1","H3K9ac_D2")

selectR = which(rowSums(countMat) > 5)
dataS = countMat[selectR,]
condition = factor(c(rep("Pi",2),rep("D",2)))
dds = DESeqDataSetFromMatrix(countData = dataS,
                             colData = DataFrame(condition),
                             design = ~ condition)
sizefactor <- c(1/(10000/15230),1/(10000/12891),1/(10000/17647),1/(10000/44351))
names(sizefactor) <- c("HЗK9ac_Pi1","HЗK9ac_Pi2","HЗK9ac_D1","HЗK9ac_D2")
sizeFactors(dds) <- sizefactor
DDS = DESeq(dds)
res = results(DDS)
resdata <- merge(as.data.frame(res), as.data.frame(counts(DDS, normalized=TRUE)),by="row.names",sort=FALSE)
sit <- as.data.frame(masterPeak)
resdata$anno <- paste0(sit$seqnames,"_",sit$start,"_",sit$end)


library(ChIPseeker)
library(ChIPpeakAnno)
library(GenomicFeatures)
library(org.Mm.eg.db)
library(TxDb.Mmusculus.UCSC.mm10.knownGene)
txdb <- TxDb.Mmusculus.UCSC.mm10.knownGene
peakAnno <- annotatePeak(masterPeak,
                         tssRegion = c(-3000,3000),
                         TxDb = txdb,
                         annoDb = "org.Mm.eg.db")
data_anno <- as.data.frame(peakAnno)
data_anno$anno <- paste0(data_anno$seqnames,"_",data_anno$start,"_",data_anno$end)
data <- data_anno[,c("SYMBOL","anno")]
result<- merge(resdata,data,by="anno")





