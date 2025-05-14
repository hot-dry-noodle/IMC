library(CATALYST)
library(viridis)
library(dittoSeq)
library(scater)
library(patchwork)
library(cowplot)
state_markers <-c("H3K36me3","H3K4me2","H4K20me2","ubiH2AK119","H3.3","p-P53","H3K79me2","H4K20me3","H3K4me3","p-PKC","H4K16ac","pSTAT3","H4K8ac","p-GSK3","H3K9me3","p-MSK1","p-mTOR","H3K27me3","H3K27ac","p-PI3K","H4K12a","H3K9ac","H3K27me1","H3K27me2","H3")
type_markers <- c("a-SMA","OPN","CD68","CD34","CD90","vWF","RUNX2","SOX9","COL1","VIM")
rowData(spe)$marker_class <- ifelse(rownames(spe) %in% type_markers, "type",
                                    ifelse(rownames(spe) %in% state_markers, "state", 
                                           "other"))
spe_cat <- spe 
spe_cat$sample_id <- factor(spe$ROI_id)
spe_cat$condition <- factor(spe$vonkossa)
spe_cat$cluster_id <- factor(spe$celltype)
metadata(spe_cat)$cluster_codes <- data.frame(celltype = factor(spe_cat$celltype))
pbMDS(spe_cat, 
      by = "cluster_id", 
      features = rownames(spe_cat)[rowData(spe_cat)$marker_class == "state"], 
      label_by = "cluster_id", 
      k = "celltype") +
  scale_color_manual(values = metadata(spe_cat)$celltype)
