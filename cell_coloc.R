roidis <- data.frame()
prop <- prop.table(table(spe$celltype,spe$ROI_id), margin = 2)
for(q in c(unique(spe$ROI_id))){
spe2 <- spe[,spe$ROI_id == q ]
if ( prop[5,q] > 0.01& prop[6,q] > 0.005){
speMSC <- spe2[,spe2$celltype == "MSC"]
speosteo <- spe2[,spe2$celltype == "osteo-like VSMC"]
meandis <- data.frame()
for (i in c(1:nrow(spatialCoords(speMSC)))){
  dis <- vector()
  for (p in c(1:nrow(spatialCoords(speosteo)))){
    dis2 <- (sqrt(((spatialCoords(speMSC)[i,1]-spatialCoords(speosteo)[p,1])^2)+((spatialCoords(speMSC)[i,2]-spatialCoords(speosteo)[p,2])^2)))
    dis <- c(dis,dis2)
    }
    meandis2 <- mean(dis)
    meandis <- rbind(meandis,meandis2)
}
colnames(meandis) <- q
roidis <- c(roidis,meandis)
}
}


roidis <- data.frame()
prop <- prop.table(table(spe$celltype,spe$ROI_id), margin = 2)
for(q in c(unique(spe$ROI_id))){
  spe2 <- spe[,spe$ROI_id == q ]
  if ( prop[5,q] > 0.01& prop[6,q] > 0.005){
    speMSC <- spe2[,spe2$celltype == "MSC"]
    speother <- spe2[,spe2$celltype != "MSC"]
    speother <- speother[,speother$celltype != "osteo-like VSMC"]
    if (ncol(speother) > 500 ) {
    cur_cells <- sample(seq_len(ncol(speother)), 500)
    speother <- speother[,cur_cells]
    }
    meandis <- data.frame()
    for (i in c(1:nrow(spatialCoords(speMSC)))){
      dis <- vector()
      for (p in c(1:nrow(spatialCoords(speother)))){
        dis2 <- (sqrt(((spatialCoords(speMSC)[i,1]-spatialCoords(speother)[p,1])^2)+((spatialCoords(speMSC)[i,2]-spatialCoords(speother)[p,2])^2)))
        dis <- c(dis,dis2)
      }
      meandis2 <- mean(dis)
      meandis <- rbind(meandis,meandis2)
    }
    colnames(meandis) <- q
    roidis <- c(roidis,meandis)
  }
}

group <- vector()
for (i in names(roidis)){
  if (length(roidis[[i]])> 0 ){
  data <- rbind(cbind(roidismscos[[i]],"a"),cbind(roidis[[i]],"b"))
  data <- as.data.frame(data)
  data$V1 <- as.numeric(data$V1)
  p <- wilcox.test(data$V1~data$V2,paired = FALSE)
  file_name <- paste0(i,".pdf")
  if (p$p.value <= 0.05){
    if (p$p.value <= 0.001){
    ggplot(data=data, aes(x=data$V1, color=data$V2)) +
      geom_density(adjust=1.5, alpha=1,lwd=1.5,linetype = 1) + labs(y='Density', x='Mean distance to osteo-like cell(pixels)')+ theme_classic()+theme(legend.position = "none",text = element_text(size = 15)) + ggtitle(str_wrap(paste0("MSC and osteo-like cell co-localized,",i,",p < 0.001")))
  ggsave(file_name, path = "/Users/wangyuyao/pic",dpi = 300,device = "pdf",width = 8,height = 5)
    pval <- "co-localized"
    names(pval) <- i
    group <- c(group,pval)
    }
    else {
      ggplot(data=data, aes(x=data$V1, color=data$V2)) +
        geom_density(adjust=1.5, alpha=1,lwd=1.5,linetype = 1) + labs(y='Density', x='Mean distance to osteo-like cell(pixels)')+ theme_classic()+theme(legend.position = "none",text = element_text(size = 15)) + ggtitle(str_wrap(paste0("MSC and osteo-like cell co-localized,",i,",p = ",round(p$p.value,3))))
      ggsave(file_name, path = "/Users/wangyuyao/pic",dpi = 300,device = "pdf",width = 8,height = 5)
      pval <- "co-localized"
      names(pval) <- i
      group <- c(group,pval)
    }
    }
  else {
    ggplot(data=data, aes(x=data$V1, color=data$V2)) +
      geom_density(adjust=1.5, alpha=1,lwd=1.5,linetype = 1) + labs(y='Density', x='Mean distance to osteo-like cell(pixels)')+ theme_classic()+theme(legend.position = "none",text = element_text(size = 15)) + ggtitle(str_wrap(paste0("MSC and osteo-like cell NOT co-localized,",i,",p = ",round(p$p.value,3))))
      ggsave(file_name, path = "/Users/wangyuyao/pic",dpi = 300,device = "pdf",width = 8,height = 5)
      pval <- "non co-localized"
      names(pval) <- i
      group <- c(group,pval)
  }
  }
}
group <- as.data.frame(group)
group$msc <- prop2$Freq[match(prop2$Var2,rownames(group))]
group$id <- rownames(group)
ggplot(group,aes(reorder(id,-msc),msc)) + geom_bar(stat = 'identity',aes(fill=group))
