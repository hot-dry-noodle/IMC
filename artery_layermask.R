library(raster)
#readthepicture
file <- list.files("/Users/wangyuyao/hyperion/mask_tiff")
spe3 <- spe[,spe$batch == "nan"]
for (i in 1:length(file)){
path <- paste("/Users/wangyuyao/hyperion/mask_tiff/",file[i],sep = "")
test <- raster(path)
test_correct <- calc(test,function(x) x)
test2 <- matrix(test_correct@data@values,ncol = test_correct@ncols,nrow = test_correct@nrows,byrow = TRUE )
spe2 <- spe[,spe$file_id == file[i]]
spatialCoords(spe2)
area <- as.data.frame(spatialCoords(spe2))
qqq <- vector()
for (s in 1:nrow(area)){
  x <- area$Pos_X[s]
  y <- area$Pos_Y[s]
  if (x == 0 ){ x <- x + 1 }
  if (y == 0 ){ y <- y + 1 }
  mask <- ifelse(test2[ceiling(y),ceiling(x)] == 0,"adv",
                 ifelse(test2[ceiling(y),ceiling(x)] > 0 & test2[ceiling(y),ceiling(x)] <= 30000,"mid",
                        ifelse(test2[ceiling(y),ceiling(x)] > 30000 & test2[ceiling(y),ceiling(x)] <= 50000,"in","none")))
  qqq <- c(qqq,mask)
}
area$Area <- qqq
area$object <- spe2$object
spe2$Area <- area$Area[match(spe2$object, area$object)]

spe3 <- cbind(spe3,spe2)
}