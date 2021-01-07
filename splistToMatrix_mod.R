splistToMatrix = function(d,sp,r)
  {
  out = matrix(0,nrow = ncell(r),ncol = length(sp))
  d2 = d[which(d$name %in% sp),]
  
  rowind = cellFromXY(r,cbind(d2$longitude,d2$latitude))
  colind = match(d2$name,sp)
  out[cbind(rowind,colind)] = 1

  out[,which(apply(out,2,sum)==0)] = NA
  rownames(out) = paste("Cell",1:ncell(r),sep ="_")
  colnames(out) = sp
  return(out)
  }
