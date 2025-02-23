#packages
library(dplyr)
library(Seurat)
library(patchwork)
library(wavelets)
library(Matrix)
library(scCATCH)
library(hdf5r)
library(readr)
library(SAVER)
library(parallel)
library(tidyverse)
#wavelet functions
# Define the wavelet transformation function
#Arguments 
#wavelet_type "2-band","3-band","4-band"
#daubechies_order "d4","d8","d12" or "d16"
#regular True if regular = 4, False regular = 2
wavelet_transform <- function(input_matrix, wavelet_type , daubechies_order = "d4", regular = FALSE) {
  
  # 2-band wavelet transform
  if (wavelet_type == "2-band") {
    # Check if the Daubechies order is valid
    if (!daubechies_order %in% c("d4","d8","d12","d16")) {
      stop("For 2-band wavelet transform, Daubechies order must be 4, 8, 12, or 16.")
    }
    
    # Perform 2-band wavelet transform
    
    result <- two_band(input_matrix,daubechies_order)
  }
  
  # 3-band wavelet transform
  else if (wavelet_type == "3-band") {
    # Perform 3-band wavelet transform
    result <- three_band(input_matrix)
  }
  
  # 4-band wavelet transform
  else if (wavelet_type == "4-band") {
    if (regular) {
      # Use 4-regular wavelet
      result <- four_band4(input_matrix)
    } else {
      # Use default 4-band wavelet
      result <- four_band(input_matrix)
    }
  }
  
  # If the wavelet type is invalid
  else {
    stop("Unsupported wavelet type. Please choose '2-band', '3-band', or '4-band'.")
  }
  
  # Return the wavelet transform result
  return(result)
}
#2-band wavelet
two_band <- function(X,daubechies_order){
  n=dim(X)[1]
  m=dim(X)[2]
  VV=matrix(NA,n,m)
  WW=matrix(NA,n,m)
  #d8 d12
  for (i in 1:m){
    p<-X[,i]
    B<-dwt(p,daubechies_order,n.level=1)
    
    B1<-B
    B2<-B
    B1@W$W1<-matrix(0,n/2,1)
    B2@V$V1<-matrix(0,n/2,1)
    VV[,i]<-idwt(B1)
    WW[,i]<-idwt(B2)
  }
  
  result= list("low" = VV, "high" = WW)
  return (result)
}
#3-band wavelet
three_band <- function(X){
  n=dim(X)[1]
  m=dim(X)[2]
  h0=c(0.33838609728386,0.53083618701374,0.72328627674361,0.23896417190576,0.04651408217589,-0.14593600755399)
  h1=c(-0.11737701613483,0.54433105395181,-0.01870574735313,-0.69911956479289,-0.13608276348796,0.42695403781698)
  h2=c(0.40363686892892,-0.62853936105471,0.46060475252131,-0.40363686892892,-0.07856742013185,0.24650202866523)
  #n=17125  ##number of genes
  #m=2359  ##number of cells
  nn=3*ceiling(n/3)
  
  vv1=matrix(0,nn,m)
  ww1=matrix(0,nn,m)
  ww2=matrix(0,nn,m)
  for (k in 1:m){
    ss=X[,k]
    if(n%%3){ss=c(ss,rep(0,3-n%%3))}
    ss=c(ss,ss[1:3])
    v1=rep(0,nn/3)
    w1=rep(0,nn/3)
    w2=rep(0,nn/3)
    j=1
    for (i in 1:(nn/3)){
      v1[i]=sum(h0*ss[j:(j+5)])
      w1[i]=sum(h1*ss[j:(j+5)])
      w2[i]=sum(h2*ss[j:(j+5)])
      j=j+3
    }
    v1=c(v1[nn/3],v1)
    w1=c(w1[nn/3],w1)
    w2=c(w2[nn/3],w2)
    
    for (i in 1:(nn/3)){
      vv1[3*i-2,k]=h0[4]*v1[i]+h0[1]*v1[i+1]
      vv1[3*i-1,k]=h0[5]*v1[i]+h0[2]*v1[i+1]
      vv1[3*i,k]=h0[6]*v1[i]+h0[3]*v1[i+1]
      ww1[3*i-2,k]=h1[4]*w1[i]+h1[1]*w1[i+1]
      ww1[3*i-1,k]=h1[5]*w1[i]+h1[2]*w1[i+1]
      ww1[3*i,k]=h1[6]*w1[i]+h1[3]*w1[i+1]
      ww2[3*i-2,k]=h2[4]*w2[i]+h2[1]*w2[i+1]
      ww2[3*i-1,k]=h2[5]*w2[i]+h2[2]*w2[i+1]
      ww2[3*i,k]=h2[6]*w2[i]+h2[3]*w2[i+1]
    }
  }
  vv1=vv1[1:n,]
  ww1=ww1[1:n,]
  ww2=ww2[1:n,]
  result= list("low" = vv1, "high1" = ww1,'high2'=ww2)
  return (result)
}

#4-band wavelet
four_band <- function(X){
  n=dim(X)[1]
  m=dim(X)[2]
  h0<-(c(-0.067371764,0.094195111,0.40580489,0.567371764,0.567371764,0.40580489,0.094195111,-0.067371764))
  h1<-(c(-0.094195111,0.067371764, 0.567371764 ,0.40580489,-0.40580489,-0.567371764,-0.067371764,0.094195111))
  h2<-(c(-0.094195111,-0.067371764,0.567371764,-0.40580489,-0.40580489,0.567371764,-0.067371764,-0.094195111))
  h3<-(c(-0.067371764,-0.094195111,0.40580489,-0.567371764,0.567371764,-0.40580489,0.094195111,0.067371764))
  
  #n=17125  ##number of genes
  #m=2359  ##number of cells
  nn=4*ceiling(n/4)
  
  vv1=matrix(0,nn,m)
  ww1=matrix(0,nn,m)
  ww2=matrix(0,nn,m)
  ww3=matrix(0,nn,m)
  for (k in 1:m){
    ss=X[,k]
    if(n%%4){ss=c(ss,rep(0,4-n%%4))}
    ss=c(ss,ss[1:4])
    v1=rep(0,nn/4)
    w1=rep(0,nn/4)
    w2=rep(0,nn/4)
    w3=rep(0,nn/4)
    j=1
    for (i in 1:(nn/4)){
      v1[i]=sum(h0*ss[j:(j+7)])
      w1[i]=sum(h1*ss[j:(j+7)])
      w2[i]=sum(h2*ss[j:(j+7)])
      w3[i]=sum(h3*ss[j:(j+7)])
      j=j+4
    }
    v1=c(v1[nn/4],v1)
    w1=c(w1[nn/4],w1)
    w2=c(w2[nn/4],w2)
    w3=c(w3[nn/4],w3)
    for (i in 1:(nn/4)){
      vv1[4*i-3,k]=h0[5]*v1[i]+h0[1]*v1[i+1]
      vv1[4*i-2,k]=h0[6]*v1[i]+h0[2]*v1[i+1]
      vv1[4*i-1,k]=h0[7]*v1[i]+h0[3]*v1[i+1]
      vv1[4*i,k]=h0[8]*v1[i]+h0[4]*v1[i+1]
      ww1[4*i-3,k]=h1[5]*w1[i]+h1[1]*w1[i+1]
      ww1[4*i-2,k]=h1[6]*w1[i]+h1[2]*w1[i+1]
      ww1[4*i-1,k]=h1[7]*w1[i]+h1[3]*w1[i+1]
      ww1[4*i,k]=h1[8]*w1[i]+h1[4]*w1[i+1]
      ww2[4*i-3,k]=h2[5]*w2[i]+h2[1]*w2[i+1]
      ww2[4*i-2,k]=h2[6]*w2[i]+h2[2]*w2[i+1]
      ww2[4*i-1,k]=h2[7]*w2[i]+h2[3]*w2[i+1]
      ww2[4*i,k]=h2[8]*w2[i]+h2[4]*w2[i+1]
      ww3[4*i-3,k]=h3[5]*w3[i]+h3[1]*w3[i+1]
      ww3[4*i-2,k]=h3[6]*w3[i]+h3[2]*w3[i+1]
      ww3[4*i-1,k]=h3[7]*w3[i]+h3[3]*w3[i+1]
      ww3[4*i,k]=h3[8]*w3[i]+h3[4]*w3[i+1]
    }
  }
  vv1=vv1[1:n,]
  ww1=ww1[1:n,]
  ww2=ww2[1:n,]
  ww3=ww3[1:n,]
  
  result= list("low" = vv1, "high1" = ww1,'high2'=ww2,'high3'=ww3)
  return (result)
}



main_function<-function(data_,ob,X,t1,t2){
  pbmc_vv<-ob
  colnames(data_)<-colnames(X)
  rownames(data_)<-rownames(X)
  vv_sprase<-Matrix(data_, sparse = TRUE)
  pbmc_vv@assays$RNA@data<-vv_sprase
  pbmc_vv <- FindVariableFeatures(pbmc_vv, selection.method = "vst", nfeatures = 5000)
  top10_vv <- head(VariableFeatures(pbmc_vv), 10)
  
  plot1 <- VariableFeaturePlot(pbmc_vv)
  plot2 <- LabelPoints(plot = plot1, points = top10_vv, repel = TRUE)
  plot1 + plot2
  
  all.genes <- rownames(pbmc_vv)
  pbmc_vv <- ScaleData(pbmc_vv, features = all.genes)
  
  pbmc_vv <- RunPCA(pbmc_vv, features = VariableFeatures(object = pbmc_vv))
  
  VizDimLoadings(pbmc_vv, dims = 1:2, reduction = "pca")
  DimPlot(pbmc_vv, reduction = "pca")
  pbmc_vv <- JackStraw(pbmc_vv, num.replicate = 100)
  pbmc_vv <- ScoreJackStraw(pbmc_vv, dims = 1:20)
  
  JackStrawPlot(pbmc_vv, dims = 1:15)
  ElbowPlot(pbmc_vv)
  
  pbmc_vv <- FindNeighbors(pbmc_vv, dims = 1:10)
  pbmc_vv <- FindClusters(pbmc_vv, resolution = 0.5)
  head(Idents(pbmc_vv), 10)
  pbmc_vv <- RunUMAP(pbmc_vv, dims = 1:10)
  
  
  pbmc_vv.markers <- FindAllMarkers(pbmc_vv, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
  

  
  result= list("data" = pbmc_vv, "markers" = pbmc_vv.markers)
  return(result)
}


##########fix for different disease:FLASE, Logfc=0.2
cell_identity<-function(data,tissue_name,match_CellMatch_,cancer){
  clu_markers <- findmarkergenes(object = data,species = 'Human',cluster = 'All',match_CellMatch = match_CellMatch_,
                                 cancer = cancer,
                                 tissue = tissue_name,
                                 cell_min_pct = 0.1,#0.25,
                                 logfc = 0.1,#0.25,
                                 pvalue = 0.01)#0.05
  
  clu_ann <- scCATCH(object = clu_markers$clu_markers,
                     species = 'Human',
                     cancer = cancer,
                     tissue = tissue_name)
  
  result= list("clu_markers " = clu_markers , "indenity" =clu_ann )
  return (result)
}
########
#############


plot_umap<-function(data,list_name){
  new.cluster.ids <- list_name
  names(new.cluster.ids) <- levels(data)
  data<- RenameIdents(data, new.cluster.ids)
  
  return (DimPlot(data, reduction = "umap", label = T, pt.size = 0.8,label.size = 2)) 
}



cell_type_name<-function(d1,d2){
  
  
  name=d1$indenity$cell_type
  d3=d2$data@meta.data$seurat_clusters
  len_=length(name)
  for(i in 1:len_){
    if (is.na(name[i])){
      name[i]='Unknown'
    }
  }
  
  a=as.numeric(unique(d3))-1
  new_name=1:max(a)
  b=as.numeric(d1$indenity$cluster)
  c=setdiff(a,b)+1
  
  if (length(c)!=0){
    for (j in c){
      new_name[j]='Unknown'
    }
    
    k=1
    for (j in (b+1)){
      new_name[j]=name[k]
      k=k+1
    }
    return (new_name)
  }
  else{
    return (name)
  }
}

#Part 4 High-resolution image

colors_input <-function(gene,names,cell_identity){
  r = c(255,210,127,64,30,218,138,227,50,255,255,255,160)
  g = c(99,105,255,224,144,112,43,23,205,0,153,255,32)
  bl = c(71,30,0,208,255,214,226,13,50,0,18,0,240)
  m = length(unique(names))
  c = as.data.frame(matrix(0,m,1))
  c[,1] = unique(names)
  a = check_cluster_cells(gene,cell_identity)
  a = as.data.frame(a)
  b = c %>% left_join(a,by= c("V1" = "names" ))
  bb = b[,2]
  bb[is.na(bb)] =1
  b[,2] = bb
  b[,2] = b[,2]/max(b[,2])
  bb = b[,2]
  bb[which(bb<0.1)] = 0.1
  b[,2]= bb*255
  colors = c(rgb(r[1],g[1],bl[1],b[1,2],maxColorValue=255))
  for(i in 2:m){
    colors = c(colors,rgb(r[i],g[i],bl[i],b[i,2],maxColorValue=255))
  }
  return (colors)
}

plot_umap1<-function(data,list_name){
  new.cluster.ids <- list_name
  names(new.cluster.ids) <- levels(data)
  data<- RenameIdents(data, new.cluster.ids)
  
  return (DimPlot(data, reduction = "umap", label = F, pt.size = 0.5,cols =colours_input))
}