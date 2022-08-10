#Author: Tong Liu, Zihuan Liu, Wenke Sun
library(dplyr)
library(Seurat)
library(patchwork)
library(wavelets)
library(Matrix)
library(scCATCH)
library(hdf5r)
library(readr)
library(UpSetR)
# Load the PBMC dataset

pbmc.data <- Read10X(data.dir = "D:/bioinfo/211110/pbmc0713",gene.column=2)

pbmc <- CreateSeuratObject(counts = pbmc.data, project = "pbmc3k", min.cells = 3, min.features = 200)
pbmc[["percent.mt"]] <- PercentageFeatureSet(pbmc, pattern = "^MT-")


VlnPlot(pbmc, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

plot1 <- FeatureScatter(pbmc, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(pbmc, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot1 + plot2

pbmc <- subset(pbmc, subset = nFeature_RNA > 200 & nFeature_RNA < 4200 & percent.mt < 5)
pbmc <- NormalizeData(pbmc, normalization.method = "LogNormalize", scale.factor = 10000)
pbmc_clean <- pbmc@assays$RNA@data
working_matrix=as.matrix(pbmc_clean)


## DWT
two_band <- function(X){
  n=dim(X)[1]
  nn=2*ceiling(n/2)
  m=dim(X)[2]
  VV=matrix(NA,nn,m)
  WW=matrix(NA,nn,m)
  XX=matrix(0,nn,m)
  XX[1:n,]=X
  
  for (i in 1:m){
    p<-XX[,i]
    B<-dwt(p,'d4',n.level=1) #d8 d12 ... for alternative 
    B1<-B
    B2<-B
    B1@W$W1<-matrix(0,nn/2,1)
    B2@V$V1<-matrix(0,nn/2,1)
    VV[,i]<-idwt(B1)
    WW[,i]<-idwt(B2)
  }
  VV=VV[1:n,]
  WW=WW[1:n,]
  result= list("low" = VV, "high" = WW)
  return (result)
}

three_band <- function(X){
  n=dim(X)[1]
  m=dim(X)[2]
  h0=c(0.33838609728386,0.53083618701374,0.72328627674361,0.23896417190576,0.04651408217589,-0.14593600755399)
  h1=c(-0.11737701613483,0.54433105395181,-0.01870574735313,-0.69911956479289,-0.13608276348796,0.42695403781698)
  h2=c(0.40363686892892,-0.62853936105471,0.46060475252131,-0.40363686892892,-0.07856742013185,0.24650202866523)
  #n=number of genes
  #m=number of cells
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


four_band <- function(X){
  n=dim(X)[1]
  m=dim(X)[2]
  h0<-(c(-0.067371764,0.094195111,0.40580489,0.567371764,0.567371764,0.40580489,0.094195111,-0.067371764))
  h1<-(c(-0.094195111,0.067371764, 0.567371764 ,0.40580489,-0.40580489,-0.567371764,-0.067371764,0.094195111))
  h2<-(c(-0.094195111,-0.067371764,0.567371764,-0.40580489,-0.40580489,0.567371764,-0.067371764,-0.094195111))
  h3<-(c(-0.067371764,-0.094195111,0.40580489,-0.567371764,0.567371764,-0.40580489,0.094195111,0.067371764))
  
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

four_band4 <- function(X){
  #4-band 4-regular
  n=dim(X)[1]
  m=dim(X)[2]
  h0<-(c(0.0857130200,0.1931394393,0.3491805097,0.5616494215,0.4955029828,0.4145647737,0.2190308939,-0.1145361261,
         -0.0952930728,-0.1306948909,-0.0827496793,0.0719795354,0.0140770701,0.0229906779,0.0145382757,-0.0190928308))
  h1<-(c(-0.1045086525,0.1183282069,-0.1011065044,-0.0115563891,0.6005913823,-0.2550401616,-0.4264277361,-0.0827398180,
         0.0722022649,0.2684936992,0.1691549718,-0.4437039320,0.0849964877,0.1388163056,0.0877812188,-0.1152813433))
  h2<-(c(0.2560950163,-0.2048089157,-0.2503433230,-0.2484277272,0.4477496752,0.0010274000,-0.0621881917,0.5562313118,
         -0.2245618041,-0.3300536827,-0.2088643503,0.2202951830,0.0207171125,0.0338351983,0.0213958651,-0.0280987676))
  h3<-(c(0.1839986022,-0.6622893130,0.6880085746,-0.1379502447,0.0446493766,-0.0823301969,-0.0923899104,-0.0233349758,
         0.0290655661,0.0702950474,0.0443561794,-0.0918374833,0.0128845052,0.0210429802,0.0133066389,-0.0174753464))
  
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
    ss=c(ss,ss[1:12])
    v1=rep(0,nn/4)
    w1=rep(0,nn/4)
    w2=rep(0,nn/4)
    w3=rep(0,nn/4)
    j=1
    for (i in 1:(nn/4)){
      v1[i]=sum(h0*ss[j:(j+15)])
      w1[i]=sum(h1*ss[j:(j+15)])
      w2[i]=sum(h2*ss[j:(j+15)])
      w3[i]=sum(h3*ss[j:(j+15)])
      j=j+4
    }
    v1=c(v1[(nn/4-2):(nn/4)],v1)
    w1=c(w1[(nn/4-2):(nn/4)],w1)
    w2=c(w2[(nn/4-2):(nn/4)],w2)
    w3=c(w3[(nn/4-2):(nn/4)],w3)
    for (i in 1:(nn/4)){
      vv1[4*i-3,k]=h0[13]*v1[i]+h0[9]*v1[i+1]+h0[5]*v1[i+2]+h0[1]*v1[i+3]
      vv1[4*i-2,k]=h0[14]*v1[i]+h0[10]*v1[i+1]+h0[6]*v1[i+2]+h0[2]*v1[i+3]
      vv1[4*i-1,k]=h0[15]*v1[i]+h0[11]*v1[i+1]+h0[7]*v1[i+2]+h0[3]*v1[i+3]
      vv1[4*i,k]=h0[16]*v1[i]+h0[12]*v1[i+1]+h0[8]*v1[i+2]+h0[4]*v1[i+3]
      ww1[4*i-3,k]=h1[13]*v1[i]+h1[9]*v1[i+1]+h1[5]*v1[i+2]+h1[1]*v1[i+3]
      ww1[4*i-2,k]=h1[14]*v1[i]+h1[10]*v1[i+1]+h1[6]*v1[i+2]+h1[2]*v1[i+3]
      ww1[4*i-1,k]=h1[15]*v1[i]+h1[11]*v1[i+1]+h1[7]*v1[i+2]+h1[3]*v1[i+3]
      ww1[4*i,k]=h1[16]*v1[i]+h1[12]*v1[i+1]+h1[8]*v1[i+2]+h1[4]*v1[i+3]
      ww2[4*i-3,k]=h2[13]*v1[i]+h2[9]*v1[i+1]+h2[5]*v1[i+2]+h2[1]*v1[i+3]
      ww2[4*i-2,k]=h2[14]*v1[i]+h2[10]*v1[i+1]+h2[6]*v1[i+2]+h2[2]*v1[i+3]
      ww2[4*i-1,k]=h2[15]*v1[i]+h2[11]*v1[i+1]+h2[7]*v1[i+2]+h2[3]*v1[i+3]
      ww2[4*i,k]=h2[16]*v1[i]+h2[12]*v1[i+1]+h2[8]*v1[i+2]+h2[4]*v1[i+3]
      ww3[4*i-3,k]=h3[13]*v1[i]+h3[9]*v1[i+1]+h3[5]*v1[i+2]+h3[1]*v1[i+3]
      ww3[4*i-2,k]=h3[14]*v1[i]+h3[10]*v1[i+1]+h3[6]*v1[i+2]+h3[2]*v1[i+3]
      ww3[4*i-1,k]=h3[15]*v1[i]+h3[11]*v1[i+1]+h3[7]*v1[i+2]+h3[3]*v1[i+3]
      ww3[4*i,k]=h3[16]*v1[i]+h3[12]*v1[i+1]+h3[8]*v1[i+2]+h3[4]*v1[i+3]
    }
  }
  vv1=vv1[1:n,]
  ww1=ww1[1:n,]
  ww2=ww2[1:n,]
  ww3=ww3[1:n,]
  
  result= list("low" = vv1, "high1" = ww1,'high2'=ww2,'high3'=ww3)
  return (result)
}


two_band_result<-two_band(working_matrix)
three_band_result<-three_band(working_matrix)
four_band_result<-four_band(working_matrix)


## PCA and UMAP
main_function<-function(data_,ob,X,t1,t2){
  pbmc_vv<-ob
  colnames(data_)<-colnames(X)
  rownames(data_)<-rownames(X)
  vv_sprase<-Matrix(data_, sparse = TRUE)
  pbmc_vv@assays$RNA@data<-vv_sprase
  pbmc_vv <- FindVariableFeatures(pbmc_vv, selection.method = "vst", nfeatures = 5000)
  top10_vv <- head(VariableFeatures(pbmc_vv), 10)
  
  # plot1 <- VariableFeaturePlot(pbmc_vv)
  # plot2 <- LabelPoints(plot = plot1, points = top10_vv, repel = TRUE)
  # plot1 + plot2
  
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
  
  #DimPlot(pbmc_vv, reduction = "umap")
  
  pbmc_vv.markers <- FindAllMarkers(pbmc_vv, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
  
  result= list("data" = pbmc_vv, "markers" = pbmc_vv.markers)
  return(result)
}


##########cell identification
cell_identity<-function(data,tissue_name,match_CellMatch_,cancer){
  clu_markers <- findmarkergenes(object = data,species = 'Human',cluster = 'All',match_CellMatch = match_CellMatch_,
                                 cancer = cancer,
                                 tissue = tissue_name,
                                 cell_min_pct = 0.1,
                                 logfc = 0.1,
                                 pvalue = 0.01)
  
  clu_ann <- scCATCH(object = clu_markers$clu_markers,
                     species = 'Human',
                     cancer = cancer,
                     tissue = tissue_name)
  
  result= list("clu_markers " = clu_markers , "indenity" =clu_ann )
  return (result)
}
#############


plot_umap<-function(data,list_name){
  new.cluster.ids <- list_name
  names(new.cluster.ids) <- levels(data)
  data<- RenameIdents(data, new.cluster.ids)
  
  return (DimPlot(data, reduction = "umap", label = TRUE, pt.size = 0.8) )
}

cell_tppe_name<-function(d1,d2){
  
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


###*************for covid data settting :****************
tissue_name<-'Lung'
match_CellMatch_<-TRUE
cancer_name<-NULL


#***********for cancer data setting:*************
tissu_name<-'Peripheral blood'
mtch_CellMatch_name<-FALSE
cancer_name<-NULL

two_band_low<-main_function(two_band_result$low,pbmc,working_matrix,'2band','low')
two_band_low_cell_identity<-cell_identity(two_band_low$data,tissu_name,mtch_CellMatch_name, cancer_name)
two_band_low_list_name<-cell_tppe_name(two_band_low_cell_identity,two_band_low)
plot_umap(two_band_low$data,two_band_low_list_name)

two_band_high<-main_function(two_band_result$high,pbmc,working_matrix,'2band','high')
two_band_high_cell_identity<-cell_identity(two_band_high$data,tissu_name,mtch_CellMatch_name, cancer_name)
two_band_high_list_name<-cell_tppe_name(two_band_high_cell_identity,two_band_high)
plot_umap(two_band_high$data,two_band_high_list_name)


three_band_low<-main_function(three_band_result$low,pbmc,working_matrix,'3band','low')
three_band_low_cell_identity<-cell_identity(three_band_low$data,tissu_name,mtch_CellMatch_name, cancer_name)
three_band_low_list_name<-cell_tppe_name(three_band_low_cell_identity,three_band_low)
plot_umap(three_band_low$data,three_band_low_list_name)

three_band_h1<-main_function(three_band_result$high1,pbmc,working_matrix,'3band','high1')
three_band_h1_cell_identity<-cell_identity(three_band_h1$data,tissu_name,mtch_CellMatch_name, cancer_name)
three_band_h1_list_name<-cell_tppe_name(three_band_h1_cell_identity,three_band_h1)
plot_umap(three_band_h1$data,three_band_h1_list_name)

three_band_h2<-main_function(three_band_result$high2,pbmc,working_matrix,'3band','high2')
three_band_h2_cell_identity<-cell_identity(three_band_h2$data,tissu_name,mtch_CellMatch_name, cancer_name)
three_band_h2_list_name<-cell_tppe_name(three_band_h2_cell_identity,three_band_h2)
plot_umap(three_band_h2$data,three_band_h2_list_name)


four_band_low<-main_function(four_band_result$low,pbmc,working_matrix,'4band','low')
four_band_low_cell_identity<-cell_identity(four_band_low$data,tissu_name,mtch_CellMatch_name, cancer_name)
four_band_low_list_name<-cell_tppe_name(four_band_low_cell_identity,four_band_low)
plot_umap(four_band_low$data,four_band_low_list_name)


four_band_h1<-main_function(four_band_result$high1,pbmc,working_matrix,'4band','high1')
four_band_h1_cell_identity<-cell_identity(four_band_h1$data,tissu_name,mtch_CellMatch_name, cancer_name)
four_band_h1_list_name<-cell_tppe_name(four_band_h1_cell_identity,four_band_h1)
plot_umap(four_band_h1$data,four_band_h1_list_name)

four_band_h2<-main_function(four_band_result$high2,pbmc,working_matrix,'4band','high2')
four_band_h2_cell_identity<-cell_identity(four_band_h2$data,tissu_name,mtch_CellMatch_name, cancer_name)
four_band_h2_list_name<-cell_tppe_name(four_band_h2_cell_identity,four_band_h2)
plot_umap(four_band_h2$data,four_band_h2_list_name)

four_band_h3<-main_function(four_band_result$high3,pbmc,working_matrix,'4band','high3')
four_band_h3_cell_identity<-cell_identity(four_band_h3$data,tissu_name,mtch_CellMatch_name, cancer_name)
four_band_h3_list_name<-cell_tppe_name(four_band_h3_cell_identity,four_band_h3)
plot_umap(four_band_h3$data,four_band_h3_list_name)

orginal<-main_function(working_matrix,pbmc,working_matrix,"orginial","na")
orginal_cell_identity<-cell_identity(orginal$data,tissu_name,mtch_CellMatch_name, cancer_name)
orginal_cell_list_name<-cell_tppe_name(orginal_cell_identity,orginal)

####calculating sparsity
spar=function(A){
  n=dim(A)
  B=as.integer(as.logical(A))
  return(sum(sum(B))/n[1]/n[2])
}

## manually check celltype

data = read.csv(file = "D:/bioinfo/211110/Book1.csv",header = T)
#celltype book refers to Cellmarker 2.0
x = "EEF1A1, IL7R, NOP53, PTGER4, RPL10, RPL19, RPL34, RPL39, RPL9, RPS12, RPS15A, RPS19, RPS21, RPS6, TOMM7, TPT1" 
cellmarksers = strsplit(x,', ')
a = c()
for(i in 1:length(cellmarksers[[1]])){
  a = c(a,which(data$Marker ==cellmarksers[[1]][i]))
}
table(data$Cell.name[a]) #CD4+ T cell for above example



####### check cells for clusters
colors = c(rgb(255, 99, 71, 40, maxColorValue=255),rgb(210, 105, 30, 80, maxColorValue=255),rgb(107, 142, 35, 60, maxColorValue=255),
           rgb(64, 224, 208, 50, maxColorValue=255),rgb(30, 144, 255, 100, maxColorValue=255),rgb(218, 112, 214, 150, maxColorValue=255),
           rgb(138, 43, 226, 200, maxColorValue=255),rgb(227, 23, 13, 100, maxColorValue=255),rgb(50, 205, 50, 100, maxColorValue=255),
           rgb(255, 0, 0, 100, maxColorValue=255))

check_cluster_cells <- function(data,data_outputs){
  
  n = length(data_outputs)
  name = data_outputs
  name[is.na(name)] ="Unknown"
  names = c()
  cell_freq = as.data.frame(table(data$data@meta.data$seurat_clusters))
  for(i in 1:n){
    l = cell_freq[i,2]
    names = append(names,rep(name[i],l))
  }
  return(table(names))
  #return(data$data@assays$RNA@counts@Dim[2]-length(which(names!='Unknown')))
}

check_cluster_cells(orginal,orginal_cell_list_name)
check_cluster_cells(four_band_h3,four_band_h3_list_name)

colors_input <-function(gene,names,cell_identity){
  r = c(255,210,127,64,30,218,138,227,50,255)
  g = c(99,105,255,224,144,112,43,23,205,0)
  bl = c(71,30,0,208,255,214,226,13,50,0)
  m = length(unique(names))
  c = as.data.frame(matrix(0,m,1))
  c[,1] = unique(names)
  a = check_cluster_cells(gene,names)
  a = as.data.frame(a)
  b = c %>% left_join(a,by= c("V1" = "names" ))
  bb=b[,2]
  bb[is.na(bb)]=1
  b[,2]=bb
  b[,2] = b[,2]/max(b[,2])
  bb=b[,2]
  bb[which(bb<0.1)]=0.1
  b[,2]=bb*255
  colors = c(rgb(r[1],g[1],bl[1],b[1,2],maxColorValue=255))
  for(i in 2:m){
    colors = c(colors,rgb(r[i],g[i],bl[i],b[i,2],maxColorValue=255))
  }
  return (colors)
}

plot_umap<-function(data,list_name){
  new.cluster.ids <- list_name
  names(new.cluster.ids) <- levels(data)
  data<- RenameIdents(data, new.cluster.ids)
  
  return (DimPlot(data, reduction = "umap", label = F, pt.size = 0.5,cols =colours_input))
}

colours_input = colors_input(two_band_low,two_band_low_list_name,two_band_low_cell_identity)
png("45232l.png",width = 6,height = 4,units = "in",res = 800,pointsize = 9)
par(mar = c(2,2,1,1),xaxs = "i",yaxs = "i",cex.axis = 3,cex.lab = 3,cex.main=3,cex.sub=3)
plot_umap(two_band_low$data,two_band_low_list_name)
dev.off()

###calculating pvalues and plotting heatmaps among clusters

numcells=dim(working_matrix)[2]
ctlist=matrix(nrow=numcells,ncol=10)
ctlist[,1]=orginal$data$seurat_clusters
ctlist[,2]=two_band_low$data$seurat_clusters
ctlist[,3]=two_band_high$data$seurat_clusters
ctlist[,4]=three_band_low$data$seurat_clusters
ctlist[,5]=three_band_h1$data$seurat_clusters
ctlist[,6]=three_band_h2$data$seurat_clusters
ctlist[,7]=four_band_low$data$seurat_clusters
ctlist[,8]=four_band_h1$data$seurat_clusters
ctlist[,9]=four_band_h2$data$seurat_clusters
ctlist[,10]=four_band_h3$data$seurat_clusters
inters=function(a,b){
  n1=max(ctlist[,a])
  n2=max(ctlist[,b])
  tab=matrix(0,n1,n2)
  for (i in 1:n1){
    for (j in 1:n2){
      tab[i,j]=length(which((ctlist[,a]==i)&(ctlist[,b]==j)))
    }
  }
  return(tab)
}
pvalu=function(X){
  n1=nrow(X)
  n2=ncol(X)
  N=sum(X)
  p=matrix(0,n1,n2)
  for (i in 1:n1){
    for (j in 1:n2){

      pro=sum(X[i,])*sum(X[,j])/N^2

      p[i,j]=phyper(X[i,j]-0.5,sum(X[i,]),sum(X[-i,]),sum(X[,j]))
      #p[i,j]=pbinom(X[i,j]-0.5,N,pro)
    }
  }
  return(p)
}
titl=c("canonical","2-Band low","2-Band high","3-Band low","3-Band high1","3-Band high2",
       "4-Band low","4-Band high1","4-Band high2","4-Band high3")



htmp=function(a,b){
  da=inters(a,b)
  A=pvalu(da)
  clu=length(which((A>0.01)&(da>=3)))
  B=format(round(A,digits=2),nsmall=3,scientific=F)
  A[which(da<5)]=1e-16
  for (i in 1:dim(da)[1]){
    for (j in 1:dim(da)[2])
      da[i,j]=paste(da[i,j],'(',B[i,j],')',sep="")
  }
  
  mn=paste(titl[a],'vs',titl[b],",",clu,'clusters')
  heatmap.2(A,Rowv=FALSE,Colv=FALSE,cellnote=da,notecol=1,trace="none",col=heat.colors(100,rev=TRUE),par(cex.main = 10,cex.lab = 10,cex.axis =10,cex.sub =10),
            main=mn,density.info="none",notecex = 2,margins = c(0.5,0.5),key = F,keysize = 0.5)#,keysize = 10,key.par=list(mgp=c(1.5, 0.5, 0),mar=c(1,1,1,0))) #cex.lab = 15,cex.axis =15,cex.main = 15,cex.sub =15)
  
}
col=heat.colors(100,rev=TRUE)


library(gplots)

cc1=c(2,4,4,5,7,7,7,8,8,9)
cc2=c(3,5,6,6,8,9,10,9,10,10)
for (i in 1:10){
  flmn=paste(titl[cc1[i]],titl[cc2[i]],".png")
  png(filename=flmn,width=1200,height=1200)
  htmp(cc1[i],cc2[i])
  dev.off()
}

##calculating numbers of significant clusters for M-band with M>=3

n1=max(ctlist[,4])
n2=max(ctlist[,5])
n3=max(ctlist[,6])
cell3=array(0,dim=c(n1,n2,n3))
for (i in 1:n1){
  for (j in 1:n2){
    for (k in 1:n3){
      cell3[i,j,k]=length(which((ctlist[,4]==i)&(ctlist[,5]==j)&(ctlist[,6]==k)))
    }
  }
}
pval3=array(0,dim=c(n1,n2,n3))
N3=sum(cell3)
sigclu3=matrix(0,N3,3)
nclu3=0
for (i in 1:n1){
  for (j in 1:n2){
    for (k in 1:n3){
      pro=sum(cell3[i,,])*sum(cell3[,j,])*sum(cell3[,,k])/N3^3
      pval3[i,j,k]=pbinom(cell3[i,j,k]-0.5,N3,pro)
      if ((pval3[i,j,k]>0.01)&(cell3[i,j,k])>=5){
        nclu3=nclu3+1
        sigclu3[nclu3,1]=i
        sigclu3[nclu3,2]=j
        sigclu3[nclu3,3]=k
      }
    }
  }
}
sigclu3=sigclu3[1:nclu3,]

n1=max(ctlist[,7])
n2=max(ctlist[,8])
n3=max(ctlist[,9])
n4=max(ctlist[,10])
cell4=array(0,dim=c(n1,n2,n3,n4))
for (i in 1:n1){
  for (j in 1:n2){
    for (k in 1:n3){
      for (l in 1:n4){
        cell4[i,j,k,l]=length(which((ctlist[,7]==i)&(ctlist[,8]==j)&(ctlist[,9]==k)&(ctlist[,10]==l)))
      }
    }
  }
}
pval4=array(0,dim=c(n1,n2,n3,n4))
N4=sum(cell4)
sigclu4=matrix(0,N4,4)
nclu4=0
for (i in 1:n1){
  for (j in 1:n2){
    for (k in 1:n3){
      for (l in 1:n4){
        pro=sum(cell4[i,,,])*sum(cell4[,j,,])*sum(cell4[,,k,])*sum(cell4[,,,l])/N4^4
        pval4[i,j,k,l]=pbinom(cell4[i,j,k,l]-0.5,N4,pro)
        if ((pval4[i,j,k,l]>0.01)&(cell4[i,j,k,l])>=5){
          nclu4=nclu4+1
          sigclu4[nclu4,1]=i
          sigclu4[nclu4,2]=j
          sigclu4[nclu4,3]=k
          sigclu4[nclu4,4]=l
        }
      }
    }
  }
}
sigclu4=sigclu4[1:nclu4,]

#intersection for allgenes
A0=orginal$markers %>% filter(avg_log2FC>=log(1.5,2)) %>% filter(p_val_adj<0.01)
A40=four_band_low$markers%>% filter(avg_log2FC>=log(1.5,2)) %>% filter(p_val_adj<0.01)
A41=four_band_h1$markers%>% filter(avg_log2FC>=log(1.5,2)) %>% filter(p_val_adj<0.01)
A42=four_band_h2$markers%>% filter(avg_log2FC>=log(1.5,2)) %>% filter(p_val_adj<0.01)
A43=four_band_h3$markers%>% filter(avg_log2FC>=log(1.5,2)) %>% filter(p_val_adj<0.01)

N=nrow(A0)+nrow(A40)+nrow(A41)+nrow(A42)+nrow(A43)
A1=matrix(NA,N,5)
grp=c(rep(1,nrow(A0)),rep(2,nrow(A40)),rep(3,nrow(A41)),rep(4,nrow(A42)),rep(5,nrow(A43)))
A1=rbind(A0,A40,A41,A42,A43)
A2=cbind(grp,A1)
A2=as.data.frame(A2)
A3=A2[!duplicated(A2[,8]),8]
A3=as.data.frame(A3)
B=matrix(0,nrow(A3),5) #2band change 5 to 3;3band change 5 to 4
for (i in 1: nrow(A3)){
  ind=which(A2[,8]==A3[i,1])
  colu=as.integer(A2[ind,1])
  B[i,colu]=1
}
C=cbind(A3,B)
colnames(C)=c("gene","canonical","4-band low","4-band h1","4-band h2","4-band h3")
#colnames(C)=c("gene","original","Daub4 low","Daub4 high")
png("3921j.png",width=15,height=10,units='in',res=600,pointsize=3)
upset(mak,order.by = c("freq"), decreasing = c(TRUE),sets.bar.color = 5, main.bar.color = 5, matrix.color = 4,text.scale=3)
dev.off()


#intersection for markers
fenli=function(x){
  obje=NULL
  for (i in 1:length(x)){
    mk=strsplit(x[i],', ')
    if(mk!='NA')
      obje=c(obje,mk[[1]])
  }
  return(obje)
}

origmk=fenli(orginal_cell_identity$indenity$celltype_related_marker)
mk4l=fenli(four_band_low_cell_identity$indenity$celltype_related_marker)
mk4h1=fenli(four_band_h1_cell_identity$indenity$celltype_related_marker)
mk4h2=fenli(four_band_h2_cell_identity$indenity$celltype_related_marker)
mk4h3=fenli(four_band_h3_cell_identity$indenity$celltype_related_marker)

allmk=c(origmk,mk4l,mk4h1,mk4h2,mk4h3)
allmk=allmk[!duplicated(allmk)]

mkdis=matrix(0,length(allmk),5) #2band change 5 to 3;3band change 5 to 4
for (i in 1:length(allmk)){
  mkdis[i,1]=length(which(origmk==allmk[i]))
  mkdis[i,2]=length(which(mk4l==allmk[i]))
  mkdis[i,3]=length(which(mk4h1==allmk[i]))
  mkdis[i,4]=length(which(mk4h2==allmk[i]))
  mkdis[i,5]=length(which(mk4h3==allmk[i]))
  
}
rownames(mkdis)=allmk
mkd=as.data.frame(mkdis)
colnames(mkd)=c("canonical","4band low","4band high1","4band high2","4band high3")

png("3921j.png",width=15,height=10,units='in',res=600,pointsize=3)
upset(mak,order.by = c("freq"), decreasing = c(TRUE),sets.bar.color = 5, main.bar.color = 5, matrix.color = 4,text.scale=3)
dev.off()