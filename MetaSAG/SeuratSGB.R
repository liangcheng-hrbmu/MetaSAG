
###############使用R包Seurat对所有KSGB的细胞进行聚类###############
###################################################################

library(Seurat)
library(dplyr)
library(getopt)

spec <- matrix(
c("CellUnirefCount","i", 1, "character",TRUE,"CellUnirefCount","This is All SGB Cell's Uniref Count Matrix File.",
  "OutPut", "o", 2, "character",TRUE,"Output","This is th OutPutFolder",
  "nfeatures","nfeatures",3,"integer",FALSE,20000,"This is argument ,nfeatures, for function FindVariableFeatures in Seurat",
  "npcs","npcs",4,"integer",FALSE,30,"This is argument ,npcs, for function RunPCA in Seurat",
  "FindNeighborsDim","FindNeighborsDim",5,"character",FALSE,"1:10","This is argument ,dims, for function FindNeighbors in Seurat",
  "FindClusterRes","FindClusterRes",6,"double",FALSE,0.5,"This is argumet ,resolution, for function FindCluster in Seurat",
  "UmapPDFWidth","UmapPDFWidth",7,"integer",FALSE,6,"This is argument ,UmapPDFWidth, for specify the width of PDF KnownSGBUmap",
  "UmapPDFHeight","UmapPDFHeight",8,"integer",FALSE,4,"This is argument ,UmapPDFHeight, for specify the heigth of PDF KnownSGBUmap",
  "DimPlotMethod","DimPlotMethod",9,"character",FALSE,"umap","This is argument ,method, for function DimPlot in Seurat",
  "FindAllMarkersPCT","FindAllMarkersPCT",10,"double",FALSE,0.25,"This is argument ,min.pct, for function FindAllMarkers in Seurat",
  "FindAllMarkerslogFC","FindAllMarkerslogFC",11,"double",FALSE,0.25,"This is argument ,logfc.threshold, for function FindAllMarkers in Seurat",
  "topN","topN",12,"integer",FALSE,10,"This is argument ,topN, for function top_n in Seurat",
  "HeatmapPDFWidth","HeatmapPDFWidth",13,"integer",FALSE,12,"This is argument ,HeatmapPDFWidth, for specify the width of PDF Heatmap",
  "HeatmapPDFHeight","HeatmapPDFHeight",14,"integer",FALSE,8,"This is argument ,HeatmapPDFHeight, for specify the Height of PDF Heatmap"),byrow=TRUE,ncol=7)

opt = getopt(spec)
CountMatrixFile = opt$CellUnirefCount
OutPutDirt = opt$OutPut
nfeatures = opt$nfeatures
npcs = opt$npcs
FindNeighborsDim = opt$FindNeighborsDim
FindNeighborsDim1 = as.integer(strsplit(FindNeighborsDim)[[1]][1])
FindNeighborsDim2 = as.integer(strsplit(FindNeighborsDim)[[1]][2])
FindClusterRes = opt$FindClusterRes
UmapPDFWidth = opt$UmapPDFWidth
UmapPDFHeight = opt$UmapPDFHeight
DimPlotMethod = opt$DimPlotMethod
FindAllMarkersPCT = opt$FindAllMarkersPCT
FindAllMarkerslogFC =opt$FindAllMarkerslogFC
topN = opt$topN
HeatmapPDFWidth = opt$HeatmapPDFWidth
HeatmapPDFHeight = opt$HeatmapPDFHeight

count <- read.delim(CountMatrixFile, row.names=1)
f=function(x){
  temp=length(which(x>0))
  if(temp<=5){
    return(FALSE)
  }else{
    return(TRUE)
  }

}

xx=apply(count,1,f)
count=count[xx,]
print(paste0('细胞数目为：',ncol(count)))
scDNA <- CreateSeuratObject(counts = count)
scDNA = NormalizeData(scDNA)
scDNA <- FindVariableFeatures(scDNA, selection.method = 'vst', nfeatures = nfeatures)
scDNA <- ScaleData(scDNA)
scDNA=RunPCA(scDNA,npcs=npcs,features=VariableFeatures(scDNA))
scDNA <- FindNeighbors(scDNA, dims = FindNeighborsDim1:FindNeighborsDim2)
scDNA <- FindClusters(scDNA, resolution = FindClusterRes)
table(scDNA$seurat_clusters)
scDNA <- RunUMAP(scDNA, dims = FindNeighborsDim1:FindNeighborsDim2)
pdf(paste0(OutPutDirt,"/KnownSGBUmap.pdf"), width = UmapPDFWidth, height = UmapPDFHeight)
DimPlot(scDNA, reduction = DimPlotMethod,label = TRUE)
dev.off()
scDNA.markers <- FindAllMarkers(scDNA, only.pos = TRUE, min.pct = FindAllMarkersPCT, logfc.threshold = FindAllMarkerslogFC)
top10markers<-scDNA.markers %>% group_by(cluster) %>% top_n(n = topN, wt = avg_log2FC)
pdf(paste0(OutPutDirt,'/KnownSGBDoHeatmap.pdf'),width=HeatmapPDFWidth,height=HeatmapPDFHeight)
DoHeatmap(scDNA , features = top10markers$gene, size = 3)
dev.off()
write.table(top10markers,paste0(OutPutDirt,'/KnownSGBCell_top10markers.txt'),sep='\t',quote=F,row.names=F)
write.table(scDNA.markers,paste0(OutPutDirt,'/KnownSGBCell_allmarkers.txt'),sep='\t',quote=F,row.names=F)
cluster_info=scDNA$seurat_clusters
df=data.frame(Cluster=cluster_info)
df$Cell=rownames(df)
write.table(df,paste0(OutPutDirt,'/KnownSGBCell_ClusterCell.txt'),sep='\t',quote=F,row.names=F)
