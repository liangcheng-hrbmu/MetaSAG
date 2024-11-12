
###############使用R包Seurat对所有KSGB的细胞进行聚类###############
###################################################################

library(Seurat)
library(dplyr)
library(getopt)

spec <- matrix(
c("CellUnirefCount","i", 2, "character", "This is All SGB Cell's Uniref Count Matrix File.",
  "OutPut", "o", 2, "character",  "This is th OutPutFolder"),byrow=TRUE,ncol=5)

opt = getopt(spec)
CountMatrixFile = opt$CellUnirefCount
OutPutDirt = opt$OutPut

print(CountMatrixFile)
print(OutPutDirt)



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
scDNA <- FindVariableFeatures(scDNA, selection.method = 'vst', nfeatures = 20000)
scDNA <- ScaleData(scDNA)
scDNA=RunPCA(scDNA,npcs=30,features=VariableFeatures(scDNA))
scDNA <- FindNeighbors(scDNA, dims = 1:10)
scDNA <- FindClusters(scDNA, resolution = 0.5)
table(scDNA$seurat_clusters)
scDNA <- RunUMAP(scDNA, dims = 1:10)
pdf(paste0(OutPutDirt,"/KnownSGBUmap.pdf"), width = 6, height = 4)
DimPlot(scDNA, reduction = 'umap',label = TRUE)
dev.off()
scDNA.markers <- FindAllMarkers(scDNA, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
top10markers<-scDNA.markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC)
pdf(paste0(OutPutDirt,'/KnownSGBDoHeatmap.pdf'),width=12,height=8)
DoHeatmap(scDNA , features = top10markers$gene, size = 3)
dev.off()
write.table(top10markers,paste0(OutPutDirt,'/KnownSGBCell_top10markers.txt'),sep='\t',quote=F,row.names=F)
write.table(scDNA.markers,paste0(OutPutDirt,'/KnownSGBCell_allmarkers.txt'),sep='\t',quote=F,row.names=F)
cluster_info=scDNA$seurat_clusters
df=data.frame(Cluster=cluster_info)
df$Cell=rownames(df)
write.table(df,paste0(OutPutDirt,'/KnownSGBCell_ClusterCell.txt'),sep='\t',quote=F,row.names=F)



