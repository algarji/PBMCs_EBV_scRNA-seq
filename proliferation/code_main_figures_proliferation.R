library(dplyr)
library(tidyverse)
library(patchwork)
library(Seurat)
library(dittoSeq)
library(SignacX)
library(RColorBrewer)
library(ggthemes)
library(sceasy)
library(pheatmap)
library(escape)
library(CrossTalkeR)

# MERGING AND GENERAL MANIPULATION OF DATA

run2 <- readRDS("D:/EBV/10x/proliferation/run2.rds")	## Demultiplexed preprocessed run containing 2 Healthy samples and 1 Healthy seronegative sample
run3 <- readRDS("D:/EBV/10x/proliferation/run3.rds")	## Demultiplexed preprocessed run containing 1 XIAP y/- and 1 Healthy sample
run4 <- readRDS("D:/EBV/10x/proliferation/run4.rds")	## Demultiplexed preprocessed run containing 3 MAGT1y/- samples

pbmc.merged <- merge(run2, y = c(run3, run4), add.cell.ids = c("run2", "run3", "run4"), project = "blood")

pbmc.merged <- NormalizeData(pbmc.merged)
pbmc.merged <- FindVariableFeatures(pbmc.merged, selection.method = "mean.var.plot")
s.genes <- cc.genes$s.genes
g2m.genes <- cc.genes$g2m.genes
pbmc.merged <- CellCycleScoring(pbmc.merged, s.features = s.genes, g2m.features = g2m.genes, set.ident = TRUE)
pbmc.merged <- ScaleData(pbmc.merged, vars.to.regress=c('S.Score', 'G2M.Score','percent.mt'), features = rownames(pbmc.merged))
pbmc.merged <- RunPCA(pbmc.merged, features = VariableFeatures(object = pbmc.merged))
pbmc.merged <- FindNeighbors(pbmc.merged, dims = 1:30)
pbmc.merged <- FindClusters(pbmc.merged, resolution = 1.2)
pbmc.merged <- RunUMAP(pbmc.merged, dims = 1:30)

pbmc.markers <- FindAllMarkers(pbmc.merged, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
pbmc.markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC) -> top10
DoHeatmap(pbmc.merged, features = top10$gene) + NoLegend()

new.cluster.ids <- c("1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11", "12", "13", "14", "15", "16", "17", "18")
names(new.cluster.ids) <- levels(pbmc.merged)
pbmc.merged <- RenameIdents(pbmc.merged, new.cluster.ids)
number <- Idents(pbmc.merged)
number <- unname(number)
pbmc.merged@meta.data$number = number

# CHECK GENERAL ANNOTATION WITH SIGNACX AND AZIMUTH

labels_fast <- SignacFast(pbmc.merged)
celltypes = GenerateLabels(labels_fast, E = pbmc.merged)
lbls = factor(celltypes$CellTypes)
levels(lbls) <- sort(unique(lbls))
pbmc.merged <- AddMetaData(pbmc.merged, metadata = lbls, col.name = "celltypes_signacX")
DimPlot(pbmc.merged, group.by = "celltypes_signacX")

azimuth <- DietSeurat(pbmc.merged)
saveRDS(azimuth, file = "../azimuth.rds")	## run azimuth and download tsv predictions
azimuth_pred <- read.table('azimuth_pred.tsv', header=TRUE, sep='\t')
identical(azimuth_pred$cell, rownames(pbmc.merged@meta.data))
pbmc.merged@meta.data$azimuth = azimuth_pred$predicted.celltype.l2
DimPlot(pbmc.merged, group.by='azimuth', label=TRUE, label.size=2, repel=TRUE, cols=col_vector)
matrix_azimuth <- data.frame(row.names=azimuth_pred$cell, pred = azimuth_pred$predicted.celltype.l2.score)
matrix_azimuth <- t(matrix_azimuth)
pbmc.merged[["azimuth_pred"]] <- CreateAssayObject(counts = matrix_azimuth)
FeaturePlot(pbmc.merged, features = "azimuthpred_pred", cols = rev(brewer.pal(n = 11, name = "RdBu")))
pbmc.merged[["azimuth_pred"]] <- NULL

# ADD INFORMATION OF EBV GENES AND EBV GENOME

m <- read.csv('metadata_ebv_genes.csv', header=TRUE, sep=',')	## Normalized count matrix created when the dataset is aligned to a partially annotated EBV genome
row.names(m) <- m$X
common_1 <- intersect(row.names(m), row.names(pbmc.merged@meta.data))
barcode <- data.frame(X=row.names(pbmc.merged@meta.data))
m <- m[match(common_1, m$X),]
m <- merge(barcode, m[, c("X", "LMP.2B", "EBNA.2", "EBNA.3A", "EBNA.3B.EBNA.3C", "EBNA.1", "EBNA.LP", "BHRF1", "EBNA.1.1", "BZLF1", "BRLF1", "RPMS1", "A73", "LMP.2A", "LMP.1")], by="X", all=TRUE)
m[is.na(m)] <- 0
pbmc.merged@meta.data$LMP2B <- m$LMP.2B
pbmc.merged@meta.data$EBNA2 <- m$EBNA.2
pbmc.merged@meta.data$EBNA3A <- m$EBNA.3A
pbmc.merged@meta.data$EBNA3B_3C <- m$EBNA.3B.EBNA.3C
pbmc.merged@meta.data$EBNA1 <- m$EBNA.1
pbmc.merged@meta.data$EBNALP <- m$EBNA.LP
pbmc.merged@meta.data$BHRF1 <- m$BHRF1
pbmc.merged@meta.data$EBNA1.1 <- m$EBNA.1.1
pbmc.merged@meta.data$BZLF1 <- m$BZLF1
pbmc.merged@meta.data$BRLF1 <- m$BRLF1
pbmc.merged@meta.data$RPMS1 <- m$RPMS1
pbmc.merged@meta.data$A73 <- m$A73
pbmc.merged@meta.data$LMP2A <- m$LMP.2A
pbmc.merged@meta.data$LMP1 <- m$LMP.1

barcodes_ebv_genome <- read.table('barcode_EBV-genome.tsv', header=TRUE)	## list of cell barcodes with evidence of transcript alignment against the whole EBV genome
barcodes_ebv_genome$ebv <- "EBV"
common_2 <- intersect(barcodes_ebv_genome$barcode, row.names(pbmc.merged@meta.data))
barcode <- data.frame(barcode=row.names(pbmc.merged@meta.data))
barcodes_ebv_genome <- barcodes_ebv_genome[match(common_2, barcodes_ebv_genome$barcode),]
barcodes_ebv_genome <- merge(barcode, barcodes_ebv_genome[, c("barcode", "ebv")], by="barcode", all=TRUE)
identical(barcodes_ebv_genome$barcode, barcode$barcode)
pbmc.merged@meta.data$ebv_genome <- barcodes_ebv_genome$ebv

# FIGURE 2A

color <- as.list(ggthemes_data$tableau$`color-palettes`$regular$`Tableau 20`)
color <- as.character(color$value)
DimPlot(pbmc.merged, label=TRUE, cols= color)

# FIGURE 2B

DimPlot(pbmc.merged, group.by='status', cols = c('#29AF7FFF', 'orange1', '#F987C5', 'grey40'), shuffle=TRUE)

# FIGURE 2C

Idents(pbmc.merged) <- 'status'
magt1 <- subset(pbmc.merged, idents=c('MAGT1'))
xiap <- subset(pbmc.merged, idents=c('XIAP'))
seronegative <- subset(pbmc.merged, idents=c('Healthy_seronegative'))
healthy <- subset(pbmc.merged, idents=c('Healthy'))
Idents(magt1) <- 'ebv_genome'
Idents(xiap) <- 'ebv_genome'
Idents(seronegative) <- 'ebv_genome'
Idents(healthy) <- 'ebv_genome'
highlight1 <- WhichCells(healthy, idents=c("EBV"))
highlight2 <- WhichCells(seronegative, idents=c("EBV"))
highlight3 <- WhichCells(magt1, idents=c("EBV"))
highlight4 <- WhichCells(xiap, idents=c("EBV"))
pbmc.merged@meta.data$ebv_genome[is.na(pbmc.merged@meta.data$ebv_genome)] <- 'NO'
DimPlot(pbmc.merged, cells.highlight= list(highlight4, highlight3, highlight2, highlight1), cols.highlight = c('#29AF7FFF', 'orange1', '#F987C5', 'grey40'), cols= "grey", sizes.highlight = 1, shape.by = 'ebv_genome', na.value = 'grey') + scale_shape_manual(values = c("EBV" = 4, "NO"=16))

# FIGURE 2D

convertFormat(pbmc.merged, from="seurat", to="anndata", outFile='scanpy_input.h5ad')
# Python in Spyder in Anaconda
import scanpy as sc
adata =sc.read_h5ad("scanpy_input.h5ad")
markers={ 'B cells': ['MS4A1', 'CD19', 'CR2', 'FCRL1', 'ADAM28'], 'Plasma cells': ['MZB1', 'IGKC', 'IGHG1', 'IGLC3'], 'Naive T cells': ['CD3E', 'IL7R', 'CCR7', 'NELL2'], 'CD4 T cells': ['CD4', 'AQP3', 'GPR15', 'TNFRSF4', 'CSF2', 'FOXP3', 'TCF7'], 'Cytotoxic cells': ['FCGR3A', 'GNLY', 'CD8A', 'GZMK', 'TRDV2', 'CCL5', 'NCAM1', 'XCL1'],  'Monocytes': ['LYZ', 'CXCL8'], }
sc.pl.stacked_violin(adata, markers, 'number', dendrogram=False, categories_order=['2','4','6','16','11','13','14','1','12','3','7','15','18','5','9','10','17','8'], save='.pdf')

# FIGURE 2E

ebv <- subset(pbmc.merged, cells = common_1)
Idents(ebv) = 'status'
convertFormat(ebv, from="seurat", to="anndata", outFile='ebv.h5ad')
# Python in Spyder in Anaconda
import scanpy as sc
adata =sc.read_h5ad("ebv.h5ad")
markers={ 'Early': ['EBNA2', 'EBNA3A', 'EBNA3B_3C'], 'Latency II': ['EBNA1.1', 'LMP1', 'LMP2A', 'LMP2B'], 'Lytic': ['BRLF1'], 'BART': ['RPMS1'], }
sc.pl.dotplot(adata, markers, 'status', dendrogram=True, cmap='seismic', save='.pdf')

# FIGURE 2F

Idents(pbmc.merged) <- 'ebv_genome'
VlnPlot(pbmc.merged, features=c('CR2', 'IGHD'), cols=c('grey', 'darkred'))

# FIGURE 2G

cluster3.markers <- FindMarkers(pbmc.merged, ident.1 = '3', ident.2=c('7', '15', '18') min.pct = 0.25)
cluster3 <- filter(cluster3.markers, avg_log2FC > 0.58 & p_val_adj < 0.05)
cluster7.markers <- FindMarkers(pbmc.merged, ident.1 = '7', ident.2=c('3', '15', '18') min.pct = 0.25)
cluster7 <- filter(cluster7.markers, avg_log2FC > 0.58 & p_val_adj < 0.05)
cluster15.markers <- FindMarkers(pbmc.merged, ident.1 = '15', ident.2=c('3', '7', '18') min.pct = 0.25)
cluster15 <- filter(cluster15.markers, avg_log2FC > 0.58 & p_val_adj < 0.05)
cluster18.markers <- FindMarkers(pbmc.merged, ident.1 = '18', ident.2=c('3', '7', '15') min.pct = 0.25)
cluster18 <- filter(cluster18.markers, avg_log2FC > 0.58 & p_val_adj < 0.05)
dittoHeatmap(cd4, genes = c(c('S100A10', 'S100A6', 'KLF2', 'LGALS3', 'AQP3', 'CLIC3', 'CXCR6', 'CCR2', 'CCR5', 'HPGD', 'PXN', 'GZMA', 'SYTL2', 'GNLY', 'NKG7', 'ZBED2', 'TNFSF4', 'RDH10', 'CD83', 'BHLHE40', 'XCL1', 'SIAH2', 'TNFRSF4', 'CD200', 'CSF2', 'DUSP2', 'IRF4', 'MYB', 'PDCD1', 'TOX2', 'TCF7', 'PLAC8', 'CCR7', 'LEF1', 'NELL2', 'FOXP3', 'RTKN2', 'IL2RA', 'BCL2', 'ADTRP', 'CXCR4', 'SOCS3')), annot.by = c("RNA_snn_res.0.5", "status"), scaled.to.max = TRUE, order.by=c('RNA_snn_res.0.5', 'status'), complex = TRUE, cluster_rows=FALSE, heatmap.colors.max.scaled = colorRampPalette(c("white", "black"))(25))

# FIGURE 2H

gene.sets <- list(Tfh <- c('BCL6', 'BTLA', 'CD200', 'CD3D', 'CD3E', 'CD3G', 'CD4', 'CD40LG', 'CXCR3', 'CXCR5', 'ICA1', 'ICOS', 'IL21R', 'IL6ST', 'MAGEH1', 'PDCD1', 'PTPN13', 'SLAM', 'STAT3', 'TNFSF4', 'TOX', 'TOX2'))
cd4 <- subset(pbmc.merged, idents = c("3", "7", "15", "18"))
ES <- enrichIt(obj = cd4, gene.sets = gene.sets, min.size = NULL)
cd4 <- AddMetaData(cd4, ES)
cd4@meta.data$active.idents <- cd4@active.ident
ES2 <- data.frame(cd4[[]], Idents(cd4))
colnames(ES2)[ncol(ES2)] <- "cluster_escape"
ridgeEnrichment(ES2, gene.set = 'output', group = "cluster_escape", add.rug = TRUE)

#FIGURE 2I

cd8 <- subset(pbmc.merged, idents = c("9"))
cd8 <- DietSeurat(cd8)
Idents(cd8) <- 'status'
cd8 <- subset(cd8, idents = c("Healthy", "MAGT1"))
cd8 <- ScaleData(cd8, vars.to.regress = c("S.Score", "G2M.Score"), features = rownames(cd8), block.size = 100)
cd8 <- RunPCA(cd8, verbose = FALSE)
ElbowPlot(cd8)
cd8 <- FindNeighbors(cd8, dims = 1:18)
cd8 <- FindClusters(cd8, resolution = 0.5)
cd8 <- RunUMAP(cd8, dims = 1:18)
new.cluster.ids <- c("1", "2", "3", "4", "5")
names(new.cluster.ids) <- levels(cd8)
cd8 <- RenameIdents(cd8, new.cluster.ids)
cluster <- Idents(cd8)
cluster <- unname(cluster)
cd8@meta.data$cluster = cluster
DimPlot(cd8, cols = c('#E31A1C', '#762A83', '#08519C', '#FF7F00', '#33A02C'), shuffle=TRUE)

# FIGURE 2J

DimPlot(cd8, group.by='status', cols = c('#29AF7FFF', '#F987C5'), shuffle=TRUE))

# FIGURE 2K

cd8.markers <- FindAllMarkers(cd8, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
cd8.markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC) -> top10
DoHeatmap(cd8, features = top10$gene, group.colors = c('#E31A1C', '#762A83', '#08519C', '#FF7F00', '#33A02C')) + NoLegend() + scale_fill_gradientn(colors = c('white',"white", "#D33682"))

# FIGURE 2L

Idents(cd8) <- 'status'
cd8.status.markers <- FindAllMarkers(cd8, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
cd8.status.markers <- filter(cd8.status.markers, avg_log2FC > 0.58)
pseudo <- AverageExpression(object = cd8, group.by='status', features=cd8.status.markers$gene, assays = 'RNA', slot='scale.data')
pseudo <- pseudo$RNA
pseudo <- filter(pseudo, Healthy > 0.1 | MAGT1 > 0.1)
pheatmap(pseudo, show_rownames = FALSE, border_color = 'black', cluster_cols = FALSE, cluster_rows=FALSE, cellwidth = 15, fontsize_row = 5)

# FIGURE 2M

Idents(cd8) <- 'status'
cd8_healthy <- subset(cd8, idents=c('Healthy'))
FeaturePlot(cd8_healthy, features=c('GNLY'), order=TRUE) + scale_color_gradientn(colours=viridis::magma(100))
cd8_magt1 <- subset(cd8, idents=c('MAGT1'))
FeaturePlot(cd8_magt1, features=c('GNLY'), order=TRUE) + scale_color_gradientn(colours=viridis::magma(100))

# FIGURE 2N

nk <- subset(pbmc.merged, idents=c('5', '17'))
dittoBarPlot(nk, 'number', 'status', color.panel = c('#B07AA1', '#59A14F')) + theme_minimal()

# FIGURE 2O

VlnPlot(nk, features=c('CX3CR1', 'FGFBP2', 'FCGR3A', 'NCAM1', 'COTL1', 'XCL1'), cols=c('#59A14F', '#B07AA1'))

# FIGURE 2P

new.cluster.ids <- c("cluster_1", "cluster_2", "cluster_3", "cluster_4", "cluster_5", "cluster_6", "cluster_7", "cluster_8", "cluster_9", "cluster_10", "cluster_11", "cluster_12", "cluster_13", "cluster_14", "cluster_15", "cluster_16", "cluster_17", "cluster_18")
names(new.cluster.ids) <- levels(pbmc.merged)
pbmc.merged <- RenameIdents(pbmc.merged, new.cluster.ids)
cluster <- Idents(pbmc.merged)
pbmc.merged@meta.data$cluster = cluster
Idents(pbmc.merged) <- 'status'
magt1_subset <- subset(pbmc.merged, idents=c('MAGT1'))
Idents(magt1_subset) <- 'number'
magt1_subset <- subset(magt1_subset, idents=c(1, 2, 4, 6, 7, 9, 10, 15, 16, 17))	## retain clusters with more than 20 cells
matrix_magt1 <- as.data.frame(magt1_subset@assays$RNA@data)
matrix_magt1 <- matrix_magt1[rowSums(matrix_magt1[,2:dim(matrix)[2]])!=0,]
xiap_subset <- subset(pbmc.merged, idents=c('XIAP'))
Idents(xiap_subset) <- 'number'
xiap_subset <- subset(xiap_subset, idents=c(2, 3, 4, 7, 11, 13, 14, 15, 18))	## retain clusters with more than 20 cells
matrix_xiap <- as.data.frame(xiap_subset@assays$RNA@data)
matrix_xiap <- matrix_xiap[rowSums(matrix_xiap[,2:dim(matrix)[2]])!=0,]
allgenes <- rownames(pbmc.merged)
matrix <- as.data.frame(pbmc.merged@assays$RNA@data)
matrix <- matrix[rowSums(matrix[,2:dim(matrix)[2]])!=0,]
s1 <- grepl('Healthy',pbmc.merged@meta.data$status)
s2 <- grepl('MAGT1',magt1_subset@meta.data$status)
s3 <- grepl('XIAP',xiap_subset@meta.data$status)
s1[match('gene',colnames(matrix))] <- TRUE
s2[match('gene',colnames(matrix_magt1))] <- TRUE
s3[match('gene',colnames(matrix_xiap))] <- TRUE
write.table(matrix[,s1], 's1_filtered_hcount.csv',row.names=T,sep=',')
write.table(matrix_magt1[,s2], 's2_filtered_hcount.csv',row.names=T,sep=',')
write.table(matrix_xiap[,s3], 's3_filtered_hcount.csv',row.names=T,sep=',')
metadata_s1 <- data.frame(cells=rownames(pbmc.merged@meta.data[grepl('Healthy',pbmc.merged@meta.data$status),]),cluster=pbmc.merged@meta.data$cluster[grepl('Healthy',pbmc.merged@meta.data$status)])
metadata_s2 <- data.frame(cells=rownames(magt1_subset@meta.data[grepl('MAGT1',pbmc.merged@meta.data$status),]),cluster=magt1_subset@meta.data$cluster[grepl('MAGT1',magt1_subset@meta.data$status)])
metadata_s3 <- data.frame(cells=rownames(xiap_subset@meta.data[grepl('XIAP',xiap_subset@meta.data$status),]),cluster=xiap_subset@meta.data$cluster[grepl('XIAP',xiap_subset@meta.data$status)])
write.csv(metadata_s1, 's1_filtered_meta.csv', row.names=FALSE)
write.csv(metadata_s2, 's2_filtered_meta.csv', row.names=FALSE)
write.csv(metadata_s3, 's3_filtered_meta.csv', row.names=FALSE)
# Linux in BASH
cd /mnt/d/EBV/10x/proliferation/
conda activate cpdb
mkdir s1 s2 s3
cellphonedb method statistical_analysis s1_filtered_meta.csv  s1_filtered_hcount.csv --counts-data=gene_name --output-path s1/
cellphonedb method statistical_analysis s2_filtered_meta.csv  s2_filtered_hcount.csv --counts-data=gene_name --output-path s2/
cellphonedb method statistical_analysis s3_filtered_meta.csv  s3_filtered_hcount.csv --counts-data=gene_name --output-path s3/
# Python in Spyder in Anaconda
def correct_lr(data):
    import pandas as pd
    def swap(a,b): return b,a
    data=data.to_dict('index')
    for k,v in data.items():
        if v['isReceptor_fst'] and v['isReceptor_scn']:
            v['isReceptor_fst'],v['isReceptor_scn'] = swap(v['isReceptor_fst'],v['isReceptor_scn'])
            v['Ligand'],v['Receptor'] = swap(v['Ligand'],v['Receptor'])
            v['Ligand.Cluster'],v['Receptor.Cluster'] = swap(v['Ligand.Cluster'],v['Receptor.Cluster'])
        elif v['isReceptor_fst'] and not v['isReceptor_scn']:
                v['isReceptor_fst'],v['isReceptor_scn'] = swap(v['isReceptor_fst'],v['isReceptor_scn'])
                v['Ligand'],v['Receptor'] = swap(v['Ligand'],v['Receptor'])
                v['Ligand.Cluster'],v['Receptor.Cluster'] = swap(v['Ligand.Cluster'],v['Receptor.Cluster'])
    res_df = pd.DataFrame.from_dict(data,orient='index')
    return(res_df)
def cpdb2df(data,clsmapping):
    data = data.fillna(0)
    df_data = {}
    df_data['Ligand'] = []
    df_data['Receptor'] = []
    df_data['Ligand.Cluster'] = []
    df_data['Receptor.Cluster'] = []
    df_data['isReceptor_fst'] = []
    df_data['isReceptor_scn'] = []
    df_data['MeanLR'] = []
    for i in range(data.shape[0]):
        pair = list(data['interacting_pair'])[i].split('_')
        for j in range(data.iloc[:,12:].shape[1]):
            c_pair = list(data.columns)[j+12].split('|')
            if float(data.iloc[i,j+12]) !=0.0:
                df_data['Ligand'].append(pair[0])
                df_data['Receptor'].append(pair[1])
                df_data['Ligand.Cluster'].append(clsmapping[c_pair[0]])
                df_data['Receptor.Cluster'].append(clsmapping[c_pair[1]])
                df_data['isReceptor_fst'].append(list(data['receptor_a'])[i])
                df_data['isReceptor_scn'].append(list(data['receptor_b'])[i])
                df_data['MeanLR'].append(data.iloc[i,j+12])
    data_final = pd.DataFrame.from_dict(df_data)
    return(data_final)
import pandas as pd
s1 = pd.read_csv('./s1/significant_means.txt',sep='\t')
s2 = pd.read_csv('./s2/significant_means.txt',sep='\t')
num_to_clust = {'cluster_1':'1', 'cluster_2':'2', 'cluster_3':'3', 'cluster_4':'4', 'cluster_5':'5', 'cluster_6':'6', 'cluster_7':'7', 'cluster_8':'8', 'cluster_9':'9', 'cluster_10':'10', 'cluster_11':'11', 'cluster_12':'12', 'cluster_13':'13', 'cluster_14':'14', 'cluster_15':'15', 'cluster_16':'16', 'cluster_17':'17', 'cluster_18':'18'}
s1_filtered = cpdb2df(s1,num_to_clust)
s2_filtered = cpdb2df(s2,num_to_clust)
s3_filtered = cpdb2df(s3,num_to_clust)
s1_filtered = correct_lr(s1_filtered)
s2_filtered = correct_lr(s2_filtered)
s3_filtered = correct_lr(s3_filtered)
s1_filtered_final =pd.DataFrame()
s1_filtered_final['source'] = s1_filtered['Ligand.Cluster']
s1_filtered_final['target'] = s1_filtered['Receptor.Cluster']
s1_filtered_final['gene_A'] = s1_filtered['Ligand']
s1_filtered_final['gene_B'] = s1_filtered['Receptor']
s1_filtered_final['type_gene_A'] = 'Ligand'
s1_filtered_final['type_gene_B']= 'Receptor'
s1_filtered_final['MeanLR']=s1_filtered['MeanLR']
s2_filtered_final =pd.DataFrame()
s2_filtered_final['source'] = s2_filtered['Ligand.Cluster']
s2_filtered_final['target'] = s2_filtered['Receptor.Cluster']
s2_filtered_final['gene_A'] = s2_filtered['Ligand']
s2_filtered_final['gene_B'] = s2_filtered['Receptor']
s2_filtered_final['type_gene_A'] = 'Ligand'
s2_filtered_final['type_gene_B']= 'Receptor'
s2_filtered_final['MeanLR']=s2_filtered['MeanLR']
s3_filtered_final =pd.DataFrame()
s3_filtered_final['source'] = s3_filtered['Ligand.Cluster']
s3_filtered_final['target'] = s3_filtered['Receptor.Cluster']
s3_filtered_final['gene_A'] = s3_filtered['Ligand']
s3_filtered_final['gene_B'] = s3_filtered['Receptor']
s3_filtered_final['type_gene_A'] = 'Ligand'
s3_filtered_final['type_gene_B']= 'Receptor'
s3_filtered_final['MeanLR']=s3_filtered['MeanLR']
s1_filtered_final.to_csv('s1_filtered_corrected.csv')
s2_filtered_final.to_csv('s2_filtered_corrected.csv')
s3_filtered_final.to_csv('s3_filtered_corrected.csv')
# R in RStudio
s1_filtered_corrected <- read.csv('s1_filtered_corrected.csv', header=TRUE)
s2_filtered_corrected <- read.csv('s2_filtered_corrected.csv', header=TRUE)
s3_filtered_corrected <- read.csv('s3_filtered_corrected.csv', header=TRUE)
paths <- c('MAGT1' = 's2_filtered_corrected.csv','Healthy' = 's1_filtered_corrected.csv')
output =  'D:/EBV/10x/proliferation/'
data <- generate_report(paths, out_path=paste0(output,'/'), threshold=0, out_file = 'trial.html', output_fmt = "html_document", report = FALSE)
plot_cci(graph = data@graphs$Healthy, colors = data@colors, plt_name = 'Differential cell-cell interaction network', coords = data@coords, emax = NULL, leg = FALSE, low = 0, high = 0, ignore_alpha = FALSE, log = FALSE, efactor = 8, vfactor = 12)
plot_cci(graph = data@graphs$MAGT1, colors = data@colors, plt_name = 'Differential cell-cell interaction network', coords = data@coords, emax = NULL, leg = FALSE, low = 0, high = 0, ignore_alpha = FALSE, log = FALSE, efactor = 8, vfactor = 12)
paths <- c('XIAP' = 's3_filtered_corrected.csv','Healthy' = 's1_filtered_corrected.csv')
output =  'D:/EBV/10x/proliferation/'
data <- generate_report(paths, out_path=paste0(output,'/'), threshold=0, out_file = 'trial.html', output_fmt = "html_document", report = FALSE)
plot_cci(graph = data@graphs$XIAP, colors = data@colors, plt_name = 'Differential cell-cell interaction network', coords = data@coords, emax = NULL, leg = FALSE, low = 0, high = 0, ignore_alpha = FALSE, log = FALSE, efactor = 8, vfactor = 12)
