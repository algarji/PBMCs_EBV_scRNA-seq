library(dplyr)
library(tidyverse)
library(patchwork)
library(Seurat)
library(dittoSeq)
library(viridis)
library(RColorBrewer)
library(ggthemes)
library(EnhancedVolcano)
library(sceasy)
library(immunarch)

# INTEGRATION AND GENERAL MANIPULATION OF DATA

run1 <- readRDS("D:/EBV/10x/activation/run1.rds")	## Demultiplexed preprocessed run containing 3 Healthy samples
run5 <- readRDS("D:/EBV/10x/activation/run5.rds")	## Demultiplexed preprocessed run containing 1 XIAPy/- and 2 PIK3CDGOF samples
run6 <- readRDS("D:/EBV/10x/activation/run6.rds")	## Demultiplexed preprocessed run containing 2 MAGT1y/- samples
run7 <- readRDS("D:/EBV/10x/activation/run7.rds")	## Demultiplexed preprocessed run containing 3 Healthy samples
run8 <- readRDS("D:/EBV/10x/activation/run8.rds")	## Demultiplexed preprocessed run containing 3 Healthy samples

pbmc.list <- list("run1" = run1, "run5" = run5, "run6" = run6, "run7" = run7, "run8" = run8)
pbmc.list <- lapply(X = pbmc.list, FUN = function(x) {
    x <- NormalizeData(x)
    x <- FindVariableFeatures(x, selection.method = "vst", nfeatures = 2000)
})
features <- SelectIntegrationFeatures(object.list = pbmc.list, nfeatures = 2000)
immune.anchors <- FindIntegrationAnchors(object.list = pbmc.list, anchor.features = features)
pbmc.combined <- IntegrateData(anchorset = immune.anchors)

DefaultAssay(pbmc.combined) <- "integrated"
pbmc.combined <- ScaleData(pbmc.combined, verbose = FALSE)
pbmc.combined <- RunPCA(pbmc.combined, npcs = 30, verbose = FALSE)
pbmc.combined <- RunUMAP(pbmc.combined, reduction = "pca", dims = 1:30)
pbmc.combined <- FindNeighbors(pbmc.combined, reduction = "pca", dims = 1:30)
pbmc.combined <- FindClusters(pbmc.combined, resolution = 0.5)

pbmc.markers <- FindAllMarkers(pbmc.combined, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
pbmc.markers %>% group_by(cluster) %>% top_n(n = 5, wt = avg_log2FC) -> top5
DoHeatmap(pbmc.combined, features = top5$gene) + NoLegend()

new.cluster.ids <- c("1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11", "12", "13", "14", "15", "16", "17", "18", "19", "20", "21")
names(new.cluster.ids) <- levels(pbmc.combined)
pbmc.combined <- RenameIdents(pbmc.combined, new.cluster.ids)
number <- Idents(pbmc.combined)
number <- unname(number)
pbmc.combined@meta.data$number = number

# FIGURE 3A

color <- c('#1A4870', '#92353C', '#567630', '#E37F07', '#6F8E84', '#C21537', '#928DC4', '#CAC27F', '#A1532E', '#7B92A8', '#2F6E67', '#9C8846', '#BFA19C', '#FFD207', '#D9E6EB', '#4DAD33', '#9A9999', '#624193', '#ADD8E5', '#CE929E', '#E84B1E')
DimPlot(pbmc.combined, label=TRUE, cols=color)

# FIGURE 3B

DimPlot(pbmc.combined, split.by='Mutation', cols=c('#31AF80', '#EC8BB8', '#F29B4B', '#688B8C', '#B6E0EE'))

# Subset CD4 T cells

cd4 <- subset(pbmc.combined,idents=c('3', '4', '5', '7', '19'))
DefaultAssay(cd4) <- 'RNA'
cd4 <- DietSeurat(cd4)
s.genes <- cc.genes$s.genes
g2m.genes <- cc.genes$g2m.genes
cd4 <- CellCycleScoring(cd4, s.features = s.genes, g2m.features = g2m.genes, set.ident = TRUE)
cd4 <- ScaleData(cd4, vars.to.regress = c("S.Score", "G2M.Score"), features = rownames(cd4), block.size = 100)
cd4 <- RunPCA(cd4, verbose = FALSE)
cd4 <- FindNeighbors(cd4, dims = 1:20)
cd4 <- FindClusters(cd4, resolution = 0.5)
cd4 <- RunUMAP(cd4, dims = 1:20)
cd4 <- subset(cd4, idents=c('1', '3', '4', '6', '7', '9', '10', '11', '12'))	## Subset only non-naive CD4 T cells (CD3+CD4+CCR7-)
cd4 <- FindClusters(cd4, resolution = 0.5)
cd4 <- RunUMAP(cd4, dims = 1:20)
new.cluster.ids <- c("1", "2", "3", "4", "5", "6", "7")
names(new.cluster.ids) <- levels(cd4)
cd4 <- RenameIdents(cd4, new.cluster.ids)
cluster <- Idents(cd4)
cluster <- unname(cluster)
cd4@meta.data$cluster = cluster
Idents(cd4) <- 'Mutation'
cd4 <- subset(cd4, idents=c('WT', 'MAGT1', 'XIAP'))	## retain only highly representative conditions
Idents(cd4) <- 'cluster'

# FIGURE 4A

DimPlot(cd4, label=TRUE, cols=c('#869FB1', '#658C69', '#EF7845', '#F8B154', '#F08CAB', '#80776F', '#8E569E'))

# FIGURE 4B

DimPlot(cd4, label=TRUE, cols=c('#869FB1', '#658C69', '#EF7845', '#F8B154', '#F08CAB', '#80776F', '#8E569E'), split.by='Mutation')

# FIGURE 4C

clusters_reg <- FindMarkers(cd4, ident.1 = "6", ident.2 = "4", min.pct = 0.25)
volcano_reg <- data.frame(row.names=row.names(clusters_reg), log2FoldChange=clusters_reg$avg_log2FC, pvalue=clusters_reg$p_val_adj)
EnhancedVolcano(volcano_reg, lab = rownames(volcano_reg), x = 'log2FoldChange', y = 'pvalue', title = 'SC4_vs_SC6', pCutoff = 10e-4, FCcutoff = 0.58, pointSize = 3.0, labSize = 2.0, col=c('black', 'black', 'black', 'red3'), colAlpha = 1)

# FIGURE 4D

markers <- c('CXCR3', 'PRF1', 'GZMA', 'GZMK', 'CCL5', 'KLRB1', 'CEBPD', 'TNF', 'CCR6', 'IL7R', 'IFITM1', 'IFIT2', 'ISG15', 'SAT1', 'TNFSF10', 'MYLIP', 'CD69', 'GNLY', 'GZMB', 'ICOS', 'OTOF', 'CCR2', 'PDCD1', 'TIGIT', 'TNFRSF18', 'CSF2', 'IFNG')
DotPlot(cd4, features = markers) + xlab('Gene') +  ylab('Cluster') + geom_point(aes(size=pct.exp), shape = 21, colour="black", stroke=0.5) + scale_colour_viridis(option="magma") + theme(axis.text.x = element_text(angle = 45, hjust=1)) + guides(size=guide_legend(override.aes=list(shape=21, colour="black", fill="white")))

# FIGURE 4E

sceasy::convertFormat(cd4, from="seurat", to="anndata", outFile='scanpy_cd4.h5ad')
# Python in Spyder in Anaconda
import numpy as np
import pandas as pd
import scanpy as sc
import scirpy as ir
from matplotlib import pyplot as plt
from matplotlib import cm as mpl_cm
from cycler import cycler
adata_tcr = ir.io.read_10x_vdj("filtered_contig_annotations_activation.csv")
adata =sc.read_h5ad("scanpy_cd4.h5ad")
ir.pp.merge_with_ir(adata, adata_tcr)
ir.tl.chain_qc(adata)
adata = adata[~adata.obs["chain_pairing"].isin(["orphan VDJ", "orphan VJ", "extra VJ", "extra VDJ", "two full chains"]), :].copy()
ir.tl.clonal_expansion(adata)
sc.pl.umap(adata, color=["clonal_expansion"], save='expansion.pdf')

# Subset CD8 T cells

cd8 <- subset(pbmc.combined,idents=c('1', '5', '8', '14', '16'))
DefaultAssay(cd8) <- 'RNA'
cd8 <- DietSeurat(cd8)
cd8 <- CellCycleScoring(cd8, s.features = s.genes, g2m.features = g2m.genes, set.ident = TRUE)
cd8 <- ScaleData(cd8, vars.to.regress = c("S.Score", "G2M.Score"), features = rownames(cd8), block.size = 100)
cd8 <- RunPCA(cd8, verbose = FALSE)
cd8 <- FindNeighbors(cd8, dims = 1:20)
cd8 <- FindClusters(cd8, resolution = 0.5)
cd8 <- RunUMAP(cd8, dims = 1:20)
cd8 <- subset(cd8, idents=c('1', '2', '5', '7', '8', '9', '10', '11', '13'))	## Subset only non-naive CD8 T cells (CD3+CD8ab+CCR7-)
cd8 <- FindClusters(cd8, resolution = 0.5)
cd8 <- RunUMAP(cd8, dims = 1:20)
new.cluster.ids <- c("1", "2", "3", "4", "5", "6", "7", "8", "9", "10")
names(new.cluster.ids) <- levels(cd8)
cd8 <- RenameIdents(cd8, new.cluster.ids)
cluster <- Idents(cd8)
cluster <- unname(cluster)
cd8@meta.data$cluster = cluster

# FIGURE 4F

DimPlot(cd8, label=TRUE, cols=c('#207AB7', '#2DA237', '#D9282A', '#AEC7E8', '#D37BB1', '#8F68AA', '#8F594E', '#F29699', '#7F8181', '#F9BB7B'))

# FIGURE 4G

dittoBarPlot(cd8, "cluster", group.by = "Mutation", scale = "count", colors=c('#207AB7', '#2DA237', '#D9282A', '#AEC7E8', '#D37BB1', '#8F68AA', '#8F594E', '#F29699', '#7F8181', '#F9BB7B'), x.reorder=c(1, 2, 3, 4, 5, 6, 7, 8, 9, 10))

# FIGURE 4H

cd8_subset <- subset(cd8, idents=c('1', '2', '3'))
DoHeatmap(cd8_subset, features = c('IL7R', 'KLRG1', 'DUSP2', 'XCL1', 'FCRL3', 'GZMK', 'CD27', 'SH2D1A', 'SELL', 'TIGIT', 'CX3CR1', 'ITGB1', 'GZMB', 'TBX21', 'CLIC3', 'ASCL2', 'GZMH', 'GNLY', 'ZNF683', 'FGFBP2', 'ADGRG1', 'FCRL6', 'KLRD1', 'KLRC2', 'KLRC1', 'CCR5', 'CD82', 'FAM111B', 'CD70', 'FAM207A', 'CD28', 'TOMM40', 'CCR1', 'IL2RA', 'NME1', 'LGALS1', 'CARHSP1', 'HLA-DRA', 'GPI', 'APOBEC3H', 'LGALS3', 'HOPX', 'CD38', 'FABP5', 'DUSP4', 'LTB'), group.colors = c('#207AB7', '#2DA237', '#D9282A')) + NoLegend()

# FIGURE 4I

proliferation <- readRDS("D:/EBV/10x/proliferation/pbmc_merged.rds")
Idents(proliferation) <- 'status'
proliferation <- subset(proliferation, idents=c('MAGT1'))
Idents(proliferation) <- 'number'
xmen_cd8_proliferation <- subset(proliferation, idents=c('9'))
xmen_cd8_proliferation@meta.data$exp <- 'PROLIFERATION'
Idents(cd8_subset) <- 'Mutation'
xmen_cd8_activation <- subset(cd8_subset, idents=c('MAGT1'))
xmen_cd8_activation@meta.data$exp <- 'ACTIVATION'
xmen_cd8 <- merge(xmen_cd8_proliferation, y = c(xmen_cd8_activation), add.cell.ids = c("proliferation", "activation"), project = "blood")
RidgePlot(xmen_cd8, features = c('EOMES', 'CD27', 'GZMK', 'KLRK1', 'GZMB', 'HLA-DRA', 'DUSP2', 'KLRD1', 'GNLY', 'KLRC1', 'KLRC2', 'ZNF683'), ncol = 6, group.by = 'exp')

# FIGURE 4J

cd8_clonal <- subset(cd8, idents=c('5', '6', '7', '8', '9', '10'))
sceasy::convertFormat(cd8_clonal, from="seurat", to="anndata", outFile='scanpy_cd8_clonal.h5ad')
# Python in Spyder in Anaconda
import scanpy as sc
adata =sc.read_h5ad("scanpy_cd8_clonal.h5ad")
markers=['TRAV6', 'TRAV26-2', 'TRAV12-1', 'TRAV8-6', 'TRAV27', 'TRAV1-1', 'TRBV20-1', 'TRBV7-6', 'TRBV2', 'TRBV7-2', 'TRBV28']
sc.pl.matrixplot(adata, markers, 'cluster', cmap='purples', standard_scale='var', colorbar_title='column scaled\nexpression', save='matrix_clonal.pdf')

# FIGURE 4K

file_path = "D:/EBV/10x/activation/clonotypes/immunarch/"	## Metadata file and independent CSV files of "filtered_contigs_annotations" per sample 
immdata <- repLoad(file_path)
exp_vol <- repExplore(immdata$data, .method = "volume")
exp_cnt <- repExplore(immdata$data, .method = "count")

# FIGURE 4L

imm_gu <- geneUsage(immdata$data, "hs.trbv", norm= T)
imm_cl_pca <- geneUsageAnalysis(imm_gu, "js+pca+kmeans", .verbose = F)
vis(imm_cl_pca, .plot = "clust")

# FIGURE 4M

sceasy::convertFormat(cd8, from="seurat", to="anndata", outFile='scanpy_cd8.h5ad')
# Python in Spyder in Anaconda
import numpy as np
import pandas as pd
import scanpy as sc
import scirpy as ir
from matplotlib import pyplot as plt
from matplotlib import cm as mpl_cm
from cycler import cycler
adata_tcr = ir.io.read_10x_vdj("filtered_contig_annotations_activation.csv")
adata =sc.read_h5ad("scanpy_cd8.h5ad")
ir.pp.merge_with_ir(adata, adata_tcr)
ir.tl.chain_qc(adata)
adata = adata[~adata.obs["chain_pairing"].isin(["orphan VDJ", "orphan VJ", "extra VJ", "extra VDJ", "two full chains"]), :].copy()
ir.tl.clonal_expansion(adata)
vdjdb = ir.datasets.vdjdb()
ir.pp.ir_dist(adata, vdjdb, metric="identity", sequence="aa")
ir.tl.ir_query(adata, vdjdb, metric="identity", sequence="aa", receptor_arms="any", dual_ir="any")
ir.tl.ir_query_annotate_df(adata, vdjdb, metric="identity", sequence="aa", include_ref_cols=["antigen.species", "antigen.gene"],).tail()
ir.tl.ir_query_annotate(adata, vdjdb, metric="identity", sequence="aa", include_ref_cols=["antigen.species"], strategy="most-frequent",)
sc.pl.umap(adata, color="antigen.species", color= 'EBV', save='umap_VDJdb.pdf')

# FIGURE 5A

# Subset NK cells and annotate clusters
Idents(nk) <- 'Mutation'
nk <- subset(nk, idents=c('WT', 'MAGT1', 'PIK3CD_2'))	## retain only highly representative conditions
Idents(nk) <- 'cluster'
table(nk@meta.data$cluster, nk@meta.data$Mutation)	## Represented in Prism

# FIGURE 5B

nk@meta.data$adapt <- 'Rest'
metadata <- nk@meta.data
metadata <- metadata %>% mutate(adapt = case_when(endsWith(cluster, "2") ~ "Adaptive"))
Idents(nk) <- 'adapt'
sceasy::convertFormat(nk, from="seurat", to="anndata", outFile='scanpy_nk.h5ad')
# Python in Spyder in Anaconda
adata =sc.read_h5ad("scanpy_nk.h5ad")
markers = ['KLRC2', 'GZMH', 'LAG3', 'HLA-DRB1', 'FCRL6', 'PRSS23', 'ITM2A', 'FCER1G', 'KLRC1', 'KLRB1', 'NCR3', 'XCL1', 'CCL3', 'LTB']
sc.tl.rank_genes_groups(adata, groupby='cluster', method='wilcoxon')
sc.pl.rank_genes_groups_violin(adata, gene_names=markers, jitter= False, strip = False)

# FIGURE 5C

Idents(nk) <- 'cluster'
nk_dim <- subset(nk, idents=c('1', '3'))
Idents(nk_dim) <- 'Mutation'
sceasy::convertFormat(nk_dim, from="seurat", to="anndata", outFile='scanpy_nk_dim.h5ad')
# Python in Spyder in Anaconda
adata =sc.read_h5ad("scanpy_nk_dim.h5ad")
markers = ['KIR2DL1', 'KIR2DL3', 'KIR3DL1', 'KIR3DL2', 'KIR2DL4', 'KIR3DL3', 'FCGR3A', 'KLRC1', 'KLRK1']
sc.pl.matrixplot(adata, markers, 'Mutation', cmap='seismic', standard_scale='group', colorbar_title='column scaled\nexpression', save='matrix_dim.pdf')

# FIGURE 5D

nk_bright <- subset(nk, idents=c('4', '5'))
Idents(nk_bright) <- 'Mutation'
VlnPlot(nk_bright, features=c('CD38', 'STMN1', 'IL7R', 'JAML'), cols=c('#EC8BB8', '#31AF80'))

# FIGURE 5E

# Subset Vd2 cells and annotate clusters
Vd2_eff <- subset(Vd2, idents=c('2', '3'))
Idents(Vd2_eff) <- 'Mutation'
DotPlot(Vd2_eff, features = c('TRDV2', 'CD27', 'CD8A', 'TIGIT', 'CD69', 'TNF', 'CCL3', 'CCL4', 'XCL1', 'FASLG', 'GZMB', 'TNFSF10', 'TNFSF14'), cols = c("darkred", "darkblue"), dot.scale = 8, split.by = "cluster") + RotatedAxis()

# FIGURE 5F

# Subset Vd1/3 cells and annotate clusters
Vd1_cyto <- subset(Vd1, idents=c('2'))
Idents(Vd1_cyto) <- 'Mutation'
sceasy::convertFormat(Vd1_cyto, from="seurat", to="anndata", outFile='scanpy_Vd1.h5ad')
# Python in Spyder in Anaconda
adata =sc.read_h5ad("scanpy_Vd1.h5ad")
markers = ['IL2RA', 'HLA-DRA', 'CD27', 'CD38', 'NCR3', 'ITGB7', 'CXCR3', 'CD82', 'TFRC', 'SLAMF1', 'PECAM1', 'CX3CR1', 'KLRB1', 'KLRD1', 'KLRG1', 'IL7R', 'GZMB', 'GNLY']
sc.pl.heatmap(adata, markers, groupby='Mutation', cmap='viridis', standard_scale='var')

# FIGURE 5G

Vd3 <- subset(Vd1, idents=c('3'))
VlnPlot(Vd3, features=c('TRDV3', 'KIR2DL3', 'KLRB1', 'FCGR3A', 'CCL3', 'CCL4L2', 'KLRF1', 'GZMB', 'GZMH'), cols=c('#99488F'))

# FIGURE 5H

# Subset MAIT cells and annotate clusters
DimPlot(mait, label=TRUE)

# FIGURE 5I

sceasy::convertFormat(mait, from="seurat", to="anndata", outFile='scanpy_mait.h5ad')
# Python in Spyder in Anaconda
import numpy as np
import pandas as pd
import scanpy as sc
import scirpy as ir
from matplotlib import pyplot as plt
from matplotlib import cm as mpl_cm
from cycler import cycler
adata_tcr = ir.io.read_10x_vdj("filtered_contig_annotations_activation.csv")
adata =sc.read_h5ad("scanpy_mait.h5ad")
ir.pp.merge_with_ir(adata, adata_tcr)
ir.tl.chain_qc(adata)
adata = adata[~adata.obs["chain_pairing"].isin(["orphan VDJ", "orphan VJ", "extra VJ", "extra VDJ", "two full chains"]), :].copy()
ir.tl.clonal_expansion(adata)
ir.pl.clonal_expansion(adata, groupby="cluster", clip_at=4, normalize=False)

# FIGURE 5J

# Python in Spyder in Anaconda
ir.pl.vdj_usage(adata, full_combination=False, max_segments=None, max_ribbons=30, fig_kws={"figsize": (8, 5)},)

# FIGURE 5K

# Python in Spyder in Anaconda
sc.pl.umap(adata, color="cc_aa_alignment", groups=["52", "81", "103"], palette=cycler(color=mpl_cm.Dark2_r.colors),)
