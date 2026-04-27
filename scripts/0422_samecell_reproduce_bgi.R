# INTO R
setwd("./multiomics_macaque_4samples/")
library(Seurat)
library(Signac)
library(Matrix)
library(ggplot2)
library(readr)
library(RColorBrewer)
library(dplyr)
library(reshape2)
library(qs)
library(patchwork)
library(ComplexHeatmap)
library(harmony) 
library(SeuratWrappers)
# library(BSgenome.Mfascicularis.NCBI.6.0)
library(rtracklayer)
library(GenomicRanges)
my36colors <-c( '#E5D2DD', '#53A85F', '#F1BB72', '#F3B1A0', '#D6E7A3', '#57C3F3', '#476D87',
                '#E95C59', '#E59CC4', '#AB3282', '#23452F', '#BD956A', '#8C549C', '#585658',
                '#9FA3A8', '#E0D4CA', '#5F3D69', '#C5DEBA', '#58A4C3', '#E4C755', '#F7F398',
                '#AA9A59', '#E63863', '#E39A35', '#C1E6F3', '#6778AE', '#91D0BE', '#B53E2B',
                '#712820', '#DCC1DD', '#CCE0F5', '#CCC9E6', '#625D9E', '#68A180', '#3A6963',
                '#968175', "#6495ED", "#FFC1C1",'#f1ac9d','#f06966','#dee2d1','#6abe83','#39BAE8','#B9EDF8','#221a12',
                '#b8d00a','#74828F','#96C0CE','#E95D22','#017890')
marker_list <- list(
  RGC = c("NEFL", "SNCG", "SYT2", "SLC17A6"),
  AC = c("GAD1", "GAD2", "TFAP2A"),
  BC = c("TRPM1", "VSX2", "ISL1","GRM6"),
  HC = c("ONECUT1", "ONECUT2","SLC12A7","RET"),
  Cone = c("GNAT2", "GNGT2","ARR3", "PDE6C"),
  Rod = c("RHO", "PDE6A", "GNAT1", "NRL","SAG"),
  Muller = c("GLUL", "RLBP1", "SLC1A3"),
  RPE = c("RPE65", "BEST1", "TTR"),
  Microglia = c("C1QC", "C1QA", "CX3CR1","PTPRC")
)

#================================================================
# 读入自己的数据
#================================================================
seu=readRDS("./results_final/macaque_multiome_seu_raw_after_qc_before_filter.rds") 
# V2(0422)：选取对方保留的细胞
seu_bgi=readRDS("./results_by_bgi/2026-04-10_filtered_integrated_rmRPE_annotated.rds")
cells_use <- colnames(seu_bgi)
length(cells_use)
cells_raw <- colnames(seu)
cells_intersect <- intersect(cells_use, cells_raw)
length(cells_use)
length(cells_raw)
length(cells_intersect)
setdiff(cells_use, cells_raw) # 查看细胞名是否一致
seu <- subset(seu, cells = cells_intersect)
all(colnames(seu) %in% cells_use)
#================================================================
# 降维聚类
#================================================================
## =========================================
## Reproduce donor workflow as closely as possible
## based on extracted seu@commands
## =========================================

set.seed(42)

## -----------------------------------------
## 1. RNA preprocessing
## fully matched
## -----------------------------------------
suppressWarnings({
  suppressMessages({
    DefaultAssay(seu) <- "RNA"
    
    seu <- NormalizeData(
      object = seu,
      assay = "RNA",
      verbose = TRUE
    )
    
    seu <- FindVariableFeatures(
      object = seu,
      assay = "RNA",
      selection.method = "vst",
      nfeatures = 2000,
      verbose = TRUE
    )
    
    seu <- ScaleData(
      object = seu,
      assay = "RNA",
      verbose = TRUE
    )
    
    seu <- RunPCA(
      object = seu,
      assay = "RNA",
      npcs = 50,
      rev.pca = FALSE,
      weight.by.var = TRUE,
      verbose = TRUE,
      reduction.name = "pca",
      reduction.key = "PC_",
      seed.use = 42
    )
  })
})

## -----------------------------------------
## 2. RNA unintegrated neighbors
## fully matched
## -----------------------------------------
suppressWarnings({
  suppressMessages({
    seu <- FindNeighbors(
      object = seu,
      reduction = "pca",
      dims = 1:10,
      assay = "RNA",
      k.param = 20,
      return.neighbor = FALSE,
      compute.SNN = TRUE,
      prune.SNN = 0.06666667,
      nn.method = "annoy",
      n.trees = 50,
      annoy.metric = "euclidean",
      nn.eps = 0,
      verbose = TRUE,
      graph.name = c("RNA_nn", "RNA_snn"),
      l2.norm = FALSE,
      cache.index = FALSE
    )
  })
})

## -----------------------------------------
## 3. RNA unintegrated UMAP
## fully matched
## -----------------------------------------
suppressWarnings({
  suppressMessages({
    seu <- RunUMAP(
      object = seu,
      reduction = "pca",
      dims = 1:10,
      assay = "RNA",
      slot = "data",
      umap.method = "uwot",
      return.model = FALSE,
      n.neighbors = 30,
      n.components = 2,
      metric = "cosine",
      learning.rate = 1,
      min.dist = 0.3,
      spread = 1,
      set.op.mix.ratio = 1,
      local.connectivity = 1,
      repulsion.strength = 1,
      negative.sample.rate = 5,
      uwot.sgd = FALSE,
      seed.use = 42,
      angular.rp.forest = FALSE,
      densmap = FALSE,
      dens.lambda = 2,
      dens.frac = 0.3,
      dens.var.shift = 0.1,
      verbose = TRUE,
      reduction.name = "umap.rna.unintegrated"
    )
  })
})

## -----------------------------------------
## 4. RNA harmony integration
## command not available in object, use best-known reconstruction
## -----------------------------------------
suppressWarnings({
  suppressMessages({
    seu <- IntegrateLayers(
      object = seu,
      assay = "RNA",
      method = HarmonyIntegration,
      orig.reduction = "pca",
      group.by = "cond",
      new.reduction = "integrated.harmony",
      verbose = FALSE
    )
  })
})

## -----------------------------------------
## 5. RNA harmony neighbors
## fully matched
## -----------------------------------------
suppressWarnings({
  suppressMessages({
    seu <- FindNeighbors(
      object = seu,
      reduction = "integrated.harmony",
      dims = 1:10,
      assay = "RNA",
      k.param = 20,
      return.neighbor = FALSE,
      compute.SNN = TRUE,
      prune.SNN = 0.06666667,
      nn.method = "annoy",
      n.trees = 50,
      annoy.metric = "euclidean",
      nn.eps = 0,
      verbose = TRUE,
      graph.name = c("RNA_nn", "RNA_snn"),
      l2.norm = FALSE,
      cache.index = FALSE
    )
  })
})

## -----------------------------------------
## 6. RNA harmony UMAP
## fully matched
## -----------------------------------------
suppressWarnings({
  suppressMessages({
    seu <- RunUMAP(
      object = seu,
      reduction = "integrated.harmony",
      dims = 1:10,
      assay = "RNA",
      slot = "data",
      umap.method = "uwot",
      return.model = FALSE,
      n.neighbors = 30,
      n.components = 2,
      metric = "cosine",
      learning.rate = 1,
      min.dist = 0.3,
      spread = 1,
      set.op.mix.ratio = 1,
      local.connectivity = 1,
      repulsion.strength = 1,
      negative.sample.rate = 5,
      uwot.sgd = FALSE,
      seed.use = 42,
      angular.rp.forest = FALSE,
      densmap = FALSE,
      dens.lambda = 2,
      dens.frac = 0.3,
      dens.var.shift = 0.1,
      verbose = TRUE,
      reduction.name = "umap.rna.harmony",
      reduction.key = "rnaHarmonyUMAP_"
    )
  })
})

## -----------------------------------------
## 7. ATAC preprocessing
## RunHarmony.ATAC command missing in object
## use best-known reconstruction based on downstream reductions
## -----------------------------------------
suppressWarnings({
  suppressMessages({
    DefaultAssay(seu) <- "ATAC"
    
    seu <- RunTFIDF(
      object = seu,
      assay = "ATAC"
    )
    
    seu <- FindTopFeatures(
      object = seu,
      assay = "ATAC",
      min.cutoff = "q0"
    )
    
    seu <- RunSVD(
      object = seu,
      assay = "ATAC"
    )
    
    seu <- RunHarmony(
      object = seu,
      group.by.vars = "cond",
      reduction.use = "lsi",
      assay.use = "ATAC",
      reduction.save = "harmony.atac",
      dim.use = 2:20,
      project.dim = FALSE
    )
  })
})

## -----------------------------------------
## 8. ATAC raw LSI UMAP
## fully matched
## -----------------------------------------
suppressWarnings({
  suppressMessages({
    seu <- RunUMAP(
      object = seu,
      reduction = "lsi",
      dims = 2:20,
      assay = "ATAC",
      slot = "data",
      umap.method = "uwot",
      return.model = FALSE,
      n.neighbors = 30,
      n.components = 2,
      metric = "cosine",
      learning.rate = 1,
      min.dist = 0.3,
      spread = 1,
      set.op.mix.ratio = 1,
      local.connectivity = 1,
      repulsion.strength = 1,
      negative.sample.rate = 5,
      uwot.sgd = FALSE,
      seed.use = 42,
      angular.rp.forest = FALSE,
      densmap = FALSE,
      dens.lambda = 2,
      dens.frac = 0.3,
      dens.var.shift = 0.1,
      verbose = TRUE,
      reduction.name = "umap.atac",
      reduction.key = "atacUMAP_"
    )
  })
})

## -----------------------------------------
## 9. ATAC harmony UMAP
## fully matched
## -----------------------------------------
suppressWarnings({
  suppressMessages({
    seu <- RunUMAP(
      object = seu,
      reduction = "harmony.atac",
      dims = 2:20,
      assay = "ATAC",
      slot = "data",
      umap.method = "uwot",
      return.model = FALSE,
      n.neighbors = 30,
      n.components = 2,
      metric = "cosine",
      learning.rate = 1,
      min.dist = 0.3,
      spread = 1,
      set.op.mix.ratio = 1,
      local.connectivity = 1,
      repulsion.strength = 1,
      negative.sample.rate = 5,
      uwot.sgd = FALSE,
      seed.use = 42,
      angular.rp.forest = FALSE,
      densmap = FALSE,
      dens.lambda = 2,
      dens.frac = 0.3,
      dens.var.shift = 0.1,
      verbose = TRUE,
      reduction.name = "umap.atac.harmony",
      reduction.key = "atacHarmonyUMAP_"
    )
  })
})

## -----------------------------------------
## 10. WNN neighbors
## fully matched at level recoverable from object
## -----------------------------------------
suppressWarnings({
  suppressMessages({
    seu <- FindMultiModalNeighbors(
      object = seu,
      reduction.list = list("integrated.harmony", "harmony.atac"),
      dims.list = list(1:10, 2:20)
    )
  })
})

## -----------------------------------------
## 11. WNN clustering
## fully matched
## -----------------------------------------
suppressWarnings({
  suppressMessages({
    seu <- FindClusters(
      object = seu,
      graph.name = "wsnn",
      cluster.name = "wsnn_res.0.2",
      modularity.fxn = 1,
      resolution = 0.2,
      method = "matrix",
      algorithm = 3,
      n.start = 10,
      n.iter = 10,
      random.seed = 0,
      group.singletons = TRUE,
      verbose = FALSE
    )
  })
})

## -----------------------------------------
## 12. WNN UMAP
## command missing in object, use standard Seurat call
## -----------------------------------------
suppressWarnings({
  suppressMessages({
    seu <- RunUMAP(
      object = seu,
      nn.name = "weighted.nn",
      reduction.name = "wnn.umap",
      reduction.key = "wnnUMAP_",
      seed.use = 42,
      verbose = TRUE
    )
  })
})

## -----------------------------------------
## 13. final cluster column
## -----------------------------------------
if ("wsnn_res.0.2" %in% colnames(seu@meta.data)) {
  seu$seurat_clusters <- seu$wsnn_res.0.2
}

#================================================================
# 结果展示
#================================================================
outdir <- './results_final/reproduce_bgi/'
dir.create(outdir, showWarnings = FALSE)
qsave(seu,file.path(outdir,"0422_reproduce_bgi_same_cell_nocallpeak_noanno.qs"))
'''
p1 <- DimPlot(seu, reduction = "umap.rna.harmony", group.by = "wsnn_res.0.2",label=T, cols = my36colors) + ggtitle("RNA")
p2 <- DimPlot(seu, reduction = "umap.atac.harmony", group.by = "wsnn_res.0.2",label=T, cols = my36colors) + ggtitle("ATAC")
p3 <- DimPlot(seu, reduction = "wnn.umap", group.by = "wsnn_res.0.2",label=T, cols = my36colors) + ggtitle("WNN")
ggsave(file.path(outdir,"0421_cluster_rna_atac_wnn_wnnres0.2.pdf"),p1+p2+p3,width=22,height=8)

final_resolution <- 0.2
seu$seurat_clusters <- seu[[paste0("wsnn_res.", final_resolution)]][, 1]
print(table(seu$seurat_clusters))

p1 <- DimPlot(seu, reduction = "umap.rna.harmony", group.by = "cond",label=F, cols = my36colors) + ggtitle("RNA")
p2 <- DimPlot(seu, reduction = "umap.atac.harmony", group.by = "cond",label=F, cols = my36colors) + ggtitle("ATAC")
p3 <- DimPlot(seu, reduction = "wnn.umap", group.by = "cond",label=F, cols = my36colors) + ggtitle("WNN")
ggsave(file.path(outdir,"0421_sample_rna_atac_wnn_wnnres0.2.pdf"),p1+p2+p3,width=22,height=8)

DefaultAssay(seu)="RNA"
p <- DotPlot(seu,assay = "RNA", 
             group.by = "wsnn_res.0.2", 
             features = marker_list,
             cols = c("grey","red"),
             dot.scale = 5)+  
  RotatedAxis() + 
  scale_x_discrete("") + 
  scale_y_discrete("") +
  theme(
    axis.text.x = element_text(size = 10,  
                               angle = 45, hjust = 1, vjust = 1)
  ) +
  ggtitle("Marker Genes Expression wnn_res.0.2") +
  labs(color = "Expression\nLevel")  # 修改图例标题
print(p)
ggsave(file.path(outdir,"0421_marker_dotplot_wnnres0.2.pdf"), p,width = 22, height = 8)

# res0.2
seu$celltype <- recode(seu$seurat_clusters,
                       "0"  = "Rod",
                       "1"  = "Rod",
                       "2"  = "Rod",
                       "3"  = "Muller",
                       "4"  = "Cone",
                       "5"  = "AC", # ?
                       "6"  = "RGC",
                       "7"  = "AC",
                       "8"  = "BC",
                       "9"  = "Muller", # ?
                       "10" = "BC",
                       "11" = "BC",
                       "12" = "BC",
                       "13" = "HC",
                       "14" = "BC",
                       "15" = "Muller", #?
                       "16" = "BC",
                       "17" = "BC", 
                       "18" = "AC",
                       "19" = "Microglia",
                       "20" = "BC", 
                       "21" = "RPE",
                       "22" = "Rod",
                       "23" = "BC",
                       "24" = "HC"
)
# 注释好的细胞类型dotplot
p <- DotPlot(seu,assay = "RNA", 
             group.by = "celltype", 
             features = marker_list,
             cols = c("grey","red"),
             dot.scale = 5)+  
  RotatedAxis() + 
  scale_x_discrete("") + 
  scale_y_discrete("") +
  theme(
    axis.text.x = element_text(size = 10,  
                               angle = 45, hjust = 1, vjust = 1)
  ) +
  ggtitle("Marker Genes Expression") +
  labs(color = "Expression\nLevel")  # 修改图例标题
print(p)
ggsave("./results_final//0_QC/filter_by_bgi_10_20_dotplot_celltype_nocallpeak_wnn_res0.2.pdf", width = 22, height = 8)

# wnn图不同res展示细胞类型
p <- DimPlot(seu, cols = my36colors, label = T,repel=T,
             pt.size = 0.3,reduction = "wnn.umap", group.by = "celltype" )+
  ggtitle(paste0("celltype wnn_res",final_resolution))
ggsave(file.path("./results_final/0_QC/",paste0("filter_by_bgi_10_20_UMAP_multiomics_by_celltype_wnn_res",final_resolution,".pdf")), p,width=8, height = 6)

table(seu$celltype)
# rna_atac_wnn 分样本展示细胞类型
p1 <- DimPlot(seu, reduction = "rnaumap", group.by = "celltype",split.by = "orig.ident",label=T, cols = my36colors) + ggtitle("RNA")
p2 <- DimPlot(seu, reduction = "atacumap", group.by = "celltype",split.by = "orig.ident",label=T, cols = my36colors) + ggtitle("ATAC")
p3 <- DimPlot(seu, reduction = "wnn.umap", group.by = "celltype",split.by = "orig.ident",label=T, cols = my36colors) + ggtitle("WNN")
p_combined <- (p1 + p2 + p3) +
  plot_annotation(title = "res0.2")
print(p_combined)
ggsave(
  "./results_final/0_QC/filter_by_bgi_10_20_UMAP_multiomics_celltype_res0.2.pdf",
  p_combined,
  width = 20,
  height = 15
)
p1 <- DimPlot(seu, reduction = "rnaumap", group.by = "orig.ident",label=F) + ggtitle("RNA")
p2 <- DimPlot(seu, reduction = "atacumap", group.by = "orig.ident",label=F) + ggtitle("ATAC")
p3 <- DimPlot(seu, reduction = "wnn.umap", group.by = "orig.ident",label=F) + ggtitle("WNN")

ggsave("results_final/0_QC/asdf by sample.pdf",p1+p2+p3, width = 20, height = 6)

p1 <- DimPlot(seu, reduction = "rnaumap", group.by = "celltype",label=T) + ggtitle("RNA")
p2 <- DimPlot(seu, reduction = "atacumap", group.by = "celltype",label=T) + ggtitle("ATAC")
p3 <- DimPlot(seu, reduction = "wnn.umap", group.by = "celltype",label=T) + ggtitle("WNN")

ggsave("results_final/0_QC/asdf_celltype.pdf",p1+p2+p3, width = 20, height = 6)
'''
# qsave(seu,file.path(outdir,'0421_reproduce_bgi_nocallpeak_noanno.qs'))
