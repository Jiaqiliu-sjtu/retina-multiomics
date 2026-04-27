
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
# V1(0421)：用对方的过滤参数
seu_filter <- subset(seu,
                     subset = nFeature_RNA > 300 &
                       nCount_RNA > 500 &
                       percent.mt < 2
)
cells_after <- ncol(seu_filter)
print(paste("过滤后细胞数量：", cells_after))
round(table(seu_filter$orig.ident)/table(seu$orig.ident)*100, 2)
seu=seu_filter
rm(seu_filter)
gc()
suppressWarnings({
  suppressMessages({
    
    library(Seurat)
    library(Signac)
    library(harmony)
    library(patchwork)
    
    ## =========================================================
    ## 0. 基本检查 / 准备分组变量
    ## =========================================================
    if (!"group" %in% colnames(seu@meta.data)) {
      if ("cond" %in% colnames(seu@meta.data)) {
        seu$group <- paste0("peri", seu$cond)
      } else if ("age" %in% colnames(seu@meta.data)) {
        seu$group <- paste0("peri", seu$age)
      } else {
        stop("seu@meta.data 里没有 group / cond / age，无法继续。")
      }
    }
    
    cat("Cells in seu:", ncol(seu), "\n")
    cat("Assays in seu:", Assays(seu), "\n")
    cat("Reductions in seu before rerun:\n")
    print(Reductions(seu))
    
    ## =========================================================
    ## 1. RNA 预处理 + Harmony
    ## =========================================================
    DefaultAssay(seu) <- "RNA"
    
    seu <- NormalizeData(seu)
    seu <- FindVariableFeatures(seu)
    seu <- ScaleData(seu)
    seu <- RunPCA(seu, npcs = 50)
    
    ## 如果你想保留旧结果，也可以把 reduction.save 改成别的名字
    seu <- RunHarmony(
      object = seu,
      group.by.vars = "cond",
      reduction.use = "pca",
      reduction.save = "integrated.harmony"
    )
    
    ## RNA 原始 UMAP（可选）
    seu <- RunUMAP(
      object = seu,
      reduction = "pca",
      dims = 1:10,
      reduction.name = "umap.rna.unintegrated",
      reduction.key = "rnaUMAP_"
    )
    
    ## RNA harmony UMAP
    seu <- RunUMAP(
      object = seu,
      reduction = "integrated.harmony",
      dims = 1:10,
      reduction.name = "umap.rna.harmony",
      reduction.key = "rnaHarmonyUMAP_"
    )
    
    ## =========================================================
    ## 2. ATAC 预处理 + Harmony(theta = 6)
    ## =========================================================
    DefaultAssay(seu) <- "ATAC"
    
    seu <- RunTFIDF(seu)
    seu <- FindTopFeatures(seu, min.cutoff = "q0")
    seu <- RunSVD(seu)
    
    seu <- RunHarmony(
      object = seu,
      group.by.vars = "cond",
      reduction.use = "lsi",
      assay.use = "ATAC",
      reduction.save = "harmony.atac",
      dim.use = 2:20,
      project.dim = FALSE,
      theta = 6,
      lambda = 1,
      sigma = 0.1,
      nclust = 50,
      max.iter.harmony = 20,
      max.iter.cluster = 100
    )
    
    ## raw ATAC UMAP
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
    
    ## harmony ATAC UMAP
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
    
    ## =========================================================
    ## 3. WNN
    ## =========================================================
    seu <- FindMultiModalNeighbors(
      object = seu,
      reduction.list = list("integrated.harmony", "harmony.atac"),
      dims.list = list(1:10, 2:20),
      k.nn = 20,
      l2.norm = TRUE,
      knn.graph.name = "wknn",
      snn.graph.name = "wsnn",
      weighted.nn.name = "weighted.nn",
      modality.weight.name = c("RNA.weight", "ATAC.weight"),
      knn.range = 200,
      prune.SNN = 0.06666667,
      sd.scale = 1,
      cross.contant.list = list(1e-4, 1e-4),
      smooth = FALSE,
      verbose = TRUE
    )
    
    ## WNN UMAP
    seu <- RunUMAP(
      object = seu,
      nn.name = "weighted.nn",
      reduction.name = "wnn.umap",
      reduction.key = "wnnUMAP_",
      umap.method = "uwot",
      metric = "cosine",
      seed.use = 42,
      verbose = TRUE
    )
    
    ## clustering
    seu <- FindClusters(
      object = seu,
      graph.name = "wsnn",
      resolution = 0.2,
      verbose = TRUE
    )
    
    ## =========================================================
    ## 4. 输出检查
    ## =========================================================
    cat("\n===== FINAL CHECK =====\n")
    cat("Reductions in seu after rerun:\n")
    print(Reductions(seu))
    
    cat("\nMetadata columns related to cluster/weights:\n")
    print(grep("cluster|res|weight", colnames(seu@meta.data), value = TRUE))
    
    cat("\nRNA.weight summary:\n")
    print(summary(seu$RNA.weight))
    
    cat("\nATAC.weight summary:\n")
    print(summary(seu$ATAC.weight))
    
  })
})
outdir <- "./results_final/reproduce_bgi"
qsave(seu,file.path(outdir,'0423_reproduce_bgi_harmonyadj_nocallpeak_noanno.qs'))
p1 <- DimPlot(seu, reduction = "umap.rna.harmony", group.by = "cond",label=F, cols = my36colors) + ggtitle("RNA")
p2 <- DimPlot(seu, reduction = "umap.atac.harmony", group.by = "cond",label=F, cols = my36colors) + ggtitle("ATAC")
p3 <- DimPlot(seu, reduction = "wnn.umap", group.by = "cond",label=F, cols = my36colors) + ggtitle("WNN")
ggsave(file.path(outdir,"0423_sample_rna_atac_wnn_wnnres0.2.pdf"),p1+p2+p3,width=22,height=8)


