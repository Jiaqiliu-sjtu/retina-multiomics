#================================================================
# 读入callpeak后的对象
seu=qread("./results_final/multiomics_after_callpeak_no_further.qs")
#================================================================
# 质控
suppressWarnings({
  suppressMessages({
    DefaultAssay(seu) <- "RNA"
    seu[["percent.mt"]] <- PercentageFeatureSet(seu, pattern = "^MT")
    DefaultAssay(seu) <- "ATAC"
    seu <- TSSEnrichment(object = seu, fast = FALSE) 
    seu <- NucleosomeSignal(object = seu)
  })
})
# 添加分组变量cond并按年龄排序
seu$cond <-paste0(gsub("yperi([0-9]+).*", "\\1", seu$orig.ident), "y")
seu$cond <- factor(seu$cond,levels=c("5y","15y","20y","29y"))
table(seu$cond)
#================================================================
# 读入scATAC对象，构建seu
scATAC=qread("./results_final/macaretina_scATAC.qs")
seu$nCount_ATAC <- scATAC$nCount_peaks
seu$nFeature_ATAC <- scATAC$nFeature_peaks
summary(seu@meta.data[, c(
  "nFeature_RNA", "nCount_RNA", "percent.mt",
  "nCount_ATAC", "TSS.enrichment", "nucleosome_signal"
)])
# 保存含有质控指标的seu
qsave(seu,"./results_final/multiomics_after_callpeak_qc_0423.qs")

cells_before <- ncol(seu)
print(paste("过滤前细胞数量：", cells_before))
#================================================================
# 设置QC过滤指标
seu_filter <- subset(seu,
                     subset = nFeature_RNA > 100 &
                       nFeature_RNA < 7000 &
                       nCount_RNA > 500 &
                       nCount_RNA < 30000 &
                       percent.mt < 20 &
                       nCount_ATAC > 100 &
                       nCount_ATAC < 100000 &
                       TSS.enrichment > 1 &
                       nucleosome_signal < 2
)
cells_after <- ncol(seu_filter)
print(paste("过滤后细胞数量：", cells_after))
print(paste("过滤掉的细胞数量：", cells_before - cells_after))
print(paste("保留比例：", round(cells_after/cells_before*100, 2), "%"))
seu=seu_filter
rm(seu_filter)
gc()
# 保存质控过滤后的seu
qsave(seu, "results_final/multiomics_after_callpeak_qc_filter_before_reduction_0423.qs")
#================================================================
# RNA降维聚类、harmony
#================================================================
dimrna <- 10
resolution_grid <- c(0.02,0.1,0.2,0.3,0.5,0.8)
suppressWarnings({
  suppressMessages({
    DefaultAssay(seu) <- "RNA"
    seu <- NormalizeData(seu, assay = "RNA")
    seu <- FindVariableFeatures(seu, assay = "RNA", selection.method = "vst", nfeatures = 2000)
    seu <- ScaleData(seu, assay = "RNA")
    seu <- RunPCA(seu, assay = "RNA", npcs = 50)
    seu <- IntegrateLayers(
      object = seu,
      assay = "RNA",
      method = HarmonyIntegration, 
      orig.reduction = "pca",group.by="cond",
      new.reduction = "integrated.rna",
      verbose = FALSE
    )
  })
})
# RNA harmony后降维聚类
suppressWarnings({
  suppressMessages({
    seu <- RunUMAP(seu, reduction = "integrated.rna", dims = 1:dimrna, assay = "RNA",reduction.name="rna.umap")
    seu <- FindNeighbors(seu, reduction = "integrated.rna", dims = 1:dimrna, assay = "RNA",graph.name = c("rna_nn", "rna_snn"))
    seu <- FindClusters(seu, resolution = resolution_grid, algorithm = 1,graph.name = "rna_snn")
  })
})
#================================================================
# ATAC降维聚类、harmony
#================================================================
dimatac <- 20
suppressWarnings({
  suppressMessages({
    seu <- RunTFIDF(seu, assay = "ATAC")
    seu <- FindTopFeatures(seu, assay = "ATAC", min.cutoff = 'q0')
    seu <- RunSVD(seu, assay = "ATAC")
    seu <- RunHarmony(
      object = seu,
      group.by.vars = "cond",
      reduction.use = "lsi",
      assay.use = "ATAC",
      reduction.save = "integrated.atac",
      dim.use = 2:dimatac,
      project.dim = FALSE,
      theta = 6,
      lambda = 1,
      sigma = 0.1,
      nclust = 50,
      max.iter.harmony = 20,
      max.iter.cluster = 100
    )
  })
})

### skip 1 dim or not
outdir <- "results_final/0423_callpeak_harmony_adjust/"
p <- DepthCor(seu, assay = "ATAC", reduction = "integrated.atac")
ggsave(file.path(outdir, "DepthCor_integrated_atac.pdf"),p, width = 6, height = 4)

#================================================================
# WNN
#================================================================
suppressWarnings({
  suppressMessages({
    seu <- RunUMAP(seu, reduction = "integrated.atac", dims = 2:dimatac, assay = "ATAC",reduction.name="atac.umap")
    seu <- FindNeighbors(seu, reduction = "integrated.atac", dims = 2:dimatac, assay = "ATAC",graph.name =c("atac_nn", "atac_snn"))
    seu <- FindClusters(seu, resolution = resolution_grid, algorithm = 1,graph.name = "atac_snn")
  })
})
suppressWarnings({
  suppressMessages({
    seu <- FindMultiModalNeighbors(seu, reduction.list = list("integrated.rna", "integrated.atac"), dims.list = list(1:dimrna, 2:dimatac)) # 会生成weighted.nn和wsnn
    # 基于WNN进行聚类
    seu <- FindClusters(seu, graph.name = "wsnn", resolution =resolution_grid)
    # WNN UMAP
    seu <- RunUMAP(seu, nn.name = "weighted.nn", reduction.name = "wnn.umap") # weighted.nn不是一个 PCA/LSI 坐标矩阵，而是rna+atac综合判断每个细胞的邻居关系
  })
})
# 保存跑完流程但是没有选res、没有注释细胞的版本
qsave(seu, "results_final/macaretina_multiomics_callpeak_d10_d20_after_reduction_harmonytheta6_0423.qs")  
#====================================
# 看样本在各模态的整合情况
p1 <- DimPlot(seu, reduction = "rna.umap", group.by = "orig.ident",label=F) + ggtitle("RNA")
p2 <- DimPlot(seu, reduction = "atac.umap", group.by = "orig.ident",label=F) + ggtitle("ATAC")
p3 <- DimPlot(seu, reduction = "wnn.umap", group.by = "orig.ident",label=F) + ggtitle("WNN")
ggsave("results_final/0423_callpeak_harmony_adjust/0423_sample_umap_10_20_harmont_theta6_res0.2_multiomics.pdf",p1+p2+p3, width = 20, height = 6)
#====================================
# 看不同res下的cluster
p1 <- DimPlot(seu, reduction = "rna.umap", group.by = "wsnn_res.0.2",label=T, cols = my36colors) + ggtitle("RNA")
p2 <- DimPlot(seu, reduction = "atac.umap", group.by = "wsnn_res.0.2",label=T, cols = my36colors) + ggtitle("ATAC")
p3 <- DimPlot(seu, reduction = "wnn.umap", group.by = "wsnn_res.0.2",label=T, cols = my36colors) + ggtitle("WNN")
p_combined <- (p1 + p2 + p3) +plot_annotation(title = "WNN_0.2")
ggsave("results_final/0423_callpeak_harmony_adjust/0423_cluster_umap_10_20_harmony_theta6res0.2_multiomics.pdf",p_combined,width = 20,height = 6)
#====================================
# 看不同res下的dotplot
DefaultAssay(seu)="RNA"
p <- DotPlot(seu,assay = "RNA", 
             group.by = "wsnn_res.0.2", 
             features = marker_list,
             cols = c("grey","red"),
             dot.scale = 5)+
  RotatedAxis() + 
  scale_x_discrete("") + 
  scale_y_discrete("") +
  theme(axis.text.x = element_text(size = 10, angle = 45, hjust = 1, vjust = 1)) +
  ggtitle("Marker Genes Expression wnn_res.0.2") + # 修改图例标题
  labs(color = "Expression\nLevel")  
ggsave("./results_final/0423_callpeak_harmony_adjust/0423_marker_dotplot_d10_d20_harmony_theta6_wnn_res0.2_multiomics.pdf", p,width = 22, height = 8)
#====================================
# 细胞注释 
final_resolution <- 0.2
seu$seurat_clusters <- seu[[paste0("wsnn_res.", final_resolution)]][, 1]
print(table(seu$seurat_clusters))
# res0.2
seu$celltype <- recode(seu$seurat_clusters,
                       "0"  = "Rod",
                       "1"  = "Muller",
                       "2"  = "BC",
                       "3"  = "Cone",
                       "4"  = "BC",
                       "5"  = "Rod",
                       "6"  = "BC",
                       "7"  = "RGC",
                       "8"  = "Muller",
                       "9"  = "BC",
                       "10" = "HC",
                       "11" = "AC",
                       "12" = "AC", #?
                       "13" = "AC",
                       "14" = "BC",
                       "15" = "BC", 
                       "16" = "Muller", # ?
                       "17" = "Muller", 
                       "18" = "AC",
                       "19" = "HC",
                       "20" = "BC", 
                       "21" = "Microglia",
                       "22" = "RPE",
                       "23" = "Rod"
)
# 生成每个cluster的信息
seu$cluster.subtype <- paste0(seu$seurat_clusters, ". ", seu$celltype)
# 注释好的细胞类型dotplot
p <- DotPlot(seu,assay = "RNA", 
             group.by = "celltype", 
             features = marker_list,
             cols = c("grey","red"),
             dot.scale = 5)+  
  RotatedAxis() + 
  scale_x_discrete("") + 
  scale_y_discrete("") +
  theme(axis.text.x = element_text(size = 10,  
                               angle = 45, hjust = 1, vjust = 1)) +
  ggtitle("Marker Genes Expression res 0.2") +
  labs(color = "Expression\nLevel")  # 修改图例标题
ggsave("./results_final/0423_callpeak_harmony_adjust/0423_dotplot_celltype_aftercallpeak_wnn_res0.2.pdf", width = 22, height = 8)
#====================================
# 展示wnn_细胞类型
p <- DimPlot(seu, cols = my36colors, label = T,repel=T,
             pt.size = 0.3,reduction = "wnn.umap", group.by = "celltype" )+
  ggtitle(paste0("celltype wnn_res",final_resolution))
ggsave(file.path(outdir,paste0("0423_anno_umap_wnn_res",final_resolution,".pdf")), p,width=8, height = 6)
# 所有样本细胞类型的比例
table(seu$celltype)
#====================================
# rna_atac_wnn celltype
p1 <- DimPlot(seu, reduction = "rna.umap", group.by = "celltype",label=T) + ggtitle("RNA")
p2 <- DimPlot(seu, reduction = "atac.umap", group.by = "celltype",label=T) + ggtitle("ATAC")
p3 <- DimPlot(seu, reduction = "wnn.umap", group.by = "celltype",label=T) + ggtitle("WNN")
ggsave("results_final/0423_callpeak_harmony_adjust/0423_anno_umap_res0.2.pdf",p1+p2+p3, width = 20, height = 6)
# rna_atac_wnn sample
p1 <- DimPlot(seu, reduction = "rna.umap", group.by = "orig.ident",label=F) + ggtitle("RNA")
p2 <- DimPlot(seu, reduction = "atac.umap", group.by = "orig.ident",label=F) + ggtitle("ATAC")
p3 <- DimPlot(seu, reduction = "wnn.umap", group.by = "orig.ident",label=F) + ggtitle("WNN")
ggsave("results_final/0423_callpeak_harmony_adjust/0423_anno_sample_res0.2.pdf",p1+p2+p3, width = 20, height = 6)

# rna_atac_wnn 分样本展示细胞类型（3 x 3 图！）
# 感觉画wnn的就可以
# p1 <- DimPlot(seu, reduction = "rna.umap", group.by = "celltype",split.by = "orig.ident",label=T, cols = my36colors) + ggtitle("RNA")
# p2 <- DimPlot(seu, reduction = "atac.umap", group.by = "celltype",split.by = "orig.ident",label=T, cols = my36colors) + ggtitle("ATAC")
p3 <- DimPlot(seu, reduction = "wnn.umap", group.by = "celltype",split.by = "orig.ident",label=T, cols = my36colors) + ggtitle("WNN")
p_combined <- (p1 + p2 + p3) +plot_annotation(title = "res0.2")
ggsave("./results_final/0423_callpeak_harmony_adjust/0423_sample_anno_umap_res0.2.pdf",p_combined,width = 20,height = 15)
#====================================

#====================================
# 保存完成所有流程的版本
qsave(seu, "results_final/macaretina_multiomics_callpeak_d10_d20_after_reduction_harmonytheta6_anno_0423.qs")  
