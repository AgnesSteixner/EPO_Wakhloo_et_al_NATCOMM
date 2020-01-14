#Analyses of scRNA matirx from EPO project using cross correlation analysis (SEURAT)
#9 May 2018
library(Seurat)
library(dplyr)
library(Matrix)
set.seed(786)


dump_all<-read.csv(file="../DATA/merged.csv")
rownames(dump_all)<-dump_all$GENE
dump_all$GENE<-NULL
dump_all[is.na(dump_all)]<-0



#CreateandInitializeSeuratObject
mm_data_all<-CreateSeuratObject(raw.data=dump_all,min.cells=3,min.genes=1000,project="AllGroups")

#CheckforMitochondrialgenesandcalculatepercentagesacrossthecells
mito.genes_all<-grep(pattern="^mt-",x=rownames(x=mm_data_all@data),value=TRUE)
percent.mito_all<-colSums(mm_data_all@raw.data[mito.genes_all,])/colSums(mm_data_all@raw.data)

#AddMetaDataandplot
mm_data_all<-AddMetaData(object=mm_data_all,metadata=percent.mito_all,col.name="percent.mito")


#FilterCellsbyremovingoutliersandmitochondodrialcellsandthenplot
mm_data_all<-FilterCells(object=mm_data_all,subset.names=c("nGene","percent.mito"),
low.thresholds=c(1000,-Inf),high.thresholds=c(8000,0.40))



dump_allV2 <- (dump_all[,(mm_data_all@cell.names)])

Group1 <- (dump_allV2[,(grep("gr1",colnames(dump_allV2)))])
Group2 <- (dump_allV2[,(grep("gr2",colnames(dump_allV2)))])

# Set up control object
ctrl <- CreateSeuratObject(raw.data = Group1, project = "MPIEM_CTRL", min.cells = 3,min.genes=1000)
ctrl@meta.data$stim <- "CTRL"
ctrl <- FilterCells(ctrl, subset.names = "nGene", low.thresholds = 1000, high.thresholds = 8000)
ctrl <- NormalizeData(ctrl)
mito.genes_all<-grep(pattern="^mt-",x=rownames(x=ctrl@data),value=TRUE)
percent.mito_all<-colSums(ctrl@raw.data[mito.genes_all,])/colSums(ctrl@raw.data)
ctrl<-AddMetaData(object=ctrl,metadata=percent.mito_all,col.name="percent.mito")

ctrl <- ScaleData(ctrl, display.progress = F,vars.to.regress=c("nUMI","percent.mito"))
# Set up stimulated object
stim <- CreateSeuratObject(raw.data = Group2, project = "MPIEM_STIM", min.cells = 3,min.genes=1000)
stim@meta.data$stim <- "STIM"
stim <- FilterCells(stim, subset.names = "nGene", low.thresholds = 1000, high.thresholds = 8000)
stim <- NormalizeData(stim)
mito.genes_all<-grep(pattern="^mt-",x=rownames(x=stim@data),value=TRUE)
percent.mito_all<-colSums(stim@raw.data[mito.genes_all,])/colSums(stim@raw.data)
stim<-AddMetaData(object=stim,metadata=percent.mito_all,col.name="percent.mito")

stim <- ScaleData(stim, display.progress = F,vars.to.regress=c("nUMI","percent.mito"))




# Gene selection for input to CCA
ctrl <- FindVariableGenes(ctrl, do.plot = F)
stim <- FindVariableGenes(stim, do.plot = F)
g.1 <- head(rownames(ctrl@hvg.info), 1000)
g.2 <- head(rownames(stim@hvg.info), 1000)
genes.use <- unique(c(g.1, g.2))
genes.use <- intersect(genes.use, rownames(ctrl@scale.data))
genes.use <- intersect(genes.use, rownames(stim@scale.data))



immune.combined <- RunCCA(ctrl, stim, genes.use = genes.use, num.cc = 30)


immune.combined <- AlignSubspace(immune.combined, reduction.type = "cca", grouping.var = "stim", 
    dims.align = 1:20)



# t-SNE and Clustering
immune.combined <- RunTSNE(immune.combined, reduction.use = "cca.aligned", dims.use = 1:20, 
    do.fast = T)
immune.combined <- FindClusters(immune.combined, reduction.type = "cca.aligned", 
    resolution = 1.1, dims.use = 1:20)


immune.tmp <- immune.combined

Ident_names2 <- immune.combined@ident


saveRDS(immune.combined, file = paste0(outDir,"immune_combined.rds"))



immune.combinedNamed <- immune.combined
new.ident <- c("Glutamatergic1","Glutamatergic1","Oligodendrocytes","Glutamatergic2","Astrocytes","Glutamatergic3","Endothelial","Gabaergic","OPCs","Glutamatergic4","Microglia","Endothelial" )
for (i in 0:11) {
    immune.combinedNamed <- RenameIdent(object = immune.combinedNamed, old.ident.name = i, 
        new.ident.name = new.ident[i + 1])
}

# find markers for every cluster compared to all remaining cells, report
# only the positive ones
All_markers <- FindAllMarkers(object = immune.combinedNamed, only.pos = TRUE, min.pct = 0.25, 
    thresh.use = 0.25)
All_markers %>% group_by(cluster) %>% top_n(2, avg_logFC)


All_markers_sig <- All_markers[(All_markers$p_val_adj<0.01),]



sessionInfo()


R version 3.4.1 (2017-06-30)
Platform: x86_64-pc-linux-gnu (64-bit)
Running under: Ubuntu 16.04.4 LTS

Matrix products: default
BLAS: /opt/conda/lib/R/lib/libRblas.so
LAPACK: /opt/conda/lib/R/lib/libRlapack.so

locale:
 [1] LC_CTYPE=en_US.UTF-8       LC_NUMERIC=C              
 [3] LC_TIME=en_US.UTF-8        LC_COLLATE=en_US.UTF-8    
 [5] LC_MONETARY=en_US.UTF-8    LC_MESSAGES=en_US.UTF-8   
 [7] LC_PAPER=en_US.UTF-8       LC_NAME=C                 
 [9] LC_ADDRESS=C               LC_TELEPHONE=C            
[11] LC_MEASUREMENT=en_US.UTF-8 LC_IDENTIFICATION=C       

attached base packages:
[1] stats     graphics  grDevices utils     datasets  methods   base     

other attached packages:
[1] bindrcpp_0.2.2 dplyr_0.7.6    Seurat_2.3.0   Matrix_1.2-12  cowplot_0.9.3 
[6] ggplot2_3.1.0 

loaded via a namespace (and not attached):
  [1] uuid_0.1-2           snow_0.4-2           backports_1.1.2     
  [4] Hmisc_4.1-1          VGAM_1.0-5           sn_1.5-2            
  [7] plyr_1.8.4           igraph_1.2.1         repr_0.13           
 [10] lazyeval_0.2.1       splines_3.4.1        digest_0.6.18       
 [13] foreach_1.4.4        htmltools_0.3.6      lars_1.2            
 [16] gdata_2.18.0         magrittr_1.5         checkmate_1.8.5     
 [19] cluster_2.0.7-1      mixtools_1.1.0       ROCR_1.0-7          
 [22] sfsmisc_1.1-1        recipes_0.1.2        gower_0.1.2         
 [25] dimRed_0.1.0         R.utils_2.6.0        colorspace_1.3-2    
 [28] crayon_1.3.4         jsonlite_1.5         bindr_0.1.1         
 [31] survival_2.40-1      zoo_1.8-1            iterators_1.0.10    
 [34] ape_5.1              glue_1.3.0           DRR_0.0.3           
 [37] gtable_0.2.0         ipred_0.9-6          kernlab_0.9-27      
 [40] ddalpha_1.3.2        prabclus_2.2-6       DEoptimR_1.0-8      
 [43] scales_1.0.0         mvtnorm_1.0-7        Rcpp_1.0.0          
 [46] metap_0.9            dtw_1.18-1           htmlTable_1.11.2    
 [49] magic_1.5-6          tclust_1.3-1         foreign_0.8-67      
 [52] proxy_0.4-22         mclust_5.4           SDMTools_1.1-221    
 [55] Formula_1.2-3        stats4_3.4.1         tsne_0.1-3          
 [58] lava_1.6.1           prodlim_2018.04.18   htmlwidgets_1.0     
 [61] FNN_1.1              gplots_3.0.1         RColorBrewer_1.1-2  
 [64] fpc_2.1-11           acepack_1.4.1        modeltools_0.2-21   
 [67] ica_1.0-1            pkgconfig_2.0.1      R.methodsS3_1.7.1   
 [70] flexmix_2.3-14       nnet_7.3-12          caret_6.0-79        
 [73] labeling_0.3         tidyselect_0.2.4     rlang_0.3.0.1       
 [76] reshape2_1.4.3       munsell_0.5.0        tools_3.4.1         
 [79] ranger_0.9.0         broom_0.4.4          ggridges_0.5.0      
 [82] evaluate_0.10.1      geometry_0.3-6       stringr_1.3.1       
 [85] ModelMetrics_1.1.0   knitr_1.20           fitdistrplus_1.0-9  
 [88] robustbase_0.93-0    caTools_1.17.1       purrr_0.2.4         
 [91] RANN_2.5.1           pbapply_1.3-4        nlme_3.1-137        
 [94] R.oo_1.22.0          RcppRoll_0.2.2       compiler_3.4.1      
 [97] rstudioapi_0.7       png_0.1-7            tibble_1.4.2        
[100] stringi_1.2.2        lattice_0.20-34      trimcluster_0.1-2   
[103] IRdisplay_0.4.4      psych_1.7.8          diffusionMap_1.1-0  
[106] pillar_1.2.1         lmtest_0.9-35        data.table_1.11.2   
[109] bitops_1.0-6         irlba_2.3.2          R6_2.3.0            
[112] latticeExtra_0.6-28  KernSmooth_2.23-15   gridExtra_2.3       
[115] codetools_0.2-15     MASS_7.3-50          gtools_3.5.0        
[118] assertthat_0.2.0     CVST_0.2-1           withr_2.1.2         
[121] mnormt_1.5-5         diptest_0.75-7       parallel_3.4.1      
[124] doSNOW_1.0.16        grid_3.4.1           rpart_4.1-13        
[127] timeDate_3012.100    IRkernel_0.8.11      tidyr_0.8.1         
[130] class_7.3-14         segmented_0.5-3.0    Cairo_1.5-9         
[133] Rtsne_0.13           pbdZMQ_0.3-2         numDeriv_2016.8-1   
[136] scatterplot3d_0.3-41 lubridate_1.7.4      base64enc_0.1-3     






