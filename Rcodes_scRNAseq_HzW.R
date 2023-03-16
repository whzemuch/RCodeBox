pp.seurat <- function(obj) {
  #singel_cells
  #@ preprocess a seurat object: Calculate the percentage of mito and filter cells
  #@ this will be used after perforing some basic EDA on seurat object and before integration

  #@ obj: seurat object
  #@ Return: a updated seurat object
  if (!inherits(obj, "Seurat")) {
    stop("obj must be a Seurat object!")
  }
  cat("\nStart to calculate the percenage of mitochondria........\n")
  obj[["percent.mt"]] <- PercentageFeatureSet(obj, pattern = "^MT")
  
  cat("\nStart to filter cells........\n")
  obj <-  subset(x = obj, 
        subset = nFeature_RNA > 200 & nFeature_RNA < 7500 & percent.mt < 12)
  cat("\nStart to normalize data with SCTtransform........\n")
  obj <-SCTransform(obj, method = "glmGamPoi", 
                    vars.to.regress = "percent.mt", verbose = FALSE)
  
  # cat("Start to regress out cell cycle genes........")
  # obj <- CellCycleScoring(obj,
  #                  s.features=cc.genes["s.genes"], 
  #                  g2m.features=cc.genes["g2m.genes"], 
  #                  set.ident=TRUE)
  # 
}


tl.run_seurat <- function(obj) {
  #single_cells
  #@ run a regular workflow on a seurat object

  #@ obj: a seurat object

  #@ Return a updated seuart object

    if (!inherits(obj, "Seurat")) {
      stop("obj must be a Seurat object!")
    }
  obj1 <- 
    obj %>% 
    RunPCA() %>% 
    FindNeighbors(dims=1:30) %>% 
    RunUMAP(dims=1:30) %>% 
    FindClusters()
} 

tl.run_monocle3 <- function(seurat_obj) {
  #single_cell
  #@ convert a seurat_obj to cell dataset and then run clustering and learing graph 
  #@ for monocle3

  #@ seurat_obj: a seuart object

  #@ Return: a cds object

  cds <- seurat_obj %>% 
    SeuratWrapper::as.cell_data_set() %>% 
    monocle3::cluster_cells(reduction_method = "UMAP") %>% 
    monocle3::learn_graph(use_partition = T) 
    
  return(cds)
}


tl.intergrate_obj_lst  <- function(seu_obj_lst) {
  #single_cell
  #@ integrated multiple seurat objects stored in a list

  #@ seu_obj_lst: a list with multiple seurat objects

  #@ Return: an integrated seurat object
  
  features <- SelectIntegrationFeatures(seu_obj_lst, nfeatures = 3000) 
  combined <- PrepSCTIntegration(seu_obj_lst, anchor.features = features) %>% 
    FindIntegrationAnchors( dims = 1:20, normalization.method = "SCT", 
                            anchor.features = features)  %>% 
    IntegrateData(normalization.method = "SCT", verbose = FALSE)
  
  return(combined)
}


remove_elems_from_vector <- function(vector, elem_lst) {
  # Remove elements from a vector
  # Deal with gene list for the feature plot when working with scRNAseq data
  
  #@ Parameters
  #@  vector:  the target vector(list)
  #@  elem_lst: the elements will be removed
  
  # return truncatedd vector
  
  has_all_elements_in_vector = all(elem_lst %in% vector)
  stopifnot(has_all_elements_in_vector )
  # https://stackoverflow.com/questions/8343509/better-error-message-for-stopifnot
  
  return (vector[!vector %in%  elem_lst])
}


plot_enrichment  <- function(marker_data, groups="cluster",
                             enrich_fun, OrgDb = org.Hs.eg.db, ...) {
  # plot the dot plot after run compareCluster function from clusterProfiler
  # https://stackoverflow.com/questions/4951442/formula-with-dynamic-number-of-variables
  # as.formula(paste("y~", paste(factors, collapse="+")))
  
  
  #@ Parameters:
  #@ OrgDb is optional, just for enrichGO function
  #@ marker_data: data frame, output of the FindAllMarker function of Seurat or 
  #  COSOG package, have gene, ENTREZID(added by full join with the output of bitr function)
  #@ enrich_fun: One of "groupGO", "enrichGO", "enrichKEGG", "enrichDO" or 
  # "enrichPathway" . Users can also supply their own function.
  #@ groups: the default is cluster, or other string or character vector
  
  #@  a compareClusterObject
  
  # separate them if there are multiple factors for grouping
  if (length(groups) > 1) {
        groups = (paste(groups, collapse = " + "))
    }
 
  geneCluster = as.formula(paste("ENTREZID ~ ", groups))
  
  geneCluster
  
  x = clusterProfiler::compareCluster(data = marker_data,
                                      OrgDb = OrgDb,
                     geneClusters = geneCluster,
                     readable=TRUE,
                     fun = enrich_fun)
 
  
  return (x)
}