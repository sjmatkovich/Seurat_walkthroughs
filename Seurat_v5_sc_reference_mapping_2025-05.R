library(Seurat)
library(SeuratData)
library(ggplot2)

InstallData("panc8")

panc8 <- LoadData("panc8")
table(panc8$tech)

# we will use data from 2 technologies for the reference
pancreas.ref <- subset(panc8, tech %in% c("celseq2", "smartseq2"))
pancreas.ref[["RNA"]] <- split(pancreas.ref[["RNA"]], f = pancreas.ref$tech)

# pre-process dataset (without integration)
pancreas.ref <- NormalizeData(pancreas.ref)
pancreas.ref <- FindVariableFeatures(pancreas.ref)
pancreas.ref <- ScaleData(pancreas.ref)
pancreas.ref <- RunPCA(pancreas.ref)
pancreas.ref <- FindNeighbors(pancreas.ref, dims = 1:30)
pancreas.ref <- FindClusters(pancreas.ref)

pancreas.ref <- RunUMAP(pancreas.ref, dims = 1:30)
DimPlot(pancreas.ref, group.by = c("celltype", "tech"))

# integrate data
pancreas.ref <- IntegrateLayers(object = pancreas.ref, method = CCAIntegration, orig.reduction = "pca", new.reduction = "integrated.cca", verbose = FALSE)
pancreas.ref <- FindNeighbors(pancreas.ref, reduction = "integrated.cca", dims = 1:30)
pancreas.ref <- FindClusters(pancreas.ref)
pancreas.ref <- RunUMAP(pancreas.ref, reduction = "integrated.cca", dims = 1:30)
DimPlot(pancreas.ref, group.by = c("tech", "celltype"))

## Cell type classification using an integrated reference

# select two technologies for the query datasets
pancreas.query <- subset(panc8, tech %in% c("fluidigmc1", "celseq"))
pancreas.query <- NormalizeData(pancreas.query)
pancreas.anchors <- FindTransferAnchors(reference = pancreas.ref, query = pancreas.query, dims = 1:30,
                                        reference.reduction = "pca")
predictions <- TransferData(anchorset = pancreas.anchors, refdata = pancreas.ref$celltype, dims = 1:30)
pancreas.query <- AddMetaData(pancreas.query, metadata = predictions)

# how close are the predictions from label transfer from different techs to the known types from other techs?

pancreas.query$prediction.match <- pancreas.query$predicted.id == pancreas.query$celltype
table(pancreas.query$prediction.match)
VlnPlot(pancreas.query, c("REG1A", "PPY", "SST", "GHRL", "VWF", "SOX10"), group.by = "predicted.id", stack=T)

## Unimodal UMAP Projection

# We also enable projection of a query onto the reference UMAP structure. This can be achieved by computing the reference UMAP model and then calling MapQuery() instead of TransferData()

pancreas.ref <- RunUMAP(pancreas.ref, dims = 1:30, reduction = "integrated.cca", return.model = TRUE)
pancreas.query <- MapQuery(anchorset = pancreas.anchors, reference = pancreas.ref, query = pancreas.query, refdata = list(celltype = "celltype"), reference.reduction = "pca", reduction.model = "umap")

## What is MapQuery doing?
#   MapQuery() is a wrapper around three functions: TransferData(), IntegrateEmbeddings(), and ProjectUMAP(). TransferData() is used to transfer cell type labels and impute the ADT values; IntegrateEmbeddings() is used to integrate reference with query by correcting the query’s projected low-dimensional embeddings; and finally ProjectUMAP() is used to project the query data onto the UMAP structure of the reference. The equivalent code for doing this with the intermediate functions is below:
#   
#   pancreas.query <- TransferData(anchorset = pancreas.anchors, reference = pancreas.ref, query = pancreas.query, refdata = list(celltype = "celltype"))
# pancreas.query <- IntegrateEmbeddings(anchorset = pancreas.anchors, reference = pancreas.ref, query = pancreas.query, new.reduction.name = "ref.pca")
# pancreas.query <- ProjectUMAP(query = pancreas.query, query.reduction = "ref.pca", reference = pancreas.ref, reference.reduction = "pca", reduction.model = "umap")

p1 <- DimPlot(pancreas.ref, reduction = "umap", group.by = "celltype", label = TRUE, label.size = 3, repel = TRUE) + NoLegend() + ggtitle("Reference annotations")
p2 <- DimPlot(pancreas.query, reduction = "ref.umap", group.by = "predicted.celltype", label = TRUE, label.size = 3, repel = TRUE) + NoLegend() + ggtitle("Query transferred labels")
p1 + p2
