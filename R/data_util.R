

################################################
#'
#' Read in peak data saved in MEX format
#'
#' Read in peak data saved in MEX format
#'
#' @param mm.file count matrix in MEX format
#' @param barcodes.file file containing cell barcodes corresponding to columns in the matrix
#' @param genes.file file containing peak or gene names corresponding to rows in the matrix
#' @return a sparseMatrix
#' @examples 
#' count.mat = readMEX(mm.file, barcodes.file, genes.file)
#'
readMEX <- function(mm.file, barcodes.file, genes.file) {
  count.mat = readMM(mm.file)
  
  barcodes = readLines(barcodes.file)
  genes = readLines(genes.file)
  
  colnames(count.mat) = barcodes
  rownames(count.mat) = genes
  
  return(count.mat)
}


################################################
#'
#' Create a polyA Seurat object using a gene-level object
#'
#' Creates a new polyA Seurat object, importing information on clustering and dimensionality reduction,
#' such as t-SNE coordinates, from a Seurat object that has been processed at the gene level. 
#' 
#' @param peak.data matrix of peak counts
#' @param genes.seurat a Seurat object
#' @param annot.info peak annotation information
#' @param project.name project name passed to the Seurat object creation
#' @param min.cells minimum number of cells for retaining a peak
#' @param min.peaks minimum number of peaks for retaining a cell
#' @param norm.scale.factor scale factor for Seurat NormalizeData function
#' 
#' @return a new polyA-level Seurat object
#' 
#' @examples
#' apa.seurat = polya_seurat_from_gene_object(apa.data, genes.seurat, annot.info)
#'
polya_seurat_from_gene_object <- function(peak.data, genes.seurat, annot.info, project.name = "PolyA",
                                   min.cells = 10, min.peaks = 200, norm.scale.factor = 10000) {
  
  ## Read in annotations to add to the Seurat object
  annot.peaks = rownames(annot.info)
  
  ## Check if there are annotations for peaks
  peaks.use = intersect(rownames(peak.data), annot.peaks)
  peak.data = peak.data[peaks.use, ]
  
  # remove any cells not in the gene-level object
  cells.keep = intersect(colnames(peak.data), genes.seurat@cell.names)
  length(cells.keep)
  
  peak.data = peak.data[, cells.keep]
  
  print(paste("Creating Seurat object with", nrow(peak.data), "peaks and", ncol(peak.data), "cells"))
  
  ## Create a Seurat object for polyA counts
  apa.seurat = CreateSeuratObject(peak.data, min.cells = min.cells, min.genes = min.peaks, project = project.name)
  
  ## Add peak annotations to the Seurat object  
  annot.info = as.data.frame(annot.info, stringsAsFactors = FALSE)
  peaks.use = intersect(annot.peaks, rownames(apa.seurat@data))
  annot.info = annot.info[peaks.use, ]
  feature.names = c("UTR3", "UTR5", "intron", "exon")
  feature.mat = annot.info[peaks.use, feature.names]
  
  features.collapsed = apply(feature.mat, 1, function(x) {
    paste(feature.names[which(x == "YES")], collapse = ";")})
  
  feature.mat$FeaturesCollapsed = features.collapsed
  feature.mat$Gene_name = annot.info$gene_id
  
  ## Store the data in the Seurat @misc slot
  apa.seurat@misc = feature.mat
  
  ## Normalise and calculate highly-variable genes
  apa.seurat <- NormalizeData(object = apa.seurat, normalization.method = "LogNormalize", 
                              scale.factor = norm.scale.factor)
  
  
  ## Add cluster identities to APA object
  cells.overlap = intersect(apa.seurat@cell.names, genes.seurat@cell.names)
  clusters.overlap = genes.seurat@ident[cells.overlap]
  clusters.overlap = clusters.overlap[apa.seurat@cell.names]
  apa.seurat = AddMetaData(object = apa.seurat, metadata = clusters.overlap, col.name = "geneLvlID")
  apa.seurat = SetAllIdent(apa.seurat, id = "geneLvlID")
  
  ## Add t-SNE coordinates to APA object
  tsne.embeddings = GetDimReduction(genes.seurat, reduction.type = 'tsne', slot = "cell.embeddings")
  tsne.embeddings = tsne.embeddings[apa.seurat@cell.names, ]
  apa.seurat = SetDimReduction(apa.seurat, reduction.type = 'tsne', slot = "cell.embeddings", new.data = tsne.embeddings)
  
  return(apa.seurat)
}

################################################
#'
#' return polyAs associated with a select gene
#' 
#' @param apa.seurat.object
#' @param gene
#' @param feature.type
#' @return a list of polyA IDs
#' @examples 
#' polya.list = select_gene_polyas(apa.seurat, "PTPRC", feature.type = c("UTR3", "exon"))
#'
select_gene_polyas <- function(apa.seurat.object, gene, feature.type = c("UTR3", "UTR5", "exon", "intron")) {
  annot.subset = subset(apa.seurat.object@misc, Gene_name == gene)
  peaks.to.use = apply(annot.subset, 1, function(x) {
    ifelse(sum(x[feature.type] == "YES") >= 1, TRUE, FALSE)
  })
  annot.subset = annot.subset[peaks.to.use, ]
  return(rownames(annot.subset))
}

