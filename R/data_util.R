

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
#' @export
readMEX <- function(mm.file, barcodes.file, genes.file) {
  count.mat = Matrix::readMM(mm.file)

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
#' @export
#'
polya_seurat_from_gene_object <- function(peak.data, genes.seurat, annot.info, project.name = "PolyA",
                                             min.cells = 10, min.peaks = 200, norm.scale.factor = 10000) {

  if (packageVersion("Seurat") < '3.0.0') {
    apa.seurat = polya_seurat_v2_from_gene_object(peak.data = peak.data, genes.seurat = genes.seurat,
                                                  annot.info = annot.info, project.name = project.name,
                                                  min.cells = min.cells, min.peaks = min.peaks,
                                                  norm.scale.factor = norm.scale.factor)
    return(apa.seurat)
  }

  ## Read in annotations to add to the Seurat object
  annot.peaks = rownames(annot.info)

  ## Check if there are annotations for peaks
  peaks.use = intersect(rownames(peak.data), annot.peaks)
  peak.data = peak.data[peaks.use, ]

  # remove any cells not in the gene-level object
  cells.keep = intersect(colnames(peak.data), colnames(genes.seurat))
  length(cells.keep)

  peak.data = peak.data[, cells.keep]

  print(paste("Creating Seurat object with", nrow(peak.data), "peaks and", ncol(peak.data), "cells"))

  ## Create a Seurat object for polyA counts
  apa.seurat = CreateSeuratObject(peak.data, min.cells = min.cells, min.features = min.peaks, project = project.name)

  ## Add peak annotations to the Seurat object
  annot.info = as.data.frame(annot.info, stringsAsFactors = FALSE)
  peaks.use = intersect(annot.peaks, rownames(Seurat::GetAssayData(apa.seurat)))
  annot.info = annot.info[peaks.use, ]
  feature.names = c("UTR3", "UTR5", "intron", "exon")
  feature.mat = annot.info[peaks.use, feature.names]

  features.collapsed = apply(feature.mat, 1, function(x) {
    paste(feature.names[which(x == "YES")], collapse = ";")})

  feature.mat$FeaturesCollapsed = features.collapsed
  feature.mat$Gene_name = annot.info$gene_id
  feature.mat$start = annot.info$start
  feature.mat$end = annot.info$end
  feature.mat$chr = annot.info$seqnames
  feature.mat$strand = annot.info$strand

  ## Add additional peak IDs for input to DEXSeq
  print("Preparing feature table for DEXSeq")
  gene.set = unique(as.character(feature.mat$Gene_name))
  dexseq.feature.table = c()
  for (this.gene in gene.set) {

    ## collect peaks
    peak.subset = subset(feature.mat, Gene_name == this.gene)
    peak.subset = peak.subset[order(peak.subset$start, decreasing = FALSE), ]

    transcript.names = paste0('transcripts "', rownames(peak.subset), '"')
    gene.ids = paste0('gene_id "', peak.subset$gene_id, '"')
    exonic.part.numbers = paste0('exonic_part_number "', 1:nrow(peak.subset), '"')
    info.part = paste(transcript.names, exonic.part.numbers, gene.ids, sep = "; ")

    dexseq.feature.set = data.frame(Gene_name = peak.subset$Gene_name,
                                   Gene_part = paste0(peak.subset$Gene_name, ":", 1:nrow(peak.subset)),
                                   Peak_number = paste0("Peak", 1:nrow(peak.subset)),
                                   Peak_name = rownames(peak.subset), stringsAsFactors = FALSE)
    rownames(dexseq.feature.set) <- dexseq.feature.set$Peak_name

    dexseq.feature.table = rbind(dexseq.feature.table, dexseq.feature.set)
  }
  dexseq.feature.table = dexseq.feature.table[rownames(feature.mat), ]

  feature.mat$Gene_part = dexseq.feature.table$Gene_part
  feature.mat$Peak_number = dexseq.feature.table$Peak_number

  ## Store the data in the Seurat @tool slot
  feature.mat.input = list(feature.mat)
  names(feature.mat.input) <- "GeneSLICER"
  apa.seurat@tools <- feature.mat.input

  ## Normalise and calculate highly-variable genes
  apa.seurat <- NormalizeData(object = apa.seurat, normalization.method = "LogNormalize",
                              scale.factor = norm.scale.factor)

  ## Add cluster identities to peak Seurat object
  cells.overlap = intersect(colnames(apa.seurat), colnames(genes.seurat))
  clusters.overlap = Idents(genes.seurat)[cells.overlap]
  clusters.overlap = clusters.overlap[colnames(apa.seurat)]
  apa.seurat = AddMetaData(object = apa.seurat, metadata = clusters.overlap, col.name = "geneLvlID")
  Idents(apa.seurat) = apa.seurat@meta.data$geneLvlID

  ## Add t-SNE coordinates to peak count object
  tryCatch({
    tsne.embeddings = Embeddings(genes.seurat, reduction = 'tsne')
    tsne.embeddings = tsne.embeddings[colnames(apa.seurat), ]
    new.embedding = CreateDimReducObject(embeddings = tsne.embeddings, key = "tSNE_", assay = "RNA")
    apa.seurat@reductions$tsne = new.embedding
    print("t-SNE coordinates added")
  }, error = function(err) {
    print("No t-SNE coodinates detected")
  })

  ## Add UMAP coordinates to peak count object
  tryCatch({
    umap.embeddings = Embeddings(genes.seurat, reduction = 'umap')
    umap.embeddings = umap.embeddings[colnames(apa.seurat), ]
    new.embedding = CreateDimReducObject(embeddings = umap.embeddings, key = "UMAP_", assay = "RNA")
    apa.seurat@reductions$umap = new.embedding
    print("UMAP coordinates added")
  }, error = function(err) {
    print("No UMAP coordinates detected")
  })

  return(apa.seurat)
}



################################################
#'
#' Create a new polyA Seurat object from the peak counts
#'
#' Creates a new polyA Seurat object from the peak counts and annotation table
#'
#' @param peak.data matrix of peak counts
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
#' @export
#'
new_polya_seurat <- function(peak.data, annot.info, project.name = "PolyA",
                                min.cells = 10, min.peaks = 200, norm.scale.factor = 10000) {

  if (packageVersion("Seurat") < '3.0.0') {
    apa.seurat = new_polya_seurat_v2(peak.data = peak.data, annot.info = annot.info,
                                                  project.name = project.name, min.cells = min.cells,
                                                  min.peaks = min.peaks, norm.scale.factor = norm.scale.factor)
    return(apa.seurat)
  }

  ## Read in annotations to add to the Seurat object
  annot.peaks = rownames(annot.info)

  ## Check if there are annotations for peaks
  peaks.use = intersect(rownames(peak.data), annot.peaks)
  peak.data = peak.data[peaks.use, ]

  print(paste("Creating Seurat object with", nrow(peak.data), "peaks and", ncol(peak.data), "cells"))

  ## Create a Seurat object for polyA counts
  apa.seurat = Seurat::CreateSeuratObject(peak.data, min.cells = min.cells, project = project.name)

  ## Add peak annotations to the Seurat object
  annot.info = as.data.frame(annot.info, stringsAsFactors = FALSE)
  #peaks.use = intersect(annot.peaks, rownames(apa.seurat@data))
  peaks.use = intersect(annot.peaks, rownames(apa.seurat))
  annot.info = annot.info[peaks.use, ]
  feature.names = c("UTR3", "UTR5", "intron", "exon")
  feature.mat = annot.info[peaks.use, feature.names]

  features.collapsed = apply(feature.mat, 1, function(x) {
    paste(feature.names[which(x == "YES")], collapse = ";")})

  feature.mat$FeaturesCollapsed = features.collapsed
  feature.mat$Gene_name = annot.info$gene_id
  feature.mat$start = annot.info$start
  feature.mat$end = annot.info$end
  feature.mat$chr = annot.info$seqnames
  feature.mat$strand = annot.info$strand

  ## Add additional peak IDs for input to DEXSeq
  print("Preparing feature table for DEXSeq")
  gene.set = unique(as.character(feature.mat$Gene_name))
  dexseq.feature.table = c()
  for (this.gene in gene.set) {

    ## collect peaks
    peak.subset = subset(feature.mat, Gene_name == this.gene)
    peak.subset = peak.subset[order(peak.subset$start, decreasing = FALSE), ]

    transcript.names = paste0('transcripts "', rownames(peak.subset), '"')
    gene.ids = paste0('gene_id "', peak.subset$gene_id, '"')
    exonic.part.numbers = paste0('exonic_part_number "', 1:nrow(peak.subset), '"')
    info.part = paste(transcript.names, exonic.part.numbers, gene.ids, sep = "; ")

    dexseq.feature.set = data.frame(Gene_name = peak.subset$Gene_name,
                                    Gene_part = paste0(peak.subset$Gene_name, ":", 1:nrow(peak.subset)),
                                    Peak_number = paste0("Peak", 1:nrow(peak.subset)),
                                    Peak_name = rownames(peak.subset), stringsAsFactors = FALSE)
    rownames(dexseq.feature.set) <- dexseq.feature.set$Peak_name

    dexseq.feature.table = rbind(dexseq.feature.table, dexseq.feature.set)
  }
  dexseq.feature.table = dexseq.feature.table[rownames(feature.mat), ]

  feature.mat$Gene_part = dexseq.feature.table$Gene_part
  feature.mat$Peak_number = dexseq.feature.table$Peak_number

  ## Store the data in the Seurat @tool slot
  feature.mat.input = list(feature.mat)
  names(feature.mat.input) <- "GeneSLICER"
  apa.seurat@tools <- feature.mat.input

  ## Normalise and calculate highly-variable genes
  apa.seurat <- Seurat::NormalizeData(object = apa.seurat, normalization.method = "LogNormalize",
                              scale.factor = norm.scale.factor)

  return(apa.seurat)
}



################################################
#'
#' Create a new peak-counts single-cell expression object from the peak counts
#'
#' Creates a new peak-counts single-cell expression object from the peak counts and annotation table
#'
#' @param peak.data matrix of peak counts
#' @param annot.info peak annotation information
#' @param cell.idents cell identities to be used for DU analysis
#' @param project.name project name passed to the Seurat object creation
#' @param min.cells minimum number of cells for retaining a peak
#' @param min.peaks minimum number of peaks for retaining a cell
#' @param norm.scale.factor scale factor for Seurat NormalizeData function
#'
#' @return a new peak-level SCE object
#'
#' @examples
#' apa.sce = new_peak_SCE(apa.data, genes.seurat, annot.info)
#'
#' @export
#'
new_peak_SCE <- function(peak.data, annot.info, cell.idents, project.name = "PolyA",
                             min.cells = 10, min.peaks = 200, norm.scale.factor = 10000) {

  ## Read in annotations to add to the SCE object
  annot.peaks = rownames(annot.info)

  names(cell.idents) <- colnames(peak.data)

  ## Check if there are annotations for peaks
  peaks.use = intersect(rownames(peak.data), annot.peaks)
  peak.data = peak.data[peaks.use, ]

  print(paste("Creating SCE object with", nrow(peak.data), "peaks and", ncol(peak.data), "cells"))

  ## filter peaks and cells
  #rows.keep <- which(rowSums(peak.data > 0) >= min.cells)
  nz.row.counts <- tabulate(peak.data@i + 1, nbins = nrow(peak.data))
  peaks.keep <- rownames(peak.data)[which(nz.row.counts >= min.cells)]

  nz.col.counts <- diff(apa.data@p)
  cells.keep <- colnames(peak.data)[which(nz.col.counts >= min.peaks)]

  ## filter the matrix and corresponding cell identities
  peak.data <- peak.data[peaks.keep, cells.keep]
  cell.idents <- cell.idents[cells.keep]

  ## Create an SCE object for peak counts
  peaks.sce <- SingleCellExperiment::SingleCellExperiment(assays = list(counts = apa.data))

  ## Add peak annotations to the SCE object
  annot.info = as.data.frame(annot.info, stringsAsFactors = FALSE)

  peaks.use = intersect(annot.peaks, rownames(peaks.sce))
  annot.info = annot.info[peaks.use, ]
  feature.names = c("UTR3", "UTR5", "intron", "exon")
  feature.mat = annot.info[peaks.use, feature.names]

  features.collapsed = apply(feature.mat, 1, function(x) {
    paste(feature.names[which(x == "YES")], collapse = ";")})

  feature.mat$FeaturesCollapsed = features.collapsed
  feature.mat$Gene_name = annot.info$gene_id
  feature.mat$start = annot.info$start
  feature.mat$end = annot.info$end
  feature.mat$chr = annot.info$seqnames
  feature.mat$strand = annot.info$strand

  ## Add additional peak IDs for input to DEXSeq
  print("Preparing feature table for DEXSeq")
  gene.set = unique(as.character(feature.mat$Gene_name))
  dexseq.feature.table = c()
  for (this.gene in gene.set) {

    ## collect peaks
    peak.subset = subset(feature.mat, Gene_name == this.gene)
    peak.subset = peak.subset[order(peak.subset$start, decreasing = FALSE), ]

    transcript.names = paste0('transcripts "', rownames(peak.subset), '"')
    gene.ids = paste0('gene_id "', peak.subset$gene_id, '"')
    exonic.part.numbers = paste0('exonic_part_number "', 1:nrow(peak.subset), '"')
    info.part = paste(transcript.names, exonic.part.numbers, gene.ids, sep = "; ")

    dexseq.feature.set = data.frame(Gene_name = peak.subset$Gene_name,
                                    Gene_part = paste0(peak.subset$Gene_name, ":", 1:nrow(peak.subset)),
                                    Peak_number = paste0("Peak", 1:nrow(peak.subset)),
                                    Peak_name = rownames(peak.subset), stringsAsFactors = FALSE)
    rownames(dexseq.feature.set) <- dexseq.feature.set$Peak_name

    dexseq.feature.table = rbind(dexseq.feature.table, dexseq.feature.set)
  }
  dexseq.feature.table = dexseq.feature.table[rownames(feature.mat), ]

  feature.mat$Gene_part = dexseq.feature.table$Gene_part
  feature.mat$Peak_number = dexseq.feature.table$Peak_number

  ## Store the data in the SCE @metadata slot
  peaks.sce@metadata$GeneSLICER <- feature.mat

  ## Add cell annotation information
  cell.data <- S4Vectors::DataFrame(CellIdent = cell.idents)
  SingleCellExperiment::colData(peaks.sce) <- cell.data

  return(peaks.sce)
}


################################################
#'
#' return polyAs associated with a select gene
#'
#' @param apa.seurat.object Seurat polyA object
#' @param gene Gene name
#' @param feature.type type of genomic features to use
#' @return a list of polyA IDs
#' @examples
#' polya.list = select_gene_polyas(apa.seurat, "PTPRC", feature.type = c("UTR3", "exon"))
#'
#' @export
#'
select_gene_polyas <- function(apa.seurat.object, gene, feature.type = c("UTR3", "UTR5", "exon", "intron")) {

  if (packageVersion("Seurat") < '3.0.0') {
    annot.subset = select_gene_polyas_v2(apa.seurat.object, gene, feature.type)
    return(annot.subset)
  }

  annot.subset = subset(Tool(apa.seurat.object, "GeneSLICER"), Gene_name == gene)
  peaks.to.use = apply(annot.subset, 1, function(x) {
    ifelse(sum(x[feature.type] == "YES") >= 1, TRUE, FALSE)
  })
  annot.subset = annot.subset[peaks.to.use, ]
  return(rownames(annot.subset))
}

