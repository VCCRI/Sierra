

################################################
#'
#' Read in peak data saved in MEX format
#'
#' Read in peak data saved in MEX format
#'
#' @param data.dir directory where output from CountPeaks is stored
#' @param mm.file count matrix in MEX format
#' @param barcodes.file file containing cell barcodes corresponding to columns in the matrix
#' @param sites.file file containing peak coordinate names corresponding to rows in the matrix
#' @return a sparseMatrix
#' @examples
#' 
#' 
#' library(Sierra)
#' extdata_path <- system.file("extdata",package = "Sierra")
#' reference.file <- paste0(extdata_path,"/Vignette_cellranger_genes_subset.gtf")
#' junctions.file <- paste0(extdata_path,"/Vignette_example_TIP_sham_junctions.bed")
#' bamfile <- paste0(extdata_path,"/Vignette_example_TIP_sham.bam")
#' 
#' whitelist.bc.file <- paste0(extdata_path,"/example_TIP_sham_whitelist_barcodes.tsv")
#'
#' ### Peak calling
#' peak.output.file <- c("Vignette_example_TIP_sham_peaks.txt")
#' 
#' FindPeaks(output.file = peak.output.file[1],   # output filename
#'          gtf.file = reference.file,           # gene model as a GTF file
#' bamfile = bamfile[1],                # BAM alignment filename.
#'          junctions.file = junctions.file,     # BED filename of splice junctions exising in BAM file.
#'          ncores = 1)                          # number of cores to use
#'
#'
#' #### Peak merging
#' peak.dataset.table = data.frame(Peak_file = peak.output.file,
#'                                Identifier = c("TIP-example-Sham"),
#'                                stringsAsFactors = FALSE)
#'
#' peak.merge.output.file = "TIP_merged_peaks.txt"
#' \dontrun{
#' MergePeakCoordinates(peak.dataset.table, output.file = peak.merge.output.file, ncores = 1)
#' 
#' count.dir <- c("example_TIP_sham_counts")
#' CountPeaks(peak.sites.file = peak.merge.output.file,  gtf.file = reference.file,
#'           bamfile = bamfile[1], whitelist.file = whitelist.bc.file[1],
#'           output.dir = count.dir,  countUMI = TRUE, ncores = 1)
#'             
#' 
#'  count.mat = ReadPeakCounts(count.dir)    # This will search directory and open required files
#'  
#'  
#'  }
#'  
#'  
#' @export
#'
ReadPeakCounts <- function(data.dir = NULL, mm.file = NULL, barcodes.file = NULL, sites.file = NULL) {

  if (!is.null(data.dir)) {
    mm.file <- paste0(data.dir, "/matrix.mtx")
    barcodes.file <- paste0(data.dir, "/barcodes.tsv")
    sites.file <- paste0(data.dir, "/sitenames.tsv")
  }

  count.mat <- Matrix::readMM(mm.file)

  barcodes <- readLines(barcodes.file)
  peaks <- readLines(sites.file)

  colnames(count.mat) <- barcodes
  rownames(count.mat) <- peaks

  return(count.mat)
}


################################################
#'
#' Create a peak count Seurat object using a gene-level object
#'
#' Creates a new peak Seurat object, importing information on clustering and dimensionality reduction,
#' such as t-SNE and UMAP coordinates, from a Seurat object that has been processed at the gene level.
#'
#' @param peak.data matrix of peak counts
#' @param genes.seurat a Seurat object
#' @param annot.info peak annotation information
#' @param project.name project name passed to the Seurat object creation
#' @param min.cells minimum number of cells for retaining a peak
#' @param min.peaks minimum number of peaks for retaining a cell
#' @param norm.scale.factor scale factor for Seurat NormalizeData function
#' @param filter.gene.mismatch whether to filter out peaks with ambiguous gene mappings
#'
#' @return a new peak-level Seurat object
#'
#' @examples
#'      
#' ## Load example data for two peaks from the Cxcl12 gene
#' extdata_path <- system.file("extdata",package = "Sierra")
#' load(paste0(extdata_path, "/Cxcl12_example.RData"))
#' load(paste0(extdata_path, "/TIP_cell_info.RData"))
#' 
#' ## Create an seurat object holding the peak data
#' peaks.seurat <- NewPeakSeurat(peak.data = peak.counts, 
#'                         annot.info = peak.annotations, 
#'                         cell.idents = tip.populations, 
#'                         tsne.coords = tip.tsne.coordinates,
#'                         min.cells = 0, min.peaks = 0)
#'                         
#' ##                         
#' peaks.seurat.transfer <- PeakSeuratFromTransfer(peak.data = peak.counts, 
#'                         genes.seurat = peaks.seurat, 
#'                         annot.info = peak.annotations)
#'      
#'
#' @export
#'
PeakSeuratFromTransfer <- function(peak.data, 
                                   genes.seurat, 
                                   annot.info, 
                                   project.name = "PolyA",
                                   min.cells = 10, 
                                   min.peaks = 200,
                                   norm.scale.factor = 10000,
                                   filter.gene.mismatch = TRUE) {

  if (packageVersion("Seurat") < '3.0.0') {
    stop("Seurat 3.0.0 or above is required for this function. Either upgrage or see ?NewPeakSCE")
  }

  # remove any cells not in the gene-level object
  cells.keep <- intersect(colnames(peak.data), colnames(genes.seurat))
  length(cells.keep)
  if (length(cells.keep) == 0) {
    stop("No cells overlapping with the Seurat object - please check that cell barcodes are matching.")
  }

  peak.data <- peak.data[, cells.keep]

  peaks.seurat <- NewPeakSeurat(peak.data = peak.data, 
                                annot.info = annot.info,
                                project.name = project.name, 
                                min.cells = min.cells,
                                min.peaks = min.peaks, 
                                norm.scale.factor = norm.scale.factor,
                                filter.gene.mismatch = filter.gene.mismatch,
                                verbose = FALSE)

  ## Add cluster identities to peak Seurat object
  cells.overlap <- intersect(colnames(peaks.seurat), colnames(genes.seurat))
  clusters.overlap <- Seurat::Idents(genes.seurat)[cells.overlap]
  clusters.overlap <- clusters.overlap[colnames(peaks.seurat)]
  peaks.seurat <- AddMetaData(object = peaks.seurat, metadata = clusters.overlap, col.name = "geneLvlID")
  Seurat::Idents(peaks.seurat) <- peaks.seurat@meta.data$geneLvlID

  ## Add t-SNE coordinates to peak count object
  tryCatch({
    tsne.embeddings <- Seurat::Embeddings(genes.seurat, reduction = 'tsne')
    tsne.embeddings <- tsne.embeddings[colnames(peaks.seurat), ]
    new.embedding <- Seurat::CreateDimReducObject(embeddings = tsne.embeddings, key = "tSNE_", assay = "RNA")
    peaks.seurat@reductions$tsne <- new.embedding
    print("t-SNE coordinates added")
  }, error = function(err) {
    print("No t-SNE coodinates detected")
  })

  ## Add UMAP coordinates to peak count object
  tryCatch({
    umap.embeddings <- Seurat::Embeddings(genes.seurat, reduction = 'umap')
    umap.embeddings <- umap.embeddings[colnames(peaks.seurat), ]
    new.embedding <- Seurat::CreateDimReducObject(embeddings = umap.embeddings, key = "UMAP_", assay = "RNA")
    peaks.seurat@reductions$umap = new.embedding
    print("UMAP coordinates added")
  }, error = function(err) {
    print("No UMAP coordinates detected")
  })

  return(peaks.seurat)
}



################################################
#'
#' Create a new peak-level Seurat object from the peak counts
#'
#' Creates a new peak-level Seurat object from the peak counts and annotation table
#'
#' @param peak.data matrix of peak counts
#' @param annot.info peak annotation information
#' @param project.name project name passed to the Seurat object creation
#' @param cell.idents a list of cell identities (optional)
#' @param tsne.coords a data-frame of t-SNE coordinates (optional)
#' @param umap.coords a data-frame of UMAP coordinates (optional)
#' @param min.cells minimum number of cells for retaining a peak
#' @param min.peaks minimum number of peaks for retaining a cell
#' @param norm.scale.factor scale factor for Seurat NormalizeData function
#' @param filter.gene.mismatch whether to filter out peaks with ambiguous gene mappings
#' @param verbose whether to print output 
#'
#' @return a new peak-level Seurat object
#'
#' @examples
#'                              
#' ## Load example data for two peaks from the Cxcl12 gene
#' extdata_path <- system.file("extdata",package = "Sierra")
#' load(paste0(extdata_path, "/Cxcl12_example.RData"))
#' load(paste0(extdata_path, "/TIP_cell_info.RData"))
#' 
#' ## Create an Seurat object holding the peak data
#' peaks.seurat <- NewPeakSeurat(peak.data = peak.counts, 
#'                         annot.info = peak.annotations, 
#'                         cell.idents = tip.populations, 
#'                         tsne.coords = tip.tsne.coordinates,
#'                         min.cells = 0, min.peaks = 0)
#'                              
#' 
#' @export
#'
NewPeakSeurat <- function(peak.data, annot.info, project.name = "PolyA", cell.idents = NULL,
                          tsne.coords = NULL, umap.coords = NULL, min.cells = 10,
                          min.peaks = 200, norm.scale.factor = 10000, 
                          filter.gene.mismatch = TRUE, verbose = TRUE) {

  if (packageVersion("Seurat") < '3.0.0') {
    stop("Seurat 3.0.0 or above is required for this function. Either upgrage or see ?NewPeakSCE")
  }

  if (filter.gene.mismatch) {
    ## Check that gene names according to the peak calling match gene
    ## names according to feature annotation - remove any discrepancies
    peak.gene.names = sub("(.*).*:.*:.*-.*:.*", "\\1", rownames(annot.info))
    annot.info$gene_coverage = annot.info$gene_id
    
    assigned.genes = strsplit(annot.info$gene_id[1], split = ',')[[1]]
    gene.checks = apply(annot.info, 1, function(x) {
      peak.gene = x["gene_id"]
      gene.cov = strsplit(x["gene_coverage"], split = ',')[[1]]
      sum.diff = sum(gene.cov != peak.gene)
      if (sum.diff > 0) {
        return(FALSE)
      } else {
        return(TRUE)
      }
    })
    peaks.keep = names(gene.checks[which(gene.checks == TRUE)])
    annot.info = annot.info[peaks.keep, ]
    
  }
  
  ## Peak names to add to the Seurat object
  annot.peaks <- rownames(annot.info)

  ## Check if there are annotations for peaks
  peaks.use <- intersect(rownames(peak.data), annot.peaks)
  peak.data <- peak.data[peaks.use, ]

  print(paste("Creating Seurat object with", nrow(peak.data), "peaks and", ncol(peak.data), "cells"))

  ## Create a Seurat object for polyA counts
  peaks.seurat <- Seurat::CreateSeuratObject(peak.data, min.cells = min.cells, min.features = min.peaks, project = project.name)

  ## Add cell annotation information if provided
  if (!is.null(cell.idents)) {
    cell.data <- data.frame(CellIdent = cell.idents)
    peaks.seurat <- Seurat::AddMetaData(peaks.seurat, cell.data, "CellIdent")
    Seurat::Idents(peaks.seurat) <- peaks.seurat$CellIdent
  }

  ## Add peak annotations to the Seurat object
  annot.info <- as.data.frame(annot.info, stringsAsFactors = FALSE)
  peaks.use <- intersect(annot.peaks, rownames(Seurat::GetAssayData(peaks.seurat)))
  annot.info <- annot.info[peaks.use, ]
  feature.names <- c("UTR3", "UTR5", "intron", "exon")
  feature.mat <- annot.info[peaks.use, feature.names]

  features.collapsed <- apply(feature.mat, 1, function(x) {
    paste(feature.names[which(x == "YES")], collapse = ";")})

  feature.mat$FeaturesCollapsed <- features.collapsed
  feature.mat$Gene_name <- annot.info$gene_id
  feature.mat$start <- annot.info$start
  feature.mat$end <- annot.info$end
  feature.mat$chr <- annot.info$seqnames
  feature.mat$strand <- annot.info$strand

  if (!is.null(annot.info$pA_motif)) {
    feature.mat$pA_motif <- annot.info$pA_motif
    feature.mat$pA_stretch <- annot.info$pA_stretch
  } else {
    warning("Motif information not found in annotation data - some Sierra functions will be unavailable.")
  }

  ## Add additional peak IDs for input to DEXSeq
  print("Preparing feature table for DEXSeq")
  gene.set <- unique(as.character(feature.mat$Gene_name))
  dexseq.feature.table <- c()
  for (this.gene in gene.set) {

    ## collect peaks
    peak.subset <- subset(feature.mat, Gene_name == this.gene)
    peak.subset <- peak.subset[order(peak.subset$start, decreasing = FALSE), ]

    transcript.names <- paste0('transcripts "', rownames(peak.subset), '"')
    gene.ids <- paste0('gene_id "', peak.subset$gene_id, '"')
    exonic.part.numbers <- paste0('exonic_part_number "', 1:nrow(peak.subset), '"')
    info.part <- paste(transcript.names, exonic.part.numbers, gene.ids, sep = "; ")

    dexseq.feature.set <- data.frame(Gene_name = peak.subset$Gene_name,
                                     Gene_part = paste0(peak.subset$Gene_name, ":", 1:nrow(peak.subset)),
                                     Peak_number = paste0("Peak", 1:nrow(peak.subset)),
                                     Peak_name = rownames(peak.subset), stringsAsFactors = FALSE)
    rownames(dexseq.feature.set) <- dexseq.feature.set$Peak_name

    dexseq.feature.table <- rbind(dexseq.feature.table, dexseq.feature.set)
  }
  dexseq.feature.table <- dexseq.feature.table[rownames(feature.mat), ]

  feature.mat$Gene_part <- dexseq.feature.table$Gene_part
  feature.mat$Peak_number <- dexseq.feature.table$Peak_number

  ## Store the data in the Seurat @tool slot
  feature.mat.input <- list(feature.mat)
  names(feature.mat.input) <- "Sierra"
  peaks.seurat@tools <- feature.mat.input

  ## Normalise and calculate highly-variable genes
  peaks.seurat <- Seurat::NormalizeData(object = peaks.seurat, normalization.method = "LogNormalize",
                              scale.factor = norm.scale.factor)

  ## Add t-SNE coordinates to peak count object
  if (!is.null(tsne.coords)) {
    tsne.coords <- tsne.coords[colnames(peaks.seurat), ]
    new.embedding <- Seurat::CreateDimReducObject(embeddings = tsne.coords, key = "tSNE_", assay = "RNA")
    peaks.seurat@reductions$tsne <- new.embedding
    if (verbose) print("t-SNE coordinates added")
  } else {
    if (verbose) print("No t-SNE coodinates included")
  }

  ## Add UMAP coordinates to peak count object
  if (!is.null(umap.coords)) {
    umap.coords <- umap.coords[colnames(peaks.seurat), ]
    new.embedding <- Seurat::CreateDimReducObject(embeddings = umap.coords, key = "UMAP_", assay = "RNA")
    peaks.seurat@reductions$umap = new.embedding
    if (verbose) print("UMAP coordinates added")
  } else {
    if (verbose) print("No UMAP coordinates included")
  }

  return(peaks.seurat)
}



################################################
#'
#' Create a new peak-counts single-cell experiment object from the peak counts
#'
#' Creates a new peak-counts single-cell experiment object from the peak counts and annotation table
#'
#' @param peak.data matrix of peak counts
#' @param annot.info peak annotation information
#' @param cell.idents named list of cell identities to be used for DU analysis
#' @param tsne.coords data-frame of t-SNE coordinates. Rownames should correspond to cell names.
#' @param umap.coords data-frame of UMAP coordinates. Rownames should correspond to cell names.
#' @param min.cells minimum number of cells for retaining a peak
#' @param min.peaks minimum number of peaks for retaining a cell
#' @param norm.scale.factor scale factor for log normalisation  function
#' @param filter.gene.mismatch whether to filter out peaks with ambiguous gene mappings
#' @param verbose whether to print output 
#'
#' @return a new peak-level SCE object
#'
#' @examples
#' 
#' 
#'  ## Load example data for two peaks from the Cxcl12 gene
#' extdata_path <- system.file("extdata",package = "Sierra")
#' load(paste0(extdata_path, "/Cxcl12_example.RData"))
#' load(paste0(extdata_path, "/TIP_cell_info.RData"))
#' 
#' ## Create an SCE object holding the peak data
#' peaks.sce <- NewPeakSCE(peak.data = peak.counts, 
#'                         annot.info = peak.annotations, 
#'                         cell.idents = tip.populations, 
#'                         tsne.coords = tip.tsne.coordinates,
#'                         min.cells = 0, min.peaks = 0)
#' @export
#'
#' @import SingleCellExperiment
#'
NewPeakSCE <- function(peak.data, annot.info, cell.idents = NULL, 
                       tsne.coords = NULL, umap.coords = NULL,
                       min.cells = 10, min.peaks = 200, norm.scale.factor = 10000, 
                       filter.gene.mismatch = TRUE, verbose = TRUE) {

  ## Check that peak.data is of dgCMatrix format
  if ( class(peak.data) != "dgcMatrix" )
    peak.data <- as(peak.data, "dgCMatrix")
  
  if (filter.gene.mismatch) {
    ## Check that gene names according to the peak calling match gene
    ## names according to feature annotation - remove any discrepancies
    peak.gene.names = sub("(.*).*:.*:.*-.*:.*", "\\1", rownames(annot.info))
    annot.info$gene_coverage = annot.info$gene_id
    
    assigned.genes = strsplit(annot.info$gene_id[1], split = ',')[[1]]
    gene.checks = apply(annot.info, 1, function(x) {
      peak.gene = x["gene_id"]
      gene.cov = strsplit(x["gene_coverage"], split = ',')[[1]]
      sum.diff = sum(gene.cov != peak.gene)
      if (sum.diff > 0) {
        return(FALSE)
      } else {
        return(TRUE)
      }
    })
    peaks.keep = names(gene.checks[which(gene.checks == TRUE)])
    annot.info = annot.info[peaks.keep, ]
  }

  ## Read in annotations to add to the SCE object
  annot.peaks = rownames(annot.info)

  ## Check if there are annotations for peaks
  peaks.use = intersect(rownames(peak.data), annot.peaks)
  peak.data = peak.data[peaks.use, ]

  if(verbose) print(paste("Creating SCE object with", nrow(peak.data), "peaks and", ncol(peak.data), "cells"))

  ## filter peaks and cells
  #rows.keep <- which(rowSums(peak.data > 0) >= min.cells)
  nz.row.counts <- tabulate(peak.data@i + 1, nbins = nrow(peak.data))
  peaks.keep <- rownames(peak.data)[which(nz.row.counts >= min.cells)]

  nz.col.counts <- diff(peak.data@p)
  cells.keep <- colnames(peak.data)[which(nz.col.counts >= min.peaks)]

  ## filter the matrix and corresponding cell identities
  peak.data <- peak.data[peaks.keep, cells.keep]

  if (!is.null(cell.idents)) {
    cell.idents <- cell.idents[cells.keep]
  }

  ## create a log-normalised matrix
  if(verbose) print("Log-normalising data")
  peak.data.norm <- peak.data
  peak.data.norm@x <- peak.data.norm@x / rep.int(Matrix::colSums(peak.data.norm), diff(peak.data.norm@p))
  peak.data.norm <- peak.data.norm * norm.scale.factor
  peak.data.norm@x <- log(peak.data.norm@x + 1)

  dim.reductions.list <- S4Vectors::SimpleList()

  ## check if t-SNE/UMAP coordinates have been provided
  if (!is.null(tsne.coords)) {
    tsne.coords <- tsne.coords[cells.keep, ]
    dim.reductions.list[['tsne']] <- tsne.coords
  }
  if (!is.null(umap.coords)) {
    umap.coords <- umap.coords[cells.keep, ]
    dim.reductions.list[['umap']] <- umap.coords
  }

  ## Create an SCE object for peak counts
  peaks.sce <- SingleCellExperiment::SingleCellExperiment(assays = list(counts = peak.data,
                                                                        lnorm_counts = peak.data.norm),
                                                          reducedDims = dim.reductions.list)

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

  if (!is.null(annot.info$pA_motif)) {
    feature.mat$pA_motif <- annot.info$pA_motif
    feature.mat$pA_stretch <- annot.info$pA_stretch
  } else {
    warning("Motif information not found in annotation data - some Sierra functions will be unavailable.")
  }

  ## Add additional peak IDs for input to DEXSeq
  if (verbose) print("Preparing feature table for DEXSeq")
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
  dexseq.feature.table <- dexseq.feature.table[rownames(feature.mat), ]

  feature.mat$Gene_part <- dexseq.feature.table$Gene_part
  feature.mat$Peak_number <- dexseq.feature.table$Peak_number

  ## Store the data in the SCE @metadata slot
  peaks.sce@metadata$Sierra <- feature.mat

  ## Add cell annotation information
  if (!is.null(cell.idents)) {
    cell.data <- S4Vectors::DataFrame(CellIdent = cell.idents)
    SummarizedExperiment::colData(peaks.sce) <- cell.data
  } else{
    warning("Cell identities not provided. DU testing will not be possible without these")
  }


  return(peaks.sce)
}


################################################
#'
#' Return peaks associated with a select gene.
#'
#' Returns peaks associated with a select gene.
#'
#' @param peaks.object Peaks SCE or Seurat object.
#' @param gene Gene name
#' @param feature.type type of genomic features to use
#' @return a list of peak IDs
#' @examples
#' 
#' 
#' #' extdata_path <- system.file("extdata",package = "Sierra")
#' load(paste0(extdata_path, "/Cxcl12_example.RData"))
#' load(paste0(extdata_path, "/TIP_cell_info.RData"))
#' 
#' ## Create an suerat object holding the peak data
#' peaks.seurat <- NewPeakSeurat(peak.data = peak.counts, 
#'                         annot.info = peak.annotations, 
#'                         cell.idents = tip.populations, 
#'                         tsne.coords = tip.tsne.coordinates,
#'                         min.cells = 0, min.peaks = 0)
#' 
#' peak.list <- SelectGenePeaks(peaks.object =  peaks.seurat ,gene = "Cxcl12")
#' 
#' @export
#'
SelectGenePeaks <- function(peaks.object, gene, feature.type = c("UTR3", "UTR5", "exon", "intron")) {

  if (class(peaks.object) == "Seurat") {
    annot.subset <- subset(Tool(peaks.object, "Sierra"), Gene_name == gene)
    peaks.to.use <- apply(annot.subset, 1, function(x) {
      ifelse(sum(x[feature.type] == "YES") >= 1, TRUE, FALSE)
    })
    annot.subset <- annot.subset[peaks.to.use, ]
    return(rownames(annot.subset))

  } else if (class(peaks.object) == "SingleCellExperiment") {
    annot.subset <- subset(peaks.object@metadata$Sierra, Gene_name == gene)
    peaks.to.use <- apply(annot.subset, 1, function(x) {
      ifelse(sum(x[feature.type] == "YES") >= 1, TRUE, FALSE)
    })
    annot.subset <- annot.subset[peaks.to.use, ]
    return(rownames(annot.subset))
  }
}

