
#######################################################################
#'
#' Apply DEXSeq to detect differential peak usage
#'
#' Apply DEXSeq to detect differential peak usage been select populations. Works by building
#' a 'pseudo-bulk' profile of cell populations by aggregating counts from individual cells
#' into a smaller number of profiles, defined by num.splits.
#'
#' @param apa.seurat.object Seurat object of peaks
#' @param population.1 a target population of cells (can be an ID/cluster label or a set of cell barcode IDs)
#' @param population.2 comparison population of cells. If NULL (default), uses all non-population.1 cells
#' @param exp.thresh minimum percent expression threshold (for a population of cells) to include a peak
#' @param fc.thresh threshold for log2 fold-change difference for returned results
#' @param adj.pval.thresh threshold for adjusted P-value for returned results
#' @param num.splits the number of pseudo-bulk profiles to create per identity class (default: 6)
#' @param feature.type genomic feature types to run analysis on (degault: all)
#' @param verbose whether to print outputs (TRUE by default)
#' @param doMAPlot make an MA plot of results (FALSE by default)
#' @param return.dexseq.res return the raw and unfiltered DEXSeq results object (FALSE by default)
#' @return a data-frame of results.
#' @examples
#' apply_DEXSeq_test(apa.seurat.object, population.1 = "1", population.2 = "2")
#'
#' @export
#'
apply_DEXSeq_test <- function(apa.seurat.object, population.1, population.2 = NULL, exp.thresh = 0.1,
                              fc.thresh=0.25, adj.pval.thresh = 0.05, num.splits = 6, seed.use = 1,
                              feature.type = c("UTR3", "UTR5", "exon", "intron"), verbose = TRUE,
                              do.MAPlot = FALSE, return.dexseq.res = FALSE, ncores = 1) {

  if (!'DEXSeq' %in% rownames(x = installed.packages())) {
    stop("Please install DEXSeq before using this function
         (http://bioconductor.org/packages/release/bioc/html/DEXSeq.html)")
  }

  ## reduce counts in a cluster to num.splits cells for genes with > 1 peak
  high.expressed.peaks <- get_highly_expressed_peaks(apa.seurat.object, population.1, population.2, threshold = exp.thresh)
  length(high.expressed.peaks)

  ## Filter peaks according to feature type
  annot.subset <- Tool(apa.seurat.object, "GeneSLICER")[high.expressed.peaks, ]
  peaks.to.use <- apply(annot.subset, 1, function(x) {
    ifelse(sum(x[feature.type] == "YES") >= 1, TRUE, FALSE)
  })
  peaks.to.use <- names(peaks.to.use[which(peaks.to.use == TRUE)])
  high.expressed.peaks <- intersect(high.expressed.peaks, peaks.to.use)
  if (verbose) print(paste(length(high.expressed.peaks), "expressed peaks in feature types", toString(feature.type)))

  gene.names <- Tool(apa.seurat.object, "GeneSLICER")[high.expressed.peaks, "Gene_name"]

  ## Identifiy genes with more than one transcript detected as expressed
  gene.table <- table(gene.names)
  multi.genes <- gene.table[gene.table > 1]
  if (verbose) print(paste(length(multi.genes), "genes detected with multiple peak sites expressed"))
  multi.gene.names <- names(multi.genes)

  peaks.use <- high.expressed.peaks[which(gene.names %in% multi.gene.names)]
  if (verbose) print(paste(length(peaks.use), "Individual peak sites to test"))

  ## make pseudo-bulk profiles out of cells
  ## set a seed to allow replication of results
  set.seed(seed.use)
  if (length(population.1) == 1) {
    cells.1 <- names(Idents(apa.seurat.object))[which(Idents(apa.seurat.object) == population.1)]
  } else{
      cells.1 <- population.1
    }

  cells.1 = sample(cells.1)
  cell.sets1 <- split(cells.1, sort(1:length(cells.1)%%num.splits))

  ## create a profile set for first cluster
  profile.set1 = matrix(, nrow = length(peaks.use), ncol = length(cell.sets1))
  for (i in 1:length(cell.sets1)) {
    this.set <- cell.sets1[[i]]
    sub.matrix <- GetAssayData(apa.seurat.object, slot = "counts", assay = "RNA")[peaks.use, this.set]
    if (length(this.set) > 1) {
      this.profile <- as.numeric(apply(sub.matrix, 1, function(x) sum(x)))
      profile.set1[, i] <- this.profile
    } else {
      profile.set1[, i] <- sub.matrix
    }
  }
  rownames(profile.set1) <- peaks.use
  colnames(profile.set1) <- paste0("Population1_", 1:length(cell.sets1))

  ## create a profile set for second cluster
  if (is.null(population.2)) {
    cells.2 <- setdiff(colnames(apa.seurat.object), cells.1)
  } else {
    if (length(population.2) == 1) {
      cells.2 <- names(Idents(apa.seurat.object))[which(Idents(apa.seurat.object) == population.2)]
    } else {
      cells.2 <- population.2
    }
  }

  cells.2 = sample(cells.2)
  cell.sets2 <- split(cells.2, sort(1:length(cells.2)%%num.splits))

  profile.set2 = matrix(, nrow = length(peaks.use), ncol = length(cell.sets2))
  for (i in 1:length(cell.sets2)) {
    this.set <- cell.sets2[[i]]
    sub.matrix <- GetAssayData(apa.seurat.object, slot = "counts", assay = "RNA")[peaks.use, this.set]
    if (length(this.set) > 1) {
      this.profile <- as.numeric(apply(sub.matrix, 1, function(x) sum(x)))
      profile.set2[, i] <- this.profile
    } else {
      profile.set2[, i] <- sub.matrix
    }
  }
  rownames(profile.set2) <- peaks.use
  colnames(profile.set2) <- paste0("Population2_", 1:length(cell.sets2))

  ## merge the count matrices together
  peak.matrix <- cbind(profile.set1, profile.set2)

  ## Create the DEXSeq sample table
  sampleTable <- data.frame(row.names = c(colnames(profile.set1), colnames(profile.set2)),
                           condition = c(rep("target", ncol(profile.set1)),
                                         rep("comparison", ncol(profile.set2))))

  dexseq.feature.table <- Tool(apa.seurat.object, "GeneSLICER")[, c("Gene_name", "Gene_part", "Peak_number")]
  dexseq.feature.table$Peak <- rownames(dexseq.feature.table)
  dexseq.feature.table <- dexseq.feature.table[rownames(peak.matrix), ]
  rownames(dexseq.feature.table) <- paste0(dexseq.feature.table$Gene_name, ":", dexseq.feature.table$Peak_number)
  rownames(peak.matrix) <- rownames(dexseq.feature.table)

  peak_ID_set = dexseq.feature.table[rownames(peak.matrix), "Peak_number"]
  gene_names = dexseq.feature.table[rownames(peak.matrix), "Gene_name"]

  ## Build the DEXSeq object
  dxd = DEXSeq::DEXSeqDataSet(peak.matrix, sampleData=sampleTable, groupID = gene_names,
                              featureID = peak_ID_set, design= ~sample+exon+condition:exon)

  ## Run DEXSeq differential exon usage
  if (verbose) print("Running DEXSeq test...")

  ## Check for parallel processing option
  if (ncores > 1) {
    BPPARAM = BiocParallel::MulticoreParam(workers = ncores)
    dxd = DEXSeq::estimateSizeFactors(dxd, locfunc = genefilter::shorth)
    dxd = DEXSeq::estimateDispersions(dxd, BPPARAM = BPPARAM)
    dxd = DEXSeq::testForDEU(dxd, BPPARAM = BPPARAM)
    dxd = DEXSeq::estimateExonFoldChanges(dxd, BPPARAM = BPPARAM)
    dxr1 = DEXSeq::DEXSeqResults(dxd)

  } else {
    dxd = DEXSeq::estimateSizeFactors(dxd, locfunc = genefilter::shorth)
    dxd = DEXSeq::estimateDispersions(dxd)
    dxd = DEXSeq::testForDEU(dxd)
    dxd = DEXSeq::estimateExonFoldChanges(dxd)
    dxr1 = DEXSeq::DEXSeqResults(dxd)
  }


  if (do.MAPlot) DEXSeq::plotMA(dxr1, alpha = adj.pval.thresh,
                                ylim = c(min(dxr1$log2fold_target_comparison), max(dxr1$log2fold_target_comparison)))

  if (return.dexseq.res) return(dxr1)

  ## subset the results according to adjusted P-value and fold-change and pull out
  ## relevant data to return
  dxrSig <- subset(as.data.frame(dxr1), padj < adj.pval.thresh & abs(log2fold_target_comparison) > fc.thresh)
  dxrSig_subset <- dxrSig[, c("groupID", "exonBaseMean", "pvalue", "padj", "log2fold_target_comparison")]
  peaks.to.add = dexseq.feature.table[rownames(dxrSig_subset), "Peak"]
  rownames(dxrSig_subset) = peaks.to.add

  population.1.pct <- get_percent_expression(apa.seurat.object, population.1, remainder=FALSE, geneSet=rownames(dxrSig_subset))
  if (is.null(population.2)) {
    population.2.pct <- get_percent_expression(apa.seurat.object, population.1, remainder=TRUE, geneSet=rownames(dxrSig_subset))
    population.2 <- "Remainder"
  } else {
    population.2.pct <- get_percent_expression(apa.seurat.object, population.2, remainder=FALSE, geneSet=rownames(dxrSig_subset))
  }

  dxrSig_subset$population1_pct <- population.1.pct
  dxrSig_subset$population2_pct <- population.2.pct

  ## Add Genomic feature type
  feature.type <- Tool(apa.seurat.object, "GeneSLICER")[rownames(dxrSig_subset), c("FeaturesCollapsed")]
  dxrSig_subset$feature_type = feature.type

  dxrSig_subset <- dxrSig_subset[, c("groupID", "feature_type", "population1_pct", "population2_pct",
                                     "pvalue", "padj", "log2fold_target_comparison")]
  colnames(dxrSig_subset) <- c("gene_name", "genomic_feature(s)", "population1_pct",
                               "population2_pct", "pvalue", "padj",  "Log2_fold_change")

  dxrSig_subset <- dxrSig_subset[order(dxrSig_subset$padj, decreasing = FALSE), ]

  return(dxrSig_subset)
}




############################################################
#'
#' Calculate log2 fold-changes for a list of genes
#'
#' Get a list of Log2 fold-change values for a gene list between a specified cluster
#' and either all remaining cells (default) or an alternative cluster
#'
#' @param seurat.object A polyA Seurat object
#' @param geneList list of genes to calculate log2 fold-changes for
#' @param cluster1 target cluster
#' @param cluster2 background cluster. If null uses all non-cluster1 cells.
#' @return a list of log2 fold-changes
#' @examples
#' foldchange.list = get_log2FC_list(seurat.object, geneList, "1")
#'
get_log2FC_list <- function(seurat.object, geneList, cluster1, cluster2=NULL) {

  ### Need to check whether the input is a cluster label or vector of cell names
  if (length(cluster1) == 1){ # cluster identity used as input
    foreground.set = names(Idents(seurat.object)[Idents(seurat.object)==cluster1])
  } else { # cell identity used as input
    foreground.set = cluster1
  }
  if (is.null(cluster2)) {
    remainder.set = setdiff(colnames(seurat.object), foreground.set)
  } else {
    if (length(cluster2) == 1) { # cluster identity used as input
      remainder.set = names(Idents(seurat.object)[Idents(seurat.object)==cluster2])
    } else { # cell identity used as input
      remainder.set = cluster2
    }
  }
  log2.fc.values = apply(GetAssayData(seurat.object)[geneList, c(foreground.set, remainder.set)], 1, function(x)
    get_log2FC(x[foreground.set], x[remainder.set]))
  return(log2.fc.values)
}

############################################################
#'
#' Calculate log2 fold-changes for a list of genes
#'
#' Get a list of Log2 fold-change values for a gene list between a specified cluster
#' and either all remaining cells (default) or an alternative cluster
#'
#' @param seurat.object A polyA Seurat object
#' @param geneList list of genes to calculate log2 fold-changes for
#' @param cluster1 target cluster
#' @param cluster2 background cluster. If null uses all non-cluster1 cells.
#' @return a list of log2 fold-changes
#' @examples
#' foldchange.list = get_log2FC_list(seurat.object, geneList, "1")
#'
get_log2FC_list_v2 <- function(seurat.object, geneList, cluster1, cluster2=NULL) {

  ### Need to check whether the input is a cluster label or vector of cell names
  if (length(cluster1) == 1){ # cluster identity used as input
    foreground.set = names(seurat.object@ident[seurat.object@ident==cluster1])
  } else { # cell identity used as input
    foreground.set = cluster1
  }
  if (is.null(cluster2)) {
    remainder.set = setdiff(seurat.object@cell.names, foreground.set)
  } else {
    if (length(cluster2) == 1) { # cluster identity used as input
      remainder.set = names(seurat.object@ident[seurat.object@ident==cluster2])
    } else { # cell identity used as input
      remainder.set = cluster2
    }
  }
  log2.fc.values = apply(seurat.object@data[geneList, c(foreground.set, remainder.set)], 1, function(x)
    get_log2FC(x[foreground.set], x[remainder.set]))
  return(log2.fc.values)
}


############################################################
#'
#' Calculate Log2 fold-change
#'
#' Calculate Log2 fold-change between two vectors of log-normalised data (default Seurat normalisation)
#'
#' @param x a vector
#' @param y a vector
#' @return log2 foldchange between x and y
get_log2FC = function(x, y) {
  return(log2(mean(exp(x))/mean(exp(y))))
}

############################################################
#'
#' Calculate average log2 expression
#'
#' Calculate average expression of a peak/gene in log2 space
#'
#' @param x a vector of log-normalised data
#' @return average value of the vector in log2-space
#'
Log2ExpMean <- function (x) {
  return(log2(x = mean(x = exp(x = x) - 1) + 1))
}

############################################################
#'
#' Identify highly expressed peaks
#'
#' Selects peaks that are considered expressed above some provided criteria within a target or
#' background cluster. Considers peaks expressed in some x\% of cells to be highly expressed. Returns the
#' union of peaks identified from the target and background cluster
#'
#' @param seurat.object the
#' @param cluster1 target cluster
#' @param cluster2 background cluster. If NULL (deafult) all non-target cells
#' @param threshold percentage threshold of detected (non-zero) expression for including a peak
#' @return an array of peak (or gene) names
#' @examples
#' get_highly_expressed_peaks(seurat.object, "1")
#' get_highly_expressed_peaks(seurat.object, cluster1 = "1", cluster2 = "2")
#'
get_highly_expressed_peaks <- function(seurat.object, cluster1, cluster2=NULL, threshold=0.05) {

  if (length(cluster1) == 1){ # cluster identity used as input
    foreground.set = names(Idents(seurat.object)[Idents(seurat.object)==cluster1])
  } else { # cell identity used as input
    foreground.set = cluster1
  }
  if (is.null(cluster2)) {
    remainder.set = names(Idents(seurat.object)[Idents(seurat.object)!=cluster1])
  } else {
    if (length(cluster2) == 1) { # cluster identity used as input
      remainder.set = names(Idents(seurat.object)[Idents(seurat.object)==cluster2])
    } else { # cell identity used as input
      remainder.set = cluster2
    }
  }

  peak.names = rownames(seurat.object)

  # Get the peaks/APA sites expressed in the foreground set based on proportion of non-zeros
  this.data <- GetAssayData(seurat.object, slot = "data", assay="RNA")
  nz.row.foreground = tabulate(this.data[, foreground.set]@i + 1, nbins = nrow(seurat.object))
  nz.prop.foreground = nz.row.foreground/length(foreground.set)
  apa.foreground = peak.names[which(nz.prop.foreground > threshold)]

  # Now identify the peaks/APA sites expressed in the background set
  nz.row.background = tabulate(this.data[, remainder.set]@i + 1, nbins = nrow(seurat.object))
  nz.prop.background = nz.row.background/length(remainder.set)
  apa.background = peak.names[which(nz.prop.background > threshold)]

  return(union(apa.foreground, apa.background))
}

############################################################
#'
#' Identify highly expressed peaks
#'
#' Selects peaks that are considered expressed above some provided criteria within a target or
#' background cluster. Considers peaks expressed in some x\% of cells to be highly expressed. Returns the
#' union of peaks identified from the target and background cluster
#'
#' @param seurat.object the
#' @param cluster1 target cluster
#' @param cluster2 background cluster. If NULL (deafult) all non-target cells
#' @param threshold percentage threshold of detected (non-zero) expression for including a peak
#' @return an array of peak (or gene) names
#' @examples
#' get_highly_expressed_peaks(seurat.object, "1")
#' get_highly_expressed_peaks(seurat.object, cluster1 = "1", cluster2 = "2")
#'
get_highly_expressed_peaks_v2 <- function(seurat.object, cluster1, cluster2=NULL, threshold=0.05) {

  if (length(cluster1) == 1){ # cluster identity used as input
    foreground.set = names(seurat.object@ident[seurat.object@ident==cluster1])
  } else { # cell identity used as input
    foreground.set = cluster1
  }
  if (is.null(cluster2)) {
    remainder.set = names(seurat.object@ident[seurat.object@ident!=cluster1])
  } else {
    if (length(cluster2) == 1) { # cluster identity used as input
      remainder.set = names(seurat.object@ident[seurat.object@ident==cluster2])
    } else { # cell identity used as input
      remainder.set = cluster2
    }
  }

  peak.names = rownames(seurat.object@data)

  # Get the peaks/APA sites expressed in the foreground set based on proportion of non-zeros
  nz.row.foreground = tabulate(seurat.object@data[, foreground.set]@i + 1, nbins = length(peak.names))
  nz.prop.foreground = nz.row.foreground/length(foreground.set)
  apa.foreground = peak.names[which(nz.prop.foreground > threshold)]

  # Now identify the peaks/APA sites expressed in the background set
  nz.row.background = tabulate(seurat.object@data[, remainder.set]@i + 1, nbins = length(peak.names))
  nz.prop.background = nz.row.background/length(remainder.set)
  apa.background = peak.names[which(nz.prop.background > threshold)]

  return(union(apa.foreground, apa.background))
}

############################################################

get_percent_expression <- function(seurat.object, this.cluster, remainder=FALSE, geneSet = rownames(seurat.object)) {

  if (length(this.cluster) == 1){ # cluster identity used as input
    foreground.set = names(Idents(seurat.object)[Idents(seurat.object)==this.cluster])
  } else { # cell identity used as input
    foreground.set = this.cluster
  }

  if (remainder) {
    cell.set <- setdiff(colnames(seurat.object), foreground.set)
  } else{
    cell.set <- foreground.set
  }

  peak.names = rownames(seurat.object)

  # Get the peaks/APA sites expressed in the foreground set based on proportion of non-zeros
  this.data <- GetAssayData(seurat.object, slot = "data", assay="RNA")
  nz.row.cells = tabulate(this.data[, cell.set]@i + 1, nbins = length(peak.names))
  nz.prop.cells = nz.row.cells/length(cell.set)
  names(nz.prop.cells) = peak.names
  nz.prop.cells = nz.prop.cells[geneSet]

  return(nz.prop.cells)
}
