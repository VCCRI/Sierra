

#######################################################################
#'
#' Apply DEXSeq to detect differential peak usage
#'
#' Apply DEXSeq to detect differential peak usage been select populations. Works by building
#' a 'pseudo-bulk' profile of cell populations by aggregating counts from individual cells
#' into a smaller number of profiles, defined by num.splits.
#'
#' @param peaks.object Either a Seurat or SCE object of peaks
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
#' DUTest(apa.seurat.object, population.1 = "1", population.2 = "2")
#'
#' @export
#'
DUTest <- function(peaks.object, population.1, population.2 = NULL, exp.thresh = 0.1,
                                     fc.thresh=0.25, adj.pval.thresh = 0.05, num.splits = 6, seed.use = 1,
                                     feature.type = c("UTR3", "exon"), verbose = TRUE,
                                     do.MAPlot = FALSE, return.dexseq.res = FALSE, ncores = 1) {

  if (class(peaks.object) == "Seurat") {
    res.table <- apply_DEXSeq_test_seurat(apa.seurat.object = peaks.object,
                                          population.1 = population.1, population.2 = population.2,
                                          exp.thresh = exp.thresh, fc.thresh = fc.thresh,
                                          adj.pval.thresh = adj.pval.thresh, num.splits = num.splits,
                                          seed.use = seed.use, feature.type = feature.type,
                                          verbose = verbose, do.MAPlot = do.MAPlot,
                                          return.dexseq.res = return.dexseq.res, ncores = ncores)
    return(res.table)

  } else if (class(peaks.object) == "SCE") {
    print ("Feature not yet implemented")
  }


}

#######################################################################
#'
#' Find alternative transcript usage between two single-cell populations
#'
#' Apply DEXSeq to detect differential peak usage been select populations. Works by building
#' a 'pseudo-bulk' profile of cell populations by aggregating counts from individual cells
#' into a smaller number of profiles, defined by num.splits.
#'
#' @param peaks.object Either a Seurat or SCE object of peaks
#' @param gtf_gr GenomicRanges object from a GTF file
#' @param gtf_TxDb TxDb from gtf file
#' @param population.1 a target population of cells (can be an ID/cluster label or a set of cell barcode IDs)
#' @param population.2 comparison population of cells. If NULL (default), uses all non-population.1 cells
#' @param exp.thresh minimum percent expression threshold (for a population of cells) to include a peak
#' @param fc.thresh threshold for log2 fold-change difference for returned results
#' @param adj.pval.thresh threshold for adjusted P-value for returned results
#' @param num.splits the number of pseudo-bulk profiles to create per identity class (default: 6)
#' @param verbose whether to print outputs (TRUE by default)
#' @param doMAPlot make an MA plot of results (FALSE by default)
#' @return a data-frame of results.
#' @examples
#' DUTest(apa.seurat.object, population.1 = "1", population.2 = "2")
#'
#' @export
#'
DetectATU <- function(peaks.object, gtf_gr, gtf_TxDb, population.1, population.2 = NULL, exp.thresh = 0.1,
                    fc.thresh=0.25, adj.pval.thresh = 0.05, num.splits = 6, seed.use = 1,
                    verbose = TRUE, do.MAPlot = FALSE, ncores = 1) {

  res.table <- DUTest(peaks.object = peaks.object, population.1 = population.1,
                      population.2 = population.2, exp.thresh = exp.thresh, fc.thresh = fc.thresh,
                      adj.pval.thresh = adj.pval.thresh, num.splits = num.splits,
                      seed.use = seed.use, feature.type = c("UTR3"),
                      verbose = verbose, do.MAPlot = do.MAPlot, ncores = ncores)

  if (verbose) print("Filtering for alternative transcript usage")

  ## format differentially used peaks from GRanges
  all.peaks = rownames(res.table)
  strand = sub(".*:.*:.*-.*:(.*)", "\\1", all.peaks)
  strand = plyr::mapvalues(x = strand, from = c("1", "-1"), to = c("+", "-"))
  peak.remainder = sub(".*:(.*:.*-.*):.*", "\\1", all.peaks)
  peaks.use = paste0(peak.remainder, ":", strand)
  res.table$granges_peaks <- peaks.use
  du.peaks.gr <- GenomicRanges::GRanges(peaks.use)

  ## format expressed peaks using GRanges
  ## get peaks that are expressed
  peaks.expressed <- get_highly_expressed_peaks(peaks.object, cluster1 = population.1,
                                                cluster2 = population.2, threshold=exp.thresh)

  genes.expressed <- sub("(.*):.*:.*-.*:.*", "\\1", peaks.expressed)
  peaks.use.idx <- which(genes.expressed %in% as.character(res.table$gene_name))
  peaks.expressed <- peaks.expressed[peaks.use.idx]

  strand = sub(".*:.*:.*-.*:(.*)", "\\1", peaks.expressed)
  strand = plyr::mapvalues(x = strand, from = c("1", "-1"), to = c("+", "-"))
  peak.remainder = sub(".*:(.*:.*-.*):.*", "\\1", peaks.expressed)
  peaks.expressed.granges = paste0(peak.remainder, ":", strand)
  expressed.peaks.gr <- GenomicRanges::GRanges(peaks.expressed.granges)

  ## make a table mapping peaks to granges peaks for later
  granges_peaks_mapping_table <- data.frame(PeakID = peaks.expressed,
                                            row.names = peaks.expressed.granges, stringsAsFactors = FALSE)

  utr3.ref <- GenomicFeatures::threeUTRsByTranscript(gtf_TxDb)
  utr3.ref <- unlist(utr3.ref)

  all_UTR_3_hits <- findOverlaps(expressed.peaks.gr , utr3.ref, type = "any")
  utr3.mappings <- as.data.frame(all_UTR_3_hits)

  query.hit.df <- as.data.frame(expressed.peaks.gr[utr3.mappings$queryHits, ])
  subject.hit.df <- as.data.frame(utr3.ref[utr3.mappings$subjectHits, ],
                                  row.names = as.character(1:nrow(utr3.mappings)))

  query.hit.df %>% mutate(granges_peak = paste0(seqnames,":",start,"-",end,":",strand)) ->
    query.hit.df

  utr3.mappings$exon_name <- subject.hit.df$exon_name
  utr3.mappings$granges_peak <- query.hit.df$granges_peak
  peak.ids <- granges_peaks_mapping_table[as.character(query.hit.df$granges_peak), 'PeakID']
  utr3.mappings$peakID <- peak.ids
  utr3.mappings %>% mutate(Gene_name = sub("(.*):.*:.*-.*:.*", "\\1", peakID)) -> utr3.mappings
  utr3.mappings %>% mutate(Start = sub(".*:(.*)-.*:.*", "\\1", granges_peak),
                           End = sub(".*:.*-(.*):.*", "\\1", granges_peak)) -> utr3.mappings

  res.table$peak_name <- rownames(res.table)
  diff.transcript.check.values <- apply(as.matrix(res.table), 1, function(x) {
    diff.site <- x["peak_name"]

    this.gene <- x["gene_name"]
    all.sites <- select_gene_polyas(peaks.object, gene = this.gene, feature.type = "UTR3")
    sites.expressed <- intersect(all.sites, peaks.expressed)

    exons.diff <- unique(subset(utr3.mappings, peakID %in% diff.site)$exon_name)

    remaining.sites <- setdiff(all.sites, diff.site)
    remaining.exons <- unique(subset(utr3.mappings, peakID %in% remaining.sites)$exon_name)

    diff.transcript.check <- ifelse(length(setdiff(remaining.exons, exons.diff)) > 0, TRUE, FALSE)
    return(diff.transcript.check)
  })

  ## Subset results table for peaks called at differential transcript usage
  res.table$Diff_transcript <- diff.transcript.check.values
  res.table <- subset(res.table, Diff_transcript == TRUE)

  ## Finally annotate the ATU peaks according to transcript name
  peaks.to.annotate <- res.table$granges_peaks
  peaks.gr <- GenomicRanges::GRanges(peaks.to.annotate)

  ## make a table mapping peaks to granges peaks for later
  granges_peaks_mapping_table <- data.frame(PeakID = rownames(res.table),
                                            row.names = peaks.to.annotate,
                                            stringsAsFactors = FALSE)

  transcripts.ref <- GenomicFeatures::transcripts(gtf_TxDb)
  transcripts.ref <- unlist(transcripts.ref)

  all_transcript_hits <- findOverlaps(peaks.gr , transcripts.ref, type = "any")
  transcript.mappings <- as.data.frame(all_transcript_hits)

  query.hit.df <- as.data.frame(peaks.gr[transcript.mappings$queryHits, ])
  subject.hit.df <- as.data.frame(transcripts.ref[transcript.mappings$subjectHits, ],
                                  row.names = as.character(1:nrow(transcript.mappings)))

  query.hit.df %>% mutate(granges_peak = paste0(seqnames,":",start,"-",end,":",strand)) ->
    query.hit.df

  transcript.mappings$transcript_name <- subject.hit.df$tx_name
  transcript.mappings$granges_peak <- query.hit.df$granges_peak
  peak.ids <- granges_peaks_mapping_table[as.character(query.hit.df$granges_peak), 'PeakID']
  transcript.mappings$peakID <- peak.ids
  transcript.mappings %>% mutate(Gene_name = sub("(.*):.*:.*-.*:.*", "\\1", peakID)) -> transcript.mappings
  transcript.mappings %>% mutate(Start = sub(".*:(.*)-.*:.*", "\\1", granges_peak),
                                 End = sub(".*:.*-(.*):.*", "\\1", granges_peak)) -> transcript.mappings

  ## Collapse the transcript names according to peak ID
  transcript.mappings %>% dplyr::group_by(peakID) %>%
    dplyr::summarize(Transcript_names = paste(transcript_name, collapse = ";")) %>%
    as.data.frame() -> peak.transcript.table
  rownames(peak.transcript.table) <- as.character(peak.transcript.table$peakID)

  ## Add the transcript names to the results table
  peak.transcript.table <- peak.transcript.table[rownames(res.table), ]
  res.table$Transcript_name <- as.character(peak.transcript.table$Transcript_names)

  res.table <- res.table[, c("gene_name", "genomic_feature(s)", "population1_pct", "population2_pct",
                             "pvalue", "padj", "Log2_fold_change", "Transcript_name")]

  if (verbose) print(paste0(nrow(res.table), " peaks detected as representing alternative transcript usage"))

  return(res.table)
}


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
#'
apply_DEXSeq_test_seurat <- function(apa.seurat.object, population.1, population.2 = NULL, exp.thresh = 0.1,
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
