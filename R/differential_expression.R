
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
#' @return a data-frame of results.
#' @examples
#' apply_DEXSeq_test(apa.seurat.object, population.1 = "1", population.2 = "2")
#'
#' @export
#'
apply_DEXSeq_test <- function(apa.seurat.object, population.1, population.2 = NULL, exp.thresh = 0.1,
                              fc.thresh=0.25, adj.pval.thresh = 0.05, num.splits = 6,
                              feature.type = c("UTR3", "UTR5", "exon", "intron"), verbose = TRUE,
                              do.MAPlot = FALSE) {

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
  set.seed(1)
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
  dxd = DEXSeq::estimateSizeFactors(dxd, locfunc = genefilter::shorth)
  dxd = DEXSeq::estimateDispersions(dxd)
  dxd = DEXSeq::testForDEU(dxd)
  dxd = DEXSeq::estimateExonFoldChanges(dxd)
  dxr1 = DEXSeq::DEXSeqResults(dxd)

  if (do.MAPlot) DEXSeq::plotMA(dxr1, alpha = adj.pval.thresh,
                                ylim = c(min(dxr1$log2fold_target_comparison), max(dxr1$log2fold_target_comparison)))

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

################################################################################################
#'
#' Apply a differential expression test to polyA sites.
#'
#' Apply a differential expression test to polyA sites, i.e. treats polyA sites as independent.
#'
#' @param apa.seurat.object a polyA Seurat object
#' @param cluster1 foreground cluster
#' @param cluster2 comparison cluster
#' @param exp.thresh threshold \% of cells expressing a peak to consider in
#' @param fc.thresh threshold log2 fold-change difference for testing a peak
#' @param adj.pval.thresh adjusted P-value threshold for retaining a peak
#' @param test.use statistical test to use. Either MAST (default) or t-test.
#' @param feature.type the feature type(s) to test for. Redundant if use.all.peaks set to TRUE.
#' @param use.all.peaks whether to use all peaks regardless of annotation. FALSE by default.
#' @return a data-frame of results.
#' @examples
#' find_de_polya(apa.seurat.object, cluster1 = "1", cluster2 = "2")
#'
find_de_polya <- function(apa.seurat.object, cluster1, cluster2, exp.thresh = 0.05,
                             fc.thresh = 0.25, adj.pval.thresh = 1e-05, test.use = "MAST",
                             feature.type = c("UTR3", "UTR5", "exon", "intron"),
                             use.all.peaks = FALSE, add.annot.info = TRUE, print.output = TRUE) {

  if (!(test.use %in% c("MAST", "t-test"))) {
    stop("Invalid test option: please choose either MAST or t-test")
  }

  if (test.use == "MAST") {
    if (!'MAST' %in% rownames(x = installed.packages())) {
      stop("Please install MAST before using this function  (https://github.com/RGLab/MAST)")
    }
  }

  if (length(cluster1) == 1){ # cluster identity used as input
    cells.fg = names(Idents(apa.seurat.object)[Idents(apa.seurat.object)==cluster1])
  } else { # cell identity used as input
    cells.fg = cluster1
  }
  if (is.null(cluster2)) {
    remainder.set = names(Idents(apa.seurat.object)[Idents(apa.seurat.object)!=cluster1])
  } else {
    if (length(cluster2) == 1) { # cluster identity used as input
      cells.bg = names(Idents(apa.seurat.object)[Idents(apa.seurat.object)==cluster2])
    } else { # cell identity used as input
      cells.bg = cluster2
    }
  }

  ## Determine the APA expressed above provided threshold proportion of cells
  if (print.output) print("Detecting expressed APA in each cluster...")
  high.expressed.peaks = get_highly_expressed_peaks(apa.seurat.object, cluster1, cluster2, threshold = exp.thresh)
  if (print.output) print(paste(length(high.expressed.peaks), "expressed peaks"))

  ## Filter peaks according to feature type
  if (use.all.peaks == FALSE) {
    annot.subset = Tool(apa.seurat.object, "GeneSLICER")[high.expressed.peaks, ]
    peaks.to.use = apply(annot.subset, 1, function(x) {
      ifelse(sum(x[feature.type] == "YES") >= 1, TRUE, FALSE)
    })
    peaks.to.use = names(peaks.to.use[which(peaks.to.use == TRUE)])
    high.expressed.peaks = intersect(high.expressed.peaks, peaks.to.use)
    if (print.output) print(paste(length(high.expressed.peaks), "expressed peaks in feature types", toString(feature.type)))
  } else {
    if (print.output) print(paste(length(high.expressed.peaks), "expressed peaks across all feature types"))
  }


  ## Calculate fold-changes for transcripts between the two populations
  all.fcs = get_log2FC_list(apa.seurat.object, high.expressed.peaks, cluster1 = cluster1, cluster2 = cluster2)

  ## Get the APAs and corresponding genes that pass a fold-change threshold
  fcs.pass = all.fcs[which(abs(all.fcs) >= fc.thresh)]
  apas.pass = names(fcs.pass)
  if (print.output) print(paste(length(apas.pass), " peaks pass LFC threshold"))

  clusters = as.character(Idents(apa.seurat.object))
  p.values = get_pvalue_list(GetAssayData(apa.seurat.object)[apas.pass, ], identities = clusters,
                             cluster1 = cluster1, cluster2 = cluster2, test.use = test.use)
  p.adj = p.adjust(p.values, method = "bonferroni", n = nrow(apa.seurat.object))


  ## Calculate the percentage of cells expressing the gene
  nz.prop.foreground <- get_percent_expression(apa.seurat.object, cells.fg)

  ## Also for the background population
  nz.prop.background <- get_percent_expression(apa.seurat.object, cells.bg)


  if (add.annot.info) {
    ## Pull out the gene names
    gene.name = Tool(apa.seurat.object, "GeneSLICER")[apas.pass, "Gene_name"]

    ## Build the results table to return
    diff.table = data.frame(Peak = apas.pass,
                            Gene = gene.name,
                            Target_pct = nz.prop.foreground[apas.pass],
                            Background_pct = nz.prop.background[apas.pass],
                            Log2FC = fcs.pass[apas.pass],
                            p_value = p.values[apas.pass],
                            P_value_adj = p.adj[apas.pass], stringsAsFactors = FALSE)

    features.add = Tool(apa.seurat.object, "GeneSLICER")[apas.pass, "FeaturesCollapsed"]
    diff.table$GenomicFeature = features.add
  } else {
    ## Don't incorporate any information from the annotation table
    diff.table = data.frame(Target_pct = nz.prop.foreground[apas.pass],
                            Background_pct = nz.prop.background[apas.pass],
                            Log2FC = fcs.pass[apas.pass],
                            p_value = p.values[apas.pass],
                            P_value_adj = p.adj[apas.pass], stringsAsFactors = FALSE)
  }
  rownames(diff.table) = apas.pass

  diff.table = diff.table[order(diff.table$P_value_adj), ]
  diff.table = subset(diff.table,  P_value_adj < adj.pval.thresh)
  if (print.output) print(paste0(nrow(diff.table), " differentialy expressed polyA sites identified"))

  return(diff.table)
}

################################################################################################
#'
#' Apply a differential expression test to polyA sites.
#'
#' Apply a differential expression test to polyA sites, i.e. treats polyA sites as independent.
#'
#' @param apa.seurat.object a polyA Seurat object
#' @param cluster1 foreground cluster
#' @param cluster2 comparison cluster
#' @param exp.thresh threshold \% of cells expressing a peak to consider in
#' @param fc.thresh threshold log2 fold-change difference for testing a peak
#' @param adj.pval.thresh adjusted P-value threshold for retaining a peak
#' @param test.use statistical test to use. Either MAST (default) or t-test.
#' @param feature.type the feature type(s) to test for. Redundant if use.all.peaks set to TRUE.
#' @param use.all.peaks whether to use all peaks regardless of annotation. FALSE by default.
#' @return a data-frame of results.
#' @examples
#' find_de_polya(apa.seurat.object, cluster1 = "1", cluster2 = "2")
#'
find_de_polya_v2 <- function(apa.seurat.object, cluster1, cluster2, exp.thresh = 0.05,
                          fc.thresh = 0.25, adj.pval.thresh = 1e-05, test.use = "MAST",
                          feature.type = c("UTR3", "UTR5", "exon", "intron"),
                          use.all.peaks = FALSE, add.annot.info = TRUE, print.output = TRUE) {

  if (!(test.use %in% c("MAST", "t-test"))) {
    stop("Invalid test option: please choose either MAST or t-test")
  }

  if (test.use == "MAST") {
    if (!'MAST' %in% rownames(x = installed.packages())) {
      stop("Please install MAST before using this function  (https://github.com/RGLab/MAST)")
    }
  }

  if (length(cluster1) == 1){ # cluster identity used as input
    cells.fg = names(apa.seurat.object@ident[apa.seurat.object@ident==cluster1])
  } else { # cell identity used as input
    cells.fg = cluster1
  }
  if (is.null(cluster2)) {
    remainder.set = names(apa.seurat.object@ident[apa.seurat.object@ident!=cluster1])
  } else {
    if (length(cluster2) == 1) { # cluster identity used as input
      cells.bg = names(apa.seurat.object@ident[apa.seurat.object@ident==cluster2])
    } else { # cell identity used as input
      cells.bg = cluster2
    }
  }

  ## Determine the APA expressed above provided threshold proportion of cells
  if (print.output) print("Detecting expressed APA in each cluster...")
  high.expressed.peaks = get_highly_expressed_peaks(apa.seurat.object, cluster1, cluster2, threshold = exp.thresh)
  if (print.output) print(paste(length(high.expressed.peaks), "expressed peaks"))

  ## Filter peaks according to feature type
  if (use.all.peaks == FALSE) {
    annot.subset = apa.seurat.object@misc[high.expressed.peaks, ]
    peaks.to.use = apply(annot.subset, 1, function(x) {
      ifelse(sum(x[feature.type] == "YES") >= 1, TRUE, FALSE)
    })
    peaks.to.use = names(peaks.to.use[which(peaks.to.use == TRUE)])
    high.expressed.peaks = intersect(high.expressed.peaks, peaks.to.use)
    if (print.output) print(paste(length(high.expressed.peaks), "expressed peaks in feature types", toString(feature.type)))
  } else {
    if (print.output) print(paste(length(high.expressed.peaks), "expressed peaks across all feature types"))
  }


  ## Calculate fold-changes for transcripts between the two populations
  all.fcs = get_log2FC_list(apa.seurat.object, high.expressed.peaks, cluster1 = cluster1, cluster2 = cluster2)

  ## Get the APAs and corresponding genes that pass a fold-change threshold
  fcs.pass = all.fcs[which(abs(all.fcs) >= fc.thresh)]
  apas.pass = names(fcs.pass)
  if (print.output) print(paste(length(apas.pass), " peaks pass LFC threshold"))

  clusters = as.character(apa.seurat.object@ident)
  p.values = get_pvalue_list(apa.seurat.object@data[apas.pass, ], identities = clusters,
                           cluster1 = cluster1, cluster2 = cluster2, test.use = test.use)
  p.adj = p.adjust(p.values, method = "bonferroni", n = nrow(apa.seurat.object@data))



  ## Calculate the percentage of cells expressing the gene
  nz.row.foreground = tabulate(apa.seurat.object@data[, cells.fg]@i + 1)
  nz.prop.foreground = nz.row.foreground/length(cells.fg)
  names(nz.prop.foreground) = rownames(apa.seurat.object@data)

  ## Also for the background population
  nz.row.background = tabulate(apa.seurat.object@data[, cells.bg]@i + 1)
  nz.prop.background = nz.row.background/length(cells.bg)
  names(nz.prop.background) = rownames(apa.seurat.object@data)

  if (add.annot.info) {
    ## Pull out the gene names
    gene.name = apa.seurat.object@misc[apas.pass, "Gene_name"]

    ## Build the results table to return
    diff.table = data.frame(Peak = apas.pass,
                            Gene = gene.name,
                            Target_pct = nz.prop.foreground[apas.pass],
                            Background_pct = nz.prop.background[apas.pass],
                            Log2FC = fcs.pass[apas.pass],
                            p_value = p.values[apas.pass],
                            P_value_adj = p.adj[apas.pass], stringsAsFactors = FALSE)

    features.add = apa.seurat.object@misc[apas.pass, "FeaturesCollapsed"]
    diff.table$GenomicFeature = features.add
  } else {
    ## Don't incorporate any information from the annotation table
    diff.table = data.frame(Target_pct = nz.prop.foreground[apas.pass],
                            Background_pct = nz.prop.background[apas.pass],
                            Log2FC = fcs.pass[apas.pass],
                            p_value = p.values[apas.pass],
                            P_value_adj = p.adj[apas.pass], stringsAsFactors = FALSE)
  }
  rownames(diff.table) = apas.pass

  diff.table = diff.table[order(diff.table$P_value_adj), ]
  diff.table = subset(diff.table,  P_value_adj < adj.pval.thresh)
  if (print.output) print(paste0(nrow(diff.table), " differentialy expressed polyA sites identified"))

  return(diff.table)
}

#############################################################
#'
#' Find genes showing differential usage of peak sites
#'
#' First run a differential expression test to identify peaks differentially expressed
#' between clusters. Selects genes that have at least one differential peak, but expressing
#' multiple peaks.
#'
#' @param apa.seurat.object a polyA Seurat object
#' @param cluster1 foreground cluster
#' @param cluster2 comparison cluster
#' @param exp.thresh threshold \% of cells expressing a peak to consider in
#' @param fc.thresh threshold log2 fold-change difference for testing a peak
#' @param adj.pval.thresh adjusted P-value threshold for retaining a peak
#' @param test.use statistical test to use. Either MAST (default) or t-test.
#' @param feature.type the feature type(s) to test for. Redundant if use.all.peaks set to TRUE.
#' @param use.all.peaks whether to use all peaks regardless of annotation. FALSE by default.
#' @param add.annot.info whether to add genomic annotation to output (TRUE by default)
#' @param print.output whether to print the output
#' @return a data-frame of results.
#' @examples
#' find_de_polya(apa.seurat.object, cluster1 = "1", cluster2 = "2")
#'
find_max_polyA_de <- function(apa.seurat.object, cluster1, cluster2 = NULL, exp.thresh = 0.1,
                                 fc.thresh = 0.5, peak.fc.thresh = 0.15, adj.pval.thresh = 1e-05,
                                 test.use = "t-test", feature.type = c("UTR3", "UTR5", "exon", "intron"),
                                 use.all.peaks = FALSE, add.annot.info = TRUE, print.output = TRUE){

  if (length(cluster1) == 1){ # cluster identity used as input
    cells.fg = names(Idents(apa.seurat.object)[Idents(apa.seurat.object)==cluster1])
  } else { # cell identity used as input
    cells.fg = cluster1
  }
  if (is.null(cluster2)) {
    cells.bg = names(Idents(apa.seurat.object)[Idents(apa.seurat.object)!=cluster1])
  } else {
    if (length(cluster2) == 1) { # cluster identity used as input
      cells.bg = names(Idents(apa.seurat.object)[Idents(apa.seurat.object)==cluster2])
    } else { # cell identity used as input
      cells.bg = cluster2
    }
  }

  ## Determine the APA expressed above provided threshold proportion of cells
  if (print.output) print("Detecting expressed APA in each cluster...")
  high.expressed.peaks = get_highly_expressed_peaks(apa.seurat.object, cluster1, cluster2, threshold = exp.thresh)
  if (print.output) print(paste(length(high.expressed.peaks), "expressed peaks"))

  ## Filter peaks according to feature type
  if (use.all.peaks == FALSE) {
    annot.subset = Tool(apa.seurat.object, "GeneSLICER")[high.expressed.peaks, ]
    peaks.to.use = apply(annot.subset, 1, function(x) {
      ifelse(sum(x[feature.type] == "YES") >= 1, TRUE, FALSE)
    })
    peaks.to.use = names(peaks.to.use[which(peaks.to.use == TRUE)])
    high.expressed.peaks = intersect(high.expressed.peaks, peaks.to.use)
    if (print.output) print(paste(length(high.expressed.peaks), "expressed peaks in feature types", toString(feature.type)))
  } else {
    if (print.output) print(paste(length(high.expressed.peaks), "expressed peaks across all feature types"))
  }

  ## Pull out the unique gene names - should be in a geneA.1, geneA.2 etc. format
  gene.names = Tool(apa.seurat.object, "GeneSLICER")[high.expressed.peaks, "Gene_name"]

  ## Identifiy genes with more than one transcript detected as expressed
  gene.table = table(gene.names)
  multi.genes = gene.table[gene.table > 1]
  print(paste(length(multi.genes), "genes detected with multiple polyA sites expressed"))
  multi.gene.names = names(multi.genes)

  peaks.use = high.expressed.peaks[which(gene.names %in% multi.gene.names)]
  print(paste(length(peaks.use), "Individual polyA sites to test"))

  ## Calculate fold-changes for transcripts between the two populations
  all.fcs = get_log2FC_list(apa.seurat.object, peaks.use, cluster1 = cluster1, cluster2 = cluster2)

  ## Calculate P-values (t-test) between the two populations
  print("Calculating differential polyA usage between populations...")
  clusters = as.character(Idents(apa.seurat.object))
  p.values = get_pvalue_list(GetAssayData(apa.seurat.object)[peaks.use, ], identities = clusters,
                             cluster1 = cluster1, cluster2 = cluster2, test.use = test.use)
  all.adj.pvalues = p.adjust(p.values, method = "bonferroni", n = nrow(apa.seurat.object))

  ## Now go through genes and determine examples by differential polyA usage
  sub.matrix = GetAssayData(apa.seurat.object)[peaks.use, ]
  sub.matrix = sub.matrix[, cells.fg]
  res.table = c()
  for (thisGene in multi.gene.names) {
    ## Pull out the transcripts corresponding to the gene
    these.peaks = high.expressed.peaks[which(gene.names == thisGene)]

    ## Obtain the adjusted P-values for the transcripts
    adj.pvalues = all.adj.pvalues[these.peaks]

    ## Obtain fold-changes between populations for the select transcripts
    fcs = all.fcs[these.peaks]

    ## Check if fold-change difference is greater than a chosen threshold
    ## And if there is significant differential expression between the two populations
    if ((max(fcs) - min(fcs)) > fc.thresh & sum(adj.pvalues < adj.pval.thresh) > 0) {

      ## Finally compare expression of the two transcripts within the population of interest
      peak1 = these.peaks[which(fcs == max(fcs))]
      peak2 = these.peaks[which(fcs == min(fcs))]
      this.pval = t.test(x = sub.matrix[peak1, cells.fg],
                         y = sub.matrix[peak2, cells.fg])$p.value
      if (is.nan(this.pval)) {
        this.pval = 1
      }
      this.adj.pval = this.pval * length(multi.gene.names)

      if (this.adj.pval < adj.pval.thresh) {
        this.res = data.frame(Gene = thisGene, Peaks = paste(these.peaks, collapse = "_"),
                              Max_peak = peak1, Peak1_Log2FC = max(fcs),
                              Peak1_Adj_Pval = adj.pvalues[peak1], Min_peak = peak2,
                              Peak2_Log2FC = min(fcs), Peak2_Adj_Pval = adj.pvalues[peak2],
                              Log2FC_diff = max(fcs) - min(fcs), peak_compare_adjPval = this.adj.pval)
        res.table = rbind(res.table, this.res)
      }
    }
  }

  ## Add usage classifications to results table
  usage.classification = rep("", nrow(res.table))
  usage.classification[which(res.table$Peak1_Log2FC > peak.fc.thresh &
                               res.table$Peak2_Log2FC < (-1 * peak.fc.thresh) &
                               res.table$Peak1_Adj_Pval < adj.pval.thresh &
                               res.table$Peak2_Adj_Pval < adj.pval.thresh)] <- "Switching"

  usage.classification[which(res.table$Peak1_Log2FC < (-1 * peak.fc.thresh) &
                               res.table$Peak2_Log2FC > peak.fc.thresh &
                               res.table$Peak1_Adj_Pval < adj.pval.thresh &
                               res.table$Peak2_Adj_Pval < adj.pval.thresh)] <- "Switching"

  usage.classification[which(res.table$Peak1_Log2FC > peak.fc.thresh &
                               res.table$Peak1_Adj_Pval < adj.pval.thresh &
                               (res.table$Peak2_Adj_Pval > adj.pval.thresh |
                                  abs(res.table$Peak2_Log2FC) < peak.fc.thresh))] <- "Peak1 up"

  usage.classification[which(res.table$Peak1_Log2FC < (-1 * peak.fc.thresh) &
                               res.table$Peak1_Adj_Pval < adj.pval.thresh &
                               (res.table$Peak2_Adj_Pval > adj.pval.thresh |
                                  abs(res.table$Peak2_Log2FC) < peak.fc.thresh))] <- "Peak1 down"

  usage.classification[which(res.table$Peak2_Log2FC > peak.fc.thresh &
                               res.table$Peak2_Adj_Pval < adj.pval.thresh &
                               (res.table$Peak1_Adj_Pval > adj.pval.thresh |
                                  abs(res.table$Peak1_Log2FC) < peak.fc.thresh))] <- "Peak2 up"

  usage.classification[which(res.table$Peak2_Log2FC < (-1 * peak.fc.thresh) &
                               res.table$Peak2_Adj_Pval < adj.pval.thresh &
                               (res.table$Peak1_Adj_Pval > adj.pval.thresh |
                                  abs(res.table$Peak1_Log2FC) < peak.fc.thresh))] <- "Peak2 down"

  usage.classification[which(res.table$Peak1_Log2FC < (-1 * peak.fc.thresh) &
                               res.table$Peak2_Log2FC < (-1 * peak.fc.thresh) &
                               res.table$Peak1_Adj_Pval < adj.pval.thresh &
                               res.table$Peak2_Adj_Pval < adj.pval.thresh)] <- "Both down"

  usage.classification[which(res.table$Peak1_Log2FC > peak.fc.thresh &
                               res.table$Peak2_Log2FC > peak.fc.thresh &
                               res.table$Peak1_Adj_Pval < adj.pval.thresh &
                               res.table$Peak2_Adj_Pval < adj.pval.thresh)] <- "Both up"

  res.table$Usage_class = usage.classification

  print(paste(nrow(res.table), " genes with differential peak expression identified"))
  ## Return the results table
  return(res.table)
}

#############################################################
#'
#' Find genes showing differential usage of peak sites
#'
#' First run a differential expression test to identify peaks differentially expressed
#' between clusters. Selects genes that have at least one differential peak, but expressing
#' multiple peaks.
#'
#' @param apa.seurat.object a polyA Seurat object
#' @param cluster1 foreground cluster
#' @param cluster2 comparison cluster
#' @param exp.thresh threshold \% of cells expressing a peak to consider in
#' @param fc.thresh threshold log2 fold-change difference for testing a peak
#' @param adj.pval.thresh adjusted P-value threshold for retaining a peak
#' @param test.use statistical test to use. Either MAST (default) or t-test.
#' @param feature.type the feature type(s) to test for. Redundant if use.all.peaks set to TRUE.
#' @param use.all.peaks whether to use all peaks regardless of annotation. FALSE by default.
#' @param add.annot.info whether to add genomic annotation to output (TRUE by default)
#' @param print.output whether to print the output
#' @return a data-frame of results.
#' @examples
#' find_de_polya(apa.seurat.object, cluster1 = "1", cluster2 = "2")
#'
find_max_polyA_de_v2 <- function(apa.seurat.object, cluster1, cluster2 = NULL, exp.thresh = 0.1,
                              fc.thresh = 0.5, peak.fc.thresh = 0.15, adj.pval.thresh = 1e-05,
                              test.use = "MAST", feature.type = c("UTR3", "UTR5", "exon", "intron"),
                              use.all.peaks = FALSE, add.annot.info = TRUE, print.output = TRUE){

  ## Pull out the set of cells for performing the comparison
  cells.fg = apa.seurat.object@cell.names[which(apa.seurat.object@ident == cluster1)]

  if (is.null(cluster2)){ ## if cluster 2 not provided just use all remaining cells
    cells.bg = apa.seurat.object@cell.names[which(apa.seurat.object@ident != cluster1)]
  } else {
    cells.bg = apa.seurat.object@cell.names[which(apa.seurat.object@ident == cluster2)]
  }

  ## Determine the APA expressed above provided threshold proportion of cells
  if (print.output) print("Detecting expressed APA in each cluster...")
  high.expressed.peaks = get_highly_expressed_peaks(apa.seurat.object, cluster1, cluster2, threshold = exp.thresh)
  if (print.output) print(paste(length(high.expressed.peaks), "expressed peaks"))

  ## Filter peaks according to feature type
  if (use.all.peaks == FALSE) {
    annot.subset = apa.seurat.object@misc[high.expressed.peaks, ]
    peaks.to.use = apply(annot.subset, 1, function(x) {
      ifelse(sum(x[feature.type] == "YES") >= 1, TRUE, FALSE)
    })
    peaks.to.use = names(peaks.to.use[which(peaks.to.use == TRUE)])
    high.expressed.peaks = intersect(high.expressed.peaks, peaks.to.use)
    if (print.output) print(paste(length(high.expressed.peaks), "expressed peaks in feature types", toString(feature.type)))
  } else {
    if (print.output) print(paste(length(high.expressed.peaks), "expressed peaks across all feature types"))
  }

  ## Pull out the unique gene names - should be in a geneA.1, geneA.2 etc. format
  gene.names = apa.seurat.object@misc[high.expressed.peaks, "Gene_name"]

  ## Identifiy genes with more than one transcript detected as expressed
  gene.table = table(gene.names)
  multi.genes = gene.table[gene.table > 1]
  print(paste(length(multi.genes), "genes detected with multiple polyA sites expressed"))
  multi.gene.names = names(multi.genes)

  peaks.use = high.expressed.peaks[which(gene.names %in% multi.gene.names)]
  print(paste(length(peaks.use), "Individual polyA sites to test"))

  ## Calculate fold-changes for transcripts between the two populations
  all.fcs = get_log2FC_list(apa.seurat.object, peaks.use, cluster1 = cluster1, cluster2 = cluster2)

  ## Calculate P-values (t-test) between the two populations
  print("Calculating differential polyA usage between populations...")
  clusters = as.character(apa.seurat.object@ident)
  p.values = get_pvalue_list(apa.seurat.object@data[peaks.use, ], identities = clusters,
                             cluster1 = cluster1, cluster2 = cluster2, test.use = test.use)
  all.adj.pvalues = p.adjust(p.values, method = "bonferroni", n = nrow(apa.seurat.object@data))

  ## Now go through genes and determine examples by differential polyA usage
  sub.matrix = apa.seurat.object@data[peaks.use, ]
  sub.matrix = sub.matrix[, cells.fg]
  res.table = c()
  for (thisGene in multi.gene.names) {
    ## Pull out the transcripts corresponding to the gene
    these.peaks = high.expressed.peaks[which(gene.names == thisGene)]

    ## Obtain the adjusted P-values for the transcripts
    adj.pvalues = all.adj.pvalues[these.peaks]

    ## Obtain fold-changes between populations for the select transcripts
    fcs = all.fcs[these.peaks]

    ## Check if fold-change difference is greater than a chosen threshold
    ## And if there is significant differential expression between the two populations
    if ((max(fcs) - min(fcs)) > fc.thresh & sum(adj.pvalues < adj.pval.thresh) > 0) {

      ## Finally compare expression of the two transcripts within the population of interest
      peak1 = these.peaks[which(fcs == max(fcs))]
      peak2 = these.peaks[which(fcs == min(fcs))]
      this.pval = t.test(x = sub.matrix[peak1, cells.fg],
                         y = sub.matrix[peak2, cells.fg])$p.value
      if (is.nan(this.pval)) {
        this.pval = 1
      }
      this.adj.pval = this.pval * length(multi.gene.names)

      if (this.adj.pval < adj.pval.thresh) {
        this.res = data.frame(Gene = thisGene, Peaks = paste(these.peaks, collapse = "_"),
                              Max_peak = peak1, Peak1_Log2FC = max(fcs),
                              Peak1_Adj_Pval = adj.pvalues[peak1], Min_peak = peak2,
                              Peak2_Log2FC = min(fcs), Peak2_Adj_Pval = adj.pvalues[peak2],
                              Log2FC_diff = max(fcs) - min(fcs), peak_compare_adjPval = this.adj.pval)
        res.table = rbind(res.table, this.res)
      }
    }
  }

  ## Add usage classifications to results table
  usage.classification = rep("", nrow(res.table))
  usage.classification[which(res.table$Peak1_Log2FC > peak.fc.thresh &
                             res.table$Peak2_Log2FC < (-1 * peak.fc.thresh) &
                             res.table$Peak1_Adj_Pval < adj.pval.thresh &
                             res.table$Peak2_Adj_Pval < adj.pval.thresh)] <- "Switching"

  usage.classification[which(res.table$Peak1_Log2FC < (-1 * peak.fc.thresh) &
                             res.table$Peak2_Log2FC > peak.fc.thresh &
                             res.table$Peak1_Adj_Pval < adj.pval.thresh &
                             res.table$Peak2_Adj_Pval < adj.pval.thresh)] <- "Switching"

  usage.classification[which(res.table$Peak1_Log2FC > peak.fc.thresh &
                             res.table$Peak1_Adj_Pval < adj.pval.thresh &
                             (res.table$Peak2_Adj_Pval > adj.pval.thresh |
                             abs(res.table$Peak2_Log2FC) < peak.fc.thresh))] <- "Peak1 up"

  usage.classification[which(res.table$Peak1_Log2FC < (-1 * peak.fc.thresh) &
                             res.table$Peak1_Adj_Pval < adj.pval.thresh &
                             (res.table$Peak2_Adj_Pval > adj.pval.thresh |
                             abs(res.table$Peak2_Log2FC) < peak.fc.thresh))] <- "Peak1 down"

  usage.classification[which(res.table$Peak2_Log2FC > peak.fc.thresh &
                             res.table$Peak2_Adj_Pval < adj.pval.thresh &
                             (res.table$Peak1_Adj_Pval > adj.pval.thresh |
                             abs(res.table$Peak1_Log2FC) < peak.fc.thresh))] <- "Peak2 up"

  usage.classification[which(res.table$Peak2_Log2FC < (-1 * peak.fc.thresh) &
                             res.table$Peak2_Adj_Pval < adj.pval.thresh &
                             (res.table$Peak1_Adj_Pval > adj.pval.thresh |
                             abs(res.table$Peak1_Log2FC) < peak.fc.thresh))] <- "Peak2 down"

  usage.classification[which(res.table$Peak1_Log2FC < (-1 * peak.fc.thresh) &
                             res.table$Peak2_Log2FC < (-1 * peak.fc.thresh) &
                             res.table$Peak1_Adj_Pval < adj.pval.thresh &
                             res.table$Peak2_Adj_Pval < adj.pval.thresh)] <- "Both down"

  usage.classification[which(res.table$Peak1_Log2FC > peak.fc.thresh &
                             res.table$Peak2_Log2FC > peak.fc.thresh &
                             res.table$Peak1_Adj_Pval < adj.pval.thresh &
                             res.table$Peak2_Adj_Pval < adj.pval.thresh)] <- "Both up"

  res.table$Usage_class = usage.classification

  print(paste(nrow(res.table), " genes with differential peak expression identified"))
  ## Return the results table
  return(res.table)
}

################################################################################################
#'
#' Apply differential expression testing to gene-level counts
#'
#' Essentially a wrapper function for find_de_polya. Applies the same approach for DE testing,
#' but without any reference to annotation information. Thus, this function can be applied to
#' a gene-level Seurat object to obtain differentially expressed genes for a direct comparison
#' with results from differential testing at the polyA level.
#'
#'
#' @param apa.seurat.object a gene-level Seurat object
#' @param cluster1 foreground cluster
#' @param cluster2 comparison cluster
#' @param exp.thresh threshold \% of cells expressing a gene to consider in
#' @param fc.thresh threshold log2 fold-change difference for testing a gene
#' @param adj.pval.thresh adjusted P-value threshold for retaining a gene
#' @param test.use statistical test to use. Either MAST (default) or t-test.
#' @return a data-frame of results.
#' @examples
#' find_de_genes(apa.seurat.object, cluster1 = "1", cluster2 = "2")
#'
find_de_genes <- function(genes.seurat.object, cluster1, cluster2, exp.thresh = 0.05,
                          fc.thresh = 0.25, adj.pval.thresh = 1e-05, test.use = "MAST") {

  ## Wrapper for find_de_polya
  res.table.genes = find_de_polya(genes.seurat.object, cluster1 = cluster1, cluster2 = cluster2,
                                  exp.thresh = exp.thresh, fc.thresh = fc.thresh,
                                  adj.pval.thresh = adj.pval.thresh, test.use = test.use,
                                  use.all.peaks = TRUE, add.annot.info = FALSE,
                                  print.output = FALSE)

  print(paste0(nrow(res.table.genes), " differentialy expressed genes identified"))
  return(res.table.genes)

}

################################################################################################
#'
#' Apply a differential expression test to gene-level collapsed polyA sites.
#'
#' Collapses polyA sites into a single gene-level count for each gene. Performs differential
#' gene expression analysis on the collapsed gene-level counts.
#'
#'
#' @param apa.seurat.object a polyA Seurat object
#' @param cluster1 foreground cluster
#' @param cluster2 comparison cluster
#' @param exp.thresh threshold \% of cells expressing a peak to consider in
#' @param fc.thresh threshold log2 fold-change difference for testing a peak
#' @param adj.pval.thresh adjusted P-value threshold for retaining a peak
#' @param test.use statistical test to use. Either MAST (default) or t-test.
#' @param feature.type the feature type(s) to test for. Redundant if use.all.peaks set to TRUE.
#' @param use.all.peaks whether to use all peaks regardless of annotation. FALSE by default.
#' @return a data-frame of results.
#' @examples
#' find_de_genes_collapsed(apa.seurat.object, cluster1 = "1", cluster2 = "2")
#'
#' @export
#'
find_de_genes_collapsed <- function(apa.seurat.object, cluster1, cluster2, exp.thresh = 0.05,
                                    fc.thresh = 0.25, adj.pval.thresh = 1e-05, test.use = "MAST",
                                    feature.type = c("UTR3", "UTR5", "exon", "intron"), use.all.peaks = FALSE) {

  if (packageVersion("Seurat") < '3.0.0') {
    res.table = find_de_genes_collapsed_v2(apa.seurat.oject, cluster1, cluster2, exp.thresh,
                                           fc.thresh, adj.pval.thresh, test.use, feature.type, use.all.peaks)
    return(res.table)
  }

  if (!(test.use %in% c("MAST", "t-test"))) {
    stop("Invalid test option: please choose either MAST or t-test")
  }

  if (test.use == "MAST") {
    if (!'MAST' %in% rownames(x = installed.packages())) {
      stop("Please install MAST before using this function  (https://github.com/RGLab/MAST)")
    }
  }

  ## Pull out the set of cells for performing the comparison
  cells.fg = colnames(apa.seurat.object)[which(Idents(apa.seurat.object) == cluster1)]
  cells.bg = colnames(apa.seurat.object)[which(Idents(apa.seurat.object) == cluster2)]

  ## Determine the APA expressed above provided threshold proportion of cells
  print("Detecting expressed APA in each cluster...")
  high.expressed.peaks = get_highly_expressed_peaks(apa.seurat.object, cluster1, cluster2, threshold = exp.thresh)
  print(paste(length(high.expressed.peaks), "expressed peaks"))

  ## Filter peaks according to feature type
  annot.subset = Tool(apa.seurat.object, "GeneSLICER")[high.expressed.peaks, ]
  if (use.all.peaks == FALSE) {
    peaks.to.use = apply(annot.subset, 1, function(x) {
      ifelse(sum(x[feature.type] == "YES") >= 1, TRUE, FALSE)
    })
    peaks.to.use = names(peaks.to.use[which(peaks.to.use == TRUE)])
    high.expressed.peaks = intersect(high.expressed.peaks, peaks.to.use)
    print(paste(length(high.expressed.peaks), "expressed peaks in feature types", toString(feature.type)))
  } else {
    print(paste(length(high.expressed.peaks), "expressed peaks across all feature types"))
  }

  ## Identifiy genes with more than one polyA detected as expressed
  gene.names = Tool(apa.seurat.object, "GeneSLICER")[high.expressed.peaks, ]$Gene_name
  gene.table = table(gene.names)
  multi.genes = gene.table[gene.table > 1]
  print(paste(length(multi.genes), "genes detected with multiple polyA expressed"))
  multi.gene.names = names(multi.genes)
  multi.gene.names = multi.gene.names[which(multi.gene.names != "")]

  ## Build a matrix of genes corresponding to summed polyA counts
  print("Collapsing gene polyA expression...")
  multi.peaks = high.expressed.peaks[which(gene.names %in% multi.gene.names)]
  multi.peak.genes = Tool(apa.seurat.object, "GeneSLICER")[multi.peaks, ]$Gene_name
  sub.matrix = GetAssayData(apa.seurat.object)[multi.peaks, c(cells.fg, cells.bg)]
  collapsed.mat = pbapply::pblapply(multi.gene.names, function(thisGene){
    these.peak.indicies = which(multi.peak.genes == thisGene)
    summed.expression = apply(sub.matrix[these.peak.indicies, ], 2, function(x) sum(x))
    return(summed.expression)
  })
  collapsed.mat = do.call(rbind, collapsed.mat)
  rownames(collapsed.mat) = multi.gene.names
  colnames(collapsed.mat) = colnames(sub.matrix)
  collapsed.mat = as(as.matrix(collapsed.mat), "sparseMatrix")

  ## Now calculate differential expression on the collaped matrix
  ## First calculate fold-changes for between the two populations
  all.fcs = apply(collapsed.mat[multi.gene.names, c(cells.fg, cells.bg)], 1, function(x)
    get_log2FC(x[cells.fg], x[cells.bg]))

  ## Get the APAs and corresponding genes that pass a fold-change threshold
  fcs.pass = all.fcs[which(abs(all.fcs) >= fc.thresh)]
  genes.pass = names(fcs.pass)
  print(paste(length(genes.pass), " collapsed genes pass LFC threshold"))

  clusters = Idents(apa.seurat.object)
  clusters = clusters[c(cells.fg, cells.bg)]
  clusters = as.character(clusters)
  p.values = get_pvalue_list(collapsed.mat[genes.pass, ], identities = clusters,
                             cluster1 = cluster1, cluster2 = cluster2, test.use = test.use)
  p.adj = p.adjust(p.values, method = "bonferroni", n = nrow(collapsed.mat))

  gene.name = Tool(apa.seurat.object, "GeneSLICER")[genes.pass, "Gene_name"]

  ## Calculate the percentage of cells expressing the gene
  nz.row.foreground = tabulate(collapsed.mat[, cells.fg]@i + 1, nbins = length(multi.gene.names))
  nz.prop.foreground = nz.row.foreground/length(cells.fg)
  names(nz.prop.foreground) = rownames(collapsed.mat)

  ## Also for the background population
  nz.row.background = tabulate(collapsed.mat[, cells.bg]@i + 1, nbins = length(multi.gene.names))
  nz.prop.background = nz.row.background/length(cells.bg)
  names(nz.prop.background) = rownames(collapsed.mat)

  diff.table = data.frame(Gene = genes.pass,
                          Target_pct = nz.prop.foreground[genes.pass],
                          Background_pct = nz.prop.background[genes.pass],
                          Log2FC = fcs.pass[genes.pass],
                          p_value = p.values[genes.pass],
                          P_value_adj = p.adj[genes.pass], stringsAsFactors = FALSE)

  diff.table = diff.table[order(diff.table$P_value_adj), ]
  diff.table = subset(diff.table,  P_value_adj < adj.pval.thresh)
  print(paste0(nrow(diff.table), " differentialy expressed genes identified"))

  return(diff.table)

}

################################################################################################
#'
#' Apply a differential expression test to gene-level collapsed polyA sites.
#'
#' Collapses polyA sites into a single gene-level count for each gene. Performs differential
#' gene expression analysis on the collapsed gene-level counts.
#'
#'
#' @param apa.seurat.object a polyA Seurat object
#' @param cluster1 foreground cluster
#' @param cluster2 comparison cluster
#' @param exp.thresh threshold \% of cells expressing a peak to consider in
#' @param fc.thresh threshold log2 fold-change difference for testing a peak
#' @param adj.pval.thresh adjusted P-value threshold for retaining a peak
#' @param test.use statistical test to use. Either MAST (default) or t-test.
#' @param feature.type the feature type(s) to test for. Redundant if use.all.peaks set to TRUE.
#' @param use.all.peaks whether to use all peaks regardless of annotation. FALSE by default.
#' @return a data-frame of results.
#' @examples
#' find_de_genes_collapsed(apa.seurat.object, cluster1 = "1", cluster2 = "2")
#'
find_de_genes_collapsed_v2 <- function(apa.seurat.object, cluster1, cluster2, exp.thresh = 0.05,
                          fc.thresh = 0.25, adj.pval.thresh = 1e-05, test.use = "MAST",
                          feature.type = c("UTR3", "UTR5", "exon", "intron"), use.all.peaks = FALSE) {

  if (!(test.use %in% c("MAST", "t-test"))) {
    stop("Invalid test option: please choose either MAST or t-test")
  }

  if (test.use == "MAST") {
    if (!'MAST' %in% rownames(x = installed.packages())) {
      stop("Please install MAST before using this function  (https://github.com/RGLab/MAST)")
    }
  }

  ## Pull out the set of cells for performing the comparison
  cells.fg = apa.seurat.object@cell.names[which(apa.seurat.object@ident == cluster1)]
  cells.bg = apa.seurat.object@cell.names[which(apa.seurat.object@ident == cluster2)]

  ## Determine the APA expressed above provided threshold proportion of cells
  print("Detecting expressed APA in each cluster...")
  high.expressed.peaks = get_highly_expressed_peaks(apa.seurat.object, cluster1, cluster2, threshold = exp.thresh)
  print(paste(length(high.expressed.peaks), "expressed peaks"))

  ## Filter peaks according to feature type
  annot.subset = apa.seurat.object@misc[high.expressed.peaks, ]
  if (use.all.peaks == FALSE) {
    peaks.to.use = apply(annot.subset, 1, function(x) {
      ifelse(sum(x[feature.type] == "YES") >= 1, TRUE, FALSE)
    })
    peaks.to.use = names(peaks.to.use[which(peaks.to.use == TRUE)])
    high.expressed.peaks = intersect(high.expressed.peaks, peaks.to.use)
    print(paste(length(high.expressed.peaks), "expressed peaks in feature types", toString(feature.type)))
  } else {
    print(paste(length(high.expressed.peaks), "expressed peaks across all feature types"))
  }

  ## Identifiy genes with more than one polyA detected as expressed
  gene.names = apa.seurat.object@misc[high.expressed.peaks, ]$Gene_name
  gene.table = table(gene.names)
  multi.genes = gene.table[gene.table > 1]
  print(paste(length(multi.genes), "genes detected with multiple polyA expressed"))
  multi.gene.names = names(multi.genes)
  multi.gene.names = multi.gene.names[which(multi.gene.names != "")]

  ## Build a matrix of genes corresponding to summed polyA counts
  print("Collapsing gene polyA expression...")
  multi.peaks = high.expressed.peaks[which(gene.names %in% multi.gene.names)]
  multi.peak.genes = apa.seurat.object@misc[multi.peaks, ]$Gene_name
  sub.matrix = apa.seurat.object@data[multi.peaks, c(cells.fg, cells.bg)]
  collapsed.mat = pbapply::pblapply(multi.gene.names, function(thisGene){
    these.peak.indicies = which(multi.peak.genes == thisGene)
    summed.expression = apply(sub.matrix[these.peak.indicies, ], 2, function(x) sum(x))
    return(summed.expression)
  })
  collapsed.mat = do.call(rbind, collapsed.mat)
  rownames(collapsed.mat) = multi.gene.names
  colnames(collapsed.mat) = colnames(sub.matrix)
  collapsed.mat = as(as.matrix(collapsed.mat), "sparseMatrix")

  ## Now calculate differential expression on the collaped matrix
  ## First calculate fold-changes for between the two populations
  all.fcs = apply(collapsed.mat[multi.gene.names, c(cells.fg, cells.bg)], 1, function(x)
    get_log2FC(x[cells.fg], x[cells.bg]))

  ## Get the APAs and corresponding genes that pass a fold-change threshold
  fcs.pass = all.fcs[which(abs(all.fcs) >= fc.thresh)]
  genes.pass = names(fcs.pass)
  print(paste(length(genes.pass), " collapsed genes pass LFC threshold"))

  clusters = apa.seurat.object@ident
  clusters = clusters[c(cells.fg, cells.bg)]
  clusters = as.character(clusters)
  p.values = get_pvalue_list(collapsed.mat[genes.pass, ], identities = clusters,
                             cluster1 = cluster1, cluster2 = cluster2, test.use = test.use)
  p.adj = p.adjust(p.values, method = "bonferroni", n = nrow(collapsed.mat))

  gene.name = apa.seurat.object@misc[genes.pass, "Gene_name"]

  ## Calculate the percentage of cells expressing the gene
  nz.row.foreground = tabulate(collapsed.mat[, cells.fg]@i + 1)
  nz.prop.foreground = nz.row.foreground/length(cells.fg)
  names(nz.prop.foreground) = rownames(collapsed.mat)

  ## Also for the background population
  nz.row.background = tabulate(collapsed.mat[, cells.bg]@i + 1)
  nz.prop.background = nz.row.background/length(cells.bg)
  names(nz.prop.background) = rownames(collapsed.mat)

  diff.table = data.frame(Gene = genes.pass,
                          Target_pct = nz.prop.foreground[genes.pass],
                          Background_pct = nz.prop.background[genes.pass],
                          Log2FC = fcs.pass[genes.pass],
                          p_value = p.values[genes.pass],
                          P_value_adj = p.adj[genes.pass], stringsAsFactors = FALSE)

  diff.table = diff.table[order(diff.table$P_value_adj), ]
  diff.table = subset(diff.table,  P_value_adj < adj.pval.thresh)
  print(paste0(nrow(diff.table), " differentialy expressed genes identified"))

  return(diff.table)

}

#########################################################################################
#'
#' Get a list of P-values for differential polyA testing
#'
#' Given an expression matrix, a list of cluster identities and cluster labels, calculate
#' a list of differential expression P-values for the polyA sites in the expression matrix.
#' Currently can choose from two tests: MAST or t-test.
#'
#' @param expressionData a matrix of expression data
#' @param identities cell identity labels (normally clusters)
#' @param cluster1 target cluster to test DE for
#' @param cluster2 optional comparison cluster (all non-cluster1 cells by default)
#' @param test.use statistical test to use. Currently either MAST or t-test.
#' @return a list of P-values
#' @examples
#' p.values = get_pvalue_list(expressionData, identities, "1")
#'
get_pvalue_list <- function(expressionData, identities, cluster1, cluster2=NULL, test.use = "MAST") {

  if (test.use == "MAST") {
    p.values = mast_de_test(expressionData, identities, cluster1, cluster2)
    p.values = p.values[rownames(expressionData)]
    return(p.values)
  } else if (test.use == "t-test") {

    if (length(cluster1) == 1){ # cluster identity used as input
      cells.1 = colnames(expressionData)[which(identities == cluster1)]
    } else { # cell identity used as input
      cells.1 = cluster1
    }
    if (is.null(cluster2)) {
      cells.2 = setdiff(colnames(expressionData), cells.1)
    } else {
      if (length(cluster2) == 1) { # cluster identity used as input
        cells.2 = colnames(expressionData)[which(identities == cluster2)]
      } else { # cell identity used as input
        cells.2 = cluster2
      }
    }

    if (is.null(dim(expressionData))) {
      p.values = t.test(expressionData[cells.1], expressionData[cells.2])
    } else {
      p.values = unlist(apply(expressionData, 1, function(x) t.test(x[cells.1], x[cells.2])$p.value))
    }
    return(p.values)
  } else {
    stop("Invalid test option: please choose either MAST or t-test")
  }
}

################################################################################################
#'
#' Do MAST differential expression test.
#'
#' Apply a MAST differential expression test to polyA sites. Code derived from the Seurat 'MASTDETest'
#' function - see 'MASTDETest' at https://github.com/satijalab/seurat for more details.
#'
#' @param expressionData a matrix of polyA site expression
#' @param identities cluster labels for cells
#' @param cluster1 target cluser ID
#' @param cluster2 comparison cluster ID
#' @return a list of P-values
#' @examples
#' doMASTDETest(expressionData, identities, cluster1, cluster2)
#'
mast_de_test <- function(expressionData, identities, cluster1, cluster2=NULL) {

  if (length(cluster1) == 1){ # cluster identity used as input
    cells.1 = colnames(expressionData)[which(identities == cluster1)]
  } else { # cell identity used as input
    cells.1 = cluster1
  }
  if (is.null(cluster2)) {
    cells.2 = setdiff(colnames(expressionData), cells.1)
  } else {
    if (length(cluster2) == 1) { # cluster identity used as input
      cells.2 = colnames(expressionData)[which(identities == cluster2)]
    } else { # cell identity used as input
      cells.2 = cluster2
    }
  }

  coldata <- data.frame(group = c(rep("Group1", length(cells.1)), rep("Group2", length(cells.2))),
                        row.names = c(cells.1, cells.2))

  coldata$group <- factor(x = coldata$group)
  coldata$wellKey <- rownames(x = coldata)
  latent.vars <- c("condition")
  countdata.test <- expressionData[, c(cells.1, cells.2)]
  fdat <- data.frame(rownames(x = countdata.test))
  colnames(x = fdat)[1] <- "primerid"
  rownames(x = fdat) <- fdat[, 1]
  sca <- MAST::FromMatrix(exprsArray = as.matrix(x = countdata.test),
                          cData = coldata, fData = fdat)
  cond <- factor(x = SummarizedExperiment::colData(sca)$group)
  cond <- relevel(x = cond, ref = "Group1")
  SummarizedExperiment::colData(sca)$condition <- cond
  fmla <- as.formula(object = " ~ condition")
  zlmCond <- MAST::zlm(formula = fmla, sca = sca)
  summaryCond <- MAST::summary(object = zlmCond, doLRT = "conditionGroup2")
  summaryDt <- summaryCond$datatable
  p_val <- summaryDt[which(summaryDt[, "component"] == "H"), 4]
  genes.return <- summaryDt[which(summaryDt[, "component"] == "H"), 1]
  #to.return <- p_val$`Pr(>Chisq)`
  to.return <- p_val
  #names(to.return) = genes.return$primerid
  names(to.return) = genes.return
  return(to.return)
}


############################################################
#'
#' Use the logistic regression method by Lior Patcher group to calculate differential polyA
#'
#' Use the logistic regression method by Lior Patcher group to calculate differential polyA
#'
#' @param apa.seurat.object a polyA Seurat object
#' @param cluster1 foreground cluster
#' @param cluster2 comparison cluster
#' @param exp.thresh threshold \% of cells expressing a peak to consider in
#' @param fc.thresh threshold log2 fold-change difference for testing a peak
#' @param adj.pval.thresh adjusted P-value threshold for retaining a peak
#' @param feature.type the feature type(s) to test for
#' @return a data-frame of results
#' @examples
#' res.table = lr_de_test(apa.seurat.object, cluster1 = "1", cluster2 = "2")
#'
de_polya_lr <- function(apa.seurat.object, cluster1, cluster2, exp.thresh = 0.05, fc.thresh = 0.25,
                                  adj.pval.thresh = 1e-05, feature.type = c("UTR3", "UTR5", "exon", "intron")) {

  ## Pull out the set of cells for performing the comparison
  cells.fg = apa.seurat.object@cell.names[which(apa.seurat.object@ident == cluster1)]
  cells.bg = apa.seurat.object@cell.names[which(apa.seurat.object@ident == cluster2)]

  ## Determine the APA expressed above provided threshold proportion of cells
  print("Detecting expressed APA in each cluster...")
  high.expressed.peaks = get_highly_expressed_peaks(apa.seurat.object, cluster1, cluster2, threshold = exp.thresh)
  print(paste(length(high.expressed.peaks), "expressed peaks"))

  ## Filter peaks according to feature type
  annot.subset = apa.seurat.object@misc[high.expressed.peaks, ]
  peaks.to.use = apply(annot.subset, 1, function(x) {
    ifelse(sum(x[feature.type] == "YES") >= 1, TRUE, FALSE)
  })
  peaks.to.use = names(peaks.to.use[which(peaks.to.use == TRUE)])

  high.expressed.peaks = intersect(high.expressed.peaks, peaks.to.use)
  print(paste(length(high.expressed.peaks), "expressed peaks in feature types", toString(feature.type)))

  ## Calculate fold-changes for transcripts between the two populations
  all.fcs = get_log2FC_list(apa.seurat.object, high.expressed.peaks, cluster1 = cluster1, cluster2 = cluster2)

  ## Get the APAs and corresponding genes that pass a fold-change threshold
  fcs.pass = all.fcs[which(abs(all.fcs) >= fc.thresh)]
  apas.pass = names(fcs.pass)
  print(paste(length(apas.pass), " peaks pass LFC threshold"))
  genes.fc.pass = apa.seurat.object@misc[apas.pass, "Gene_name"]

  ## Pull out the unique gene names - should be in a geneA.1, geneA.2 etc. format
  gene.names = apa.seurat.object@misc[high.expressed.peaks, "Gene_name"]

  ## Pull out genes with more than one site
  genes.multi.apa = names(table(gene.names)[which(table(gene.names) > 1)])

  ## Finally take the intersection of multi genes and genes containing peak passing FC threshold
  genes.use = intersect(genes.fc.pass, genes.multi.apa)

  ## Filter the apa sites for genes that will be tested
  all.apa.sites = high.expressed.peaks
  all.genes = apa.seurat.object@misc[all.apa.sites, "Gene_name"]
  apa.sites.test = all.apa.sites[which(all.genes %in% genes.use)]

  exp.matrix.full = Matrix::t(apa.seurat.object@raw.data[apa.sites.test, c(cells.fg, cells.bg)])

  value.labels = c(rep(1, length(cells.fg)), rep(0, length(cells.bg)))

  print("Calculating p-values...")
  p.values = rep(NA, length(genes.use))
  pb <- progress_bar$new(format = "Processing [:bar] :percent eta: :eta",
                         total = length(genes.use), clear=FALSE)
  pb$tick(0)
  for (i in seq(1:length(genes.use))) {

    thisGene = genes.use[i]
    apa.sites = all.apa.sites[which(gene.names == thisGene)]

    exp.matrix = exp.matrix.full[, apa.sites]
    exp.matrix = as.data.frame(as.matrix(exp.matrix))
    exp.matrix = cbind(value.labels, exp.matrix)

    f.labels = paste0("Gene", "_", as.character(1:length(apa.sites)))

    p.val = lr_test(f.labels, exp.matrix, value.labels, length(cells.fg), length(cells.bg), length(apa.sites))
    p.values[i] = p.val
    pb$tick()
  }

  names(p.values) = genes.use
  pval.adjusted = p.adjust(p.values, n=nrow(apa.seurat.object@data), method = "bonferroni")

  ## Build a table of genes, p-values, peaks & fold-changes
  res.table = data.frame(Gene = genes.use, P_val = as.numeric(p.values), P_adj = as.numeric(pval.adjusted))

  ## Add the set of peaks associated with the gene, and the peak with highest fold-change between clusters
  #fcs.pass = all.fcs[which(abs(all.fcs) > fc.thresh)]
  genes.fc.pass = apa.seurat.object@misc[apas.pass, "Gene_name"]

  mean.exp.all = apply(apa.seurat.object@data[apa.sites.test, c(cells.fg, cells.bg)], 1, function(x) Log2ExpMean(x))

  res.table.updated = c()
  for (i in 1:nrow(res.table)) {
    this.res = res.table[i, ]
    this.gene = as.character(this.res$Gene)
    apa.sites = apas.pass[which(genes.fc.pass == this.gene)]
    apa.fcs = fcs.pass[apa.sites]
    apa.exp = mean.exp.all[apa.sites]
    new.res = data.frame(Gene = rep(this.gene, length(apa.sites)), P_val = rep(this.res$P_val, length(apa.sites)),
                         P_adj = rep(this.res$P_adj, length(apa.sites)), Peak = apa.sites, Log2_FC = apa.fcs,
                         Log2AveExp = apa.exp, stringsAsFactors = FALSE)
    res.table.updated = rbind(res.table.updated, new.res)
  }

  ## Add the percentage of cells expressing the peak

  ## cluster 1
  nz.row.foreground = tabulate(apa.seurat.object@data[, cells.fg]@i + 1)
  nz.prop.foreground = nz.row.foreground/length(cells.fg)
  names(nz.prop.foreground) = rownames(apa.seurat.object@data)
  nz.prop.foreground.add = nz.prop.foreground[res.table.updated$Peak]
  res.table.updated$Target_pct = nz.prop.foreground.add

  ## Also for the background population
  nz.row.background = tabulate(apa.seurat.object@data[, cells.bg]@i + 1)
  nz.prop.background = nz.row.background/length(cells.bg)
  names(nz.prop.background) = rownames(apa.seurat.object@data)
  nz.prop.background.add = nz.prop.background[res.table.updated$Peak]
  res.table.updated$Background_pct = nz.prop.background.add

  ## Add genomic feature
  features.add = apa.seurat.object@misc[res.table.updated$Peak, "FeaturesCollapsed"]
  res.table.updated$GenomicFeature = features.add

  res.table.updated = subset(res.table.updated, P_adj < adj.pval.thresh)
  res.table.updated = res.table.updated[order(abs(res.table.updated$Log2_FC), decreasing = TRUE), ]

  return(res.table.updated)
}

############################################################
#'
#' logistic-regression test for differential expression
#'
#' Perform the logistic-regression test described by Lior Patcher and coleagues, which has
#' been adapted here for use on polyA sites. Builds a logistic regression
#' function using the polyA sites for a gene as the predictors. Performance
#' for correctly classifying class labels is compared to a null model to derive
#' a P-value for the gene.
#'
#' @param f.labels class (cluster) labels for cells
#' @param exp.matrix matrix of counts
#' @param value.labels values
#' @param num_cells_fg number of 'positive' class cells
#' @param num_cells_bg number of background or 'negative' class cells
#' @param num_apa number of polyA sites used as predictors
#' @return a P-value
#' @examples
#' p.val = lr_test(f.labels, exp.matrix, value.labels, num_cells_fg, num_cells_bg, num_apa)
#'
lr_test <- function(f.labels, exp.matrix, value.labels, num_cells_fg, num_cells_bg, num_apa) {

  this.fml = as.formula(paste("Cluster ~ ", paste(f.labels, collapse = "+")))
  colnames(exp.matrix) = c("Cluster", f.labels)

  model <- glm(formula = this.fml, family = binomial(link = "logit"), data = exp.matrix)
  summary(model)

  N1 = num_cells_fg
  N2 = num_cells_bg
  k = num_apa

  p_of_1 = N1/(N1+N2)
  llnull=(N1+N2)*(p_of_1*log(p_of_1) + (1-p_of_1)*log(1-p_of_1))

  fitted.results <- predict(model, newdata=exp.matrix[, 2:ncol(exp.matrix)], type='response')

  gene_score = MLmetrics::LogLoss(fitted.results, as.numeric(value.labels))

  llf=-gene_score*(N1+N2)
  llr=llf-llnull
  llr_pval = pchisq(2*llr, k, lower.tail = FALSE)

  return(llr_pval)
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
