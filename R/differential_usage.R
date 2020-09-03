

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
#' @param seed.use seed to set the randomised assignment of cells to pseudo-bulk profiles
#' @param feature.type genomic feature types to run analysis on (default: UTR3, exon)
#' @param include.annotations whether to include junction, polyA motif and stretch annotations in output (default: FALSE)
#' @param filter.pA.stretch whether to filter out peaks annotated as proximal to an A-rich region (default: FALSE)
#' @param verbose whether to print outputs (TRUE by default)
#' @param do.MAPlot make an MA plot of results (FALSE by default)
#' @param return.dexseq.res return the raw and unfiltered DEXSeq results object (FALSE by default)
#' @param ncores number of cores to run DEXSeq with 
#' @return a data-frame of results.
#' @examples
#' 
#' 
#' 
#' extdata_path <- system.file("extdata",package = "Sierra")
#' load(paste0(extdata_path,"/TIP_cell_info.RData"))
#' \dontrun{
#' peak.annotations <- read.table("TIP_merged_peak_annotations.txt", header = TRUE,sep = "\t",
#'                                       row.names = 1,stringsAsFactors = FALSE)
#' peaks.seurat <- NewPeakSeurat(peak.data = peak.counts, 
#'                              annot.info = peak.annotations, 
#'                              cell.idents = tip.populations,
#'                              tsne.coords = tip.tsne.coordinates,
#'                              min.cells = 0, min.peaks = 0)
#' 
#' res.table = DUTest(peaks.seurat, population.1 = "F-SL", population.2 = "EC1",
#'                          exp.thresh = 0.1,  feature.type = c("UTR3", "exon"))
#' }
#'
#' @export
#'
DUTest <- function(peaks.object, 
                   population.1, 
                   population.2 = NULL, 
                   exp.thresh = 0.1,
                   fc.thresh=0.25, 
                   adj.pval.thresh = 0.05, 
                   num.splits = 6, 
                   seed.use = 1,
                   feature.type = c("UTR3", "exon"), 
                   include.annotations = FALSE,
                   filter.pA.stretch = FALSE,
                   verbose = TRUE, 
                   do.MAPlot = FALSE,
                   return.dexseq.res = FALSE, 
                   ncores = 1) {

  if (class(peaks.object) == "Seurat") {
    res.table <- apply_DEXSeq_test_seurat(apa.seurat.object = peaks.object,
                                          population.1 = population.1, 
                                          population.2 = population.2,
                                          exp.thresh = exp.thresh, 
                                          fc.thresh = fc.thresh,
                                          adj.pval.thresh = adj.pval.thresh, 
                                          num.splits = num.splits,
                                          seed.use = seed.use, 
                                          feature.type = feature.type,
                                          include.annotations = include.annotations,
                                          filter.pA.stretch = filter.pA.stretch,
                                          verbose = verbose, 
                                          do.MAPlot = do.MAPlot,
                                          return.dexseq.res = return.dexseq.res, 
                                          ncores = ncores)
    return(res.table)

  } else if (class(peaks.object) == "SingleCellExperiment") {
    res.table <- apply_DEXSeq_test_sce(peaks.sce.object = peaks.object,
                                       population.1 = population.1, 
                                       population.2 = population.2,
                                       exp.thresh = exp.thresh, 
                                       fc.thresh = fc.thresh,
                                       adj.pval.thresh = adj.pval.thresh,
                                       num.splits = num.splits,
                                       seed.use = seed.use, 
                                       feature.type = feature.type,
                                       include.annotations = include.annotations,
                                       filter.pA.stretch = filter.pA.stretch,
                                       verbose = verbose, 
                                       do.MAPlot = do.MAPlot,
                                       return.dexseq.res = return.dexseq.res, 
                                       ncores = ncores)
  } else{
    stop("Invalid data object provided.")
  }


}

#######################################################################
#'
#' Detect shifts in 3'UTR length usage between cell populations
#' 
#' Detect global shifts in 3'UTR length usage between defined cell populations.
#' Firsts applies the DUTest function to detect differential usage (DU) peaks on 3'UTRs, 
#' after filtering out peaks annotated as proximal to A-rich regions. Identifies peaks
#' on the same 3'UTR as each DU peak, and determines a position of the DU peak on the
#' 3'UTR relative to the terminating exon. Returns a table of DU results, with the location 
#' of each peak relative to the total number of peaks on the corresponding 3'UTR. Results 
#' table can be input to the PlotUTRLengthShift function to visualise the results, 
#' and evaluate global shifts. 
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
#' @param seed.use seed to set the randomised assignment of cells to pseudo-bulk profiles
#' @param verbose whether to print outputs (TRUE by default)
#' @param do.MAPlot make an MA plot of results (FALSE by default)
#' @param ncores Number of cores for multithreading
#' @return a data-frame of results.
#' 
#' 
#' 
#' 
#' @importFrom magrittr "%>%"
#' 
#' @export
#' 
DetectUTRLengthShift <- function(peaks.object, 
                                   gtf_gr,
                                   gtf_TxDb,
                                   population.1,
                                   population.2 = NULL,
                                   exp.thresh = 0.1,
                                   fc.thresh = 0.25,
                                   adj.pval.thresh = 0.05,
                                   num.splits = 6,
                                   seed.use = 1,
                                   verbose = TRUE,
                                   do.MAPlot = FALSE,
                                   ncores = 1) {
  
  res.table = DUTest(peaks.object = peaks.object, 
                     population.1 = population.1, 
                     population.2 = population.2, 
                     exp.thresh = exp.thresh,
                     feature.type = c("UTR3"), 
                     filter.pA.stretch = TRUE,
                     fc.thresh = fc.thresh, 
                     adj.pval.thresh = adj.pval.thresh,
                     num.splits = num.splits,
                     seed.use = seed.use,
                     verbose = verbose,
                     do.MAPlot = do.MAPlot, 
                     ncores = ncores)
  
  if (nrow(res.table) == 0) {
    print("No DU peaks identified")
    return(NULL)
  }
  
  if (verbose) print("Detecting shifts in 3'UTR length usage")
  
  ## retrieve the annotation information
  annot.df <- Tool(peaks.object, "Sierra")
  
  ## format differentially used peaks from GRanges
  all.peaks = rownames(res.table)
  res.table$peak_ID <- all.peaks
  strand = sub(".*:.*:.*-.*:(.*)", "\\1", all.peaks)
  strand = plyr::mapvalues(x = strand, from = c("1", "-1"), to = c("+", "-"))
  peak.remainder = sub(".*:(.*:.*-.*):.*", "\\1", all.peaks)
  peaks.use = paste0(peak.remainder, ":", strand)
  res.table$granges_peaks <- peaks.use
  res.table %>% dplyr::distinct(granges_peaks, .keep_all = TRUE) -> res.table
  rownames(res.table) <- res.table$peak_ID
  du.peaks.gr <- GenomicRanges::GRanges(peaks.use)
  
  ## format expressed peaks using GRanges
  ## get peaks that are expressed
  peaks.expressed <- GetExpressedPeaks(peaks.object, population.1 = population.1,
                                       population.2 = population.2, threshold=0.1)
  
  ## cross-reference expressed peaks with peaks not tagged as a-rich
  expressed.arich.annot <- annot.df[peaks.expressed, "pA_stretch"]
  peaks.non.arich <- rownames(subset(annot.df, pA_stretch == FALSE))
  peaks.expressed <- intersect(peaks.expressed, peaks.non.arich)
  
  genes.expressed <- sub("(.*):.*:.*-.*:.*", "\\1", peaks.expressed)
  peaks.use.idx <- which(genes.expressed %in% res.table$gene_name)
  peaks.expressed <- peaks.expressed[peaks.use.idx]
  
  strand = sub(".*:.*:.*-.*:(.*)", "\\1", peaks.expressed)
  strand = plyr::mapvalues(x = strand, from = c("1", "-1"), to = c("+", "-"))
  peak.remainder = sub(".*:(.*:.*-.*):.*", "\\1", peaks.expressed)
  peaks.expressed.granges = paste0(peak.remainder, ":", strand)
  expressed.peaks.gr <- GenomicRanges::GRanges(peaks.expressed.granges)
  
  ## make a table mapping peaks to granges peaks for later
  granges_peaks_mapping_table <- data.frame(PeakID = peaks.expressed,
                                            row.names = peaks.expressed.granges, 
                                            stringsAsFactors = FALSE)
  
  utr3.ref <- GenomicFeatures::threeUTRsByTranscript(gtf_TxDb)
  utr3.ref <- unlist(utr3.ref)
  
  
  all_UTR_3_hits <- GenomicRanges::findOverlaps(expressed.peaks.gr , utr3.ref, type = "any")
  utr3.mappings <- as.data.frame(all_UTR_3_hits)
  
  query.hit.df <- as.data.frame(expressed.peaks.gr[utr3.mappings$queryHits, ])
  subject.hit.df <- as.data.frame(utr3.ref[utr3.mappings$subjectHits, ],
                                  row.names = as.character(1:nrow(utr3.mappings)))
  
  query.hit.df %>% dplyr::mutate(granges_peak = paste0(seqnames,":",start,"-",end,":",strand)) ->
    query.hit.df
  
  utr3.mappings$exon_name <- subject.hit.df$exon_name
  utr3.mappings$granges_peak <- query.hit.df$granges_peak
  peak.ids <- granges_peaks_mapping_table[as.character(query.hit.df$granges_peak), 'PeakID']
  utr3.mappings$peakID <- peak.ids
  utr3.mappings %>% dplyr::mutate(Gene_name = sub("(.*):.*:.*-.*:.*", "\\1", peakID)) -> utr3.mappings
  utr3.mappings %>% dplyr::mutate(Start = sub(".*:(.*)-.*:.*", "\\1", granges_peak), 
                                  End = sub(".*:.*-(.*):.*", "\\1", granges_peak)) -> utr3.mappings
  
  #### Here calculate a 'relative peak location' to the start of the UTR ####
  
  ## Go through the peaks in a strand-specific manner
  strand <- sub(".*:.*:.*-.*:(.*)", "\\1", rownames(res.table))
  res.table$Strand <- strand
  res.table.pos.strand <- subset(res.table, Strand == "1")
  res.table.neg.strand <- subset(res.table, Strand == "-1")
  
  ## Upregulated positive-strand peaks 
  peaks.res.pos.up <- subset(res.table.pos.strand, Log2_fold_change > 0)
  locations.res.table.pos.up <- make_utr3_peak_location_table(peaks.object, peaks.res.pos.up, "1", 
                                                              peaks.expressed, utr3.mappings)
  
  ## Downregulated positive-strand peaks 
  peaks.res.pos.down <- subset(res.table.pos.strand, Log2_fold_change < 0)
  locations.res.table.pos.down <- make_utr3_peak_location_table(peaks.object, peaks.res.pos.down, "1", 
                                                                peaks.expressed, utr3.mappings)
  
  ## Upregulated negative-strand peaks 
  peaks.res.neg.up <- subset(res.table.neg.strand, Log2_fold_change > 0)
  locations.res.table.neg.up <- make_utr3_peak_location_table(peaks.object, peaks.res.neg.up, "-1", 
                                                              peaks.expressed, utr3.mappings)
  
  ## Downregulated negative-strand peaks 
  peaks.res.neg.down <- subset(res.table.neg.strand, Log2_fold_change < 0)
  locations.res.table.neg.down <- make_utr3_peak_location_table(peaks.object, peaks.res.neg.down, "-1", 
                                                                peaks.expressed, utr3.mappings)
  if (!is.null(locations.res.table.pos.up)) {
    locations.res.table.pos.up$FC_direction <- rep("Up", nrow(locations.res.table.pos.up))
    } else {locations.res.table.pos.up <- c()}
  
  if (!is.null(locations.res.table.pos.down)) {
    locations.res.table.pos.down$FC_direction <- rep("Down", nrow(locations.res.table.pos.down))
    } else {locations.res.table.pos.down <- c()}
  
  if (!is.null(locations.res.table.neg.up)) {
    locations.res.table.neg.up$FC_direction <- rep("Up", nrow(locations.res.table.neg.up))
  } else {locations.res.table.neg.up <- c()}
  
  if (!is.null(locations.res.table.neg.down)) {
    locations.res.table.neg.down$FC_direction <- rep("Down", nrow(locations.res.table.neg.down))
  } else {locations.res.table.neg.down <- c()}
   
  tables.combined <- do.call(rbind, list(locations.res.table.pos.up, locations.res.table.pos.down,
                                         locations.res.table.neg.up, locations.res.table.neg.down))
  
  ## Add output information from DEXSEq back to the table
  dexseq.res <- res.table[rownames(tables.combined), ]
  info.add <- c("genomic_feature(s)", "population1_pct", "population2_pct", "pvalue", "padj", "Log2_fold_change")
  tables.combined <- cbind(dexseq.res[, info.add], tables.combined)
  
  return(tables.combined)
}



#######################################################################
#'
#' Find alternative 3' end usage between two single-cell populations
#'
#' Wrapper function to DUTest for detecting differential 3' end use. First applies DUTest to
#' test for differential usage between 3'UTRs. For DU 3'UTR peaks, evaluates whether the DU peaks
#' fall in different 3'UTRs.
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
#' @param seed.use seed to set the randomised assignment of cells to pseudo-bulk profiles
#' @param verbose whether to print outputs (TRUE by default)
#' @param do.MAPlot make an MA plot of results (FALSE by default)
#' @param ncores Number of cores for multithreading
#' @return a data-frame of results.
#' @examples
#' \dontrun{
#'      DetectAEU(apa.seurat.object, population.1 = "1", population.2 = "2")
#'  }
#' @export
#' 
#' @importFrom magrittr "%>%"
#'
DetectAEU <- function(peaks.object, gtf_gr, gtf_TxDb, population.1, population.2 = NULL, exp.thresh = 0.1,
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
  peaks.expressed <- GetExpressedPeaks(peaks.object, population.1 = population.1,
                                       population.2 = population.2, threshold=exp.thresh)

  genes.expressed <- sub("(.*):.*:.*-.*:.*", "\\1", peaks.expressed)
  peaks.use.idx <- which(genes.expressed %in% as.character(res.table$gene_name))
  peaks.expressed <- peaks.expressed[peaks.use.idx]

  strand = sub(".*:.*:.*-.*:(.*)", "\\1", peaks.expressed)
  strand = plyr::mapvalues(x = strand, from = c("1", "-1"), to = c("+", "-"))
  peak.remainder = sub(".*:(.*:.*-.*):.*", "\\1", peaks.expressed)
  peaks.expressed.granges = paste0(peak.remainder, ":", strand)

  ## filter out duplicate coordinates mapping to different genes
  duplicate.peaks <- names(table(peaks.expressed.granges))[which( table(peaks.expressed.granges) > 1)]
  if (length(duplicate.peaks) > 0) {
    peaks.remove.idx <- which(peaks.expressed.granges %in% duplicate.peaks)
    peaks.expressed.granges <- peaks.expressed.granges[-peaks.remove.idx]
    peaks.expressed <- peaks.expressed[-peaks.remove.idx]
    res.table <- res.table[intersect(rownames(res.table), peaks.expressed)]
  }

  expressed.peaks.gr <- GenomicRanges::GRanges(peaks.expressed.granges)

  ## make a table mapping peaks to granges peaks for later
  granges_peaks_mapping_table <- data.frame(PeakID = peaks.expressed,
                                            row.names = peaks.expressed.granges, 
                                            stringsAsFactors = FALSE)

  utr3.ref <- GenomicFeatures::threeUTRsByTranscript(gtf_TxDb)
  utr3.ref <- unlist(utr3.ref)

  all_UTR_3_hits <- GenomicRanges::findOverlaps(expressed.peaks.gr , utr3.ref, type = "any")
  utr3.mappings <- as.data.frame(all_UTR_3_hits)

  query.hit.df <- as.data.frame(expressed.peaks.gr[utr3.mappings$queryHits, ])
  subject.hit.df <- as.data.frame(utr3.ref[utr3.mappings$subjectHits, ],
                                  row.names = as.character(1:nrow(utr3.mappings)))

  query.hit.df %>% dplyr::mutate(granges_peak = paste0(seqnames,":",start,"-",end,":",strand)) ->
    query.hit.df

  utr3.mappings$exon_name <- subject.hit.df$exon_name
  utr3.mappings$granges_peak <- query.hit.df$granges_peak
  peak.ids <- granges_peaks_mapping_table[as.character(query.hit.df$granges_peak), 'PeakID']
  utr3.mappings$peakID <- peak.ids

  res.table$peak_name <- rownames(res.table)

  ## For each DU peak, identify all remaining expressed peaks falling on 3'UTRs within the
  ## relevant gene. If the 3'UTR ID/s of the DU peak are different to the remaining peaks,
  ## mark the DU peak as differential transcript usage.
  diff.transcript.check.values <- apply(as.matrix(res.table), 1, function(x) {
    diff.site <- x["peak_name"]

    this.gene <- x["gene_name"]
    all.sites <- SelectGenePeaks(peaks.object, gene = this.gene, feature.type = "UTR3")
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

  all_transcript_hits <- GenomicRanges::findOverlaps(peaks.gr , transcripts.ref, type = "any")
  transcript.mappings <- as.data.frame(all_transcript_hits)

  query.hit.df <- as.data.frame(peaks.gr[transcript.mappings$queryHits, ])
  subject.hit.df <- as.data.frame(transcripts.ref[transcript.mappings$subjectHits, ],
                                  row.names = as.character(1:nrow(transcript.mappings)))

  query.hit.df %>% dplyr::mutate(granges_peak = paste0(seqnames,":",start,"-",end,":",strand)) ->
    query.hit.df

  transcript.mappings$transcript_name <- subject.hit.df$tx_name
  transcript.mappings$granges_peak <- query.hit.df$granges_peak
  peak.ids <- granges_peaks_mapping_table[as.character(query.hit.df$granges_peak), 'PeakID']
  transcript.mappings$peakID <- peak.ids

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
#' @param seed.use seed
#' @param feature.type genomic feature types to run analysis on (default: all)
#' @param include.annotations whether to include junction, polyA motif and stretch annotations in output (default: FALSE)
#' @param filter.pA.stretch whether to filter out peaks annotated as proximal to an A-rich region (default: FALSE)
#' @param verbose whether to print outputs (TRUE by default)
#' @param do.MAPlot make an MA plot of results (FALSE by default)
#' @param return.dexseq.res return the raw and unfiltered DEXSeq results object (FALSE by default)
#' @param ncores Number of cores to use for multithreading
#' @return a data-frame of results.
#' @examples
#' 
#' \dontrun{
#'    apply_DEXSeq_test(apa.seurat.object, population.1 = "1", population.2 = "2")
#'  }
#'
apply_DEXSeq_test_seurat <- function(apa.seurat.object, 
                                     population.1, 
                                     population.2 = NULL, 
                                     exp.thresh = 0.1,
                                     fc.thresh=0.25, 
                                     adj.pval.thresh = 0.05, 
                                     num.splits = 6, 
                                     seed.use = 1,
                                     feature.type = c("UTR3", "UTR5", "exon", "intron"),
                                     include.annotations = FALSE,
                                     filter.pA.stretch = FALSE, 
                                     verbose = TRUE, 
                                     do.MAPlot = FALSE,
                                     return.dexseq.res = FALSE, 
                                     ncores = 1) {

  if (!'DEXSeq' %in% rownames(x = installed.packages())) {
    stop("Please install DEXSeq before using this function
         (http://bioconductor.org/packages/release/bioc/html/DEXSeq.html)")
  }

  ## reduce counts in a cluster to num.splits cells for genes with > 1 peak
  high.expressed.peaks <- GetExpressedPeaks(apa.seurat.object, population.1, population.2, threshold = exp.thresh)
  length(high.expressed.peaks)


  ## Filter peaks according to feature type
  annot.subset <- Tool(apa.seurat.object, "Sierra")[high.expressed.peaks, ]
  peaks.to.use <- apply(annot.subset, 1, function(x) {
    ifelse(sum(x[feature.type] == "YES") >= 1, TRUE, FALSE)
  })
  peaks.to.use <- names(peaks.to.use[which(peaks.to.use == TRUE)])
  high.expressed.peaks <- intersect(high.expressed.peaks, peaks.to.use)
  if (verbose) print(paste(length(high.expressed.peaks), "expressed peaks in feature types", toString(feature.type)))

  ## Check if A-rich peaks are to be filtered out
  if (filter.pA.stretch) {
    if (is.null(Tool(apa.seurat.object, "Sierra")$pA_stretch)) {
      stop("pA_stretch not in annotation data: please run nnotate_gr_from_gtf with
           an input genome to provide required annotation.")
    } else{
      annot.subset <- Tool(apa.seurat.object, "Sierra")[high.expressed.peaks, ]
      peaks.non.arich <- rownames(subset(annot.subset, pA_stretch == FALSE))
      high.expressed.peaks <- intersect(high.expressed.peaks, peaks.non.arich)
      if (verbose) print(paste(length(high.expressed.peaks), "peaks after filtering out A-rich annotations"))
    }
  }

  annotations.highly.expressed <- Tool(apa.seurat.object, "Sierra")[high.expressed.peaks, ]
  
  ## Make sure that gene annotations are not empty
  annotations.highly.expressed <- subset(annotations.highly.expressed, Gene_name != "")
  high.expressed.peaks <- rownames(annotations.highly.expressed)
  gene.names <- annotations.highly.expressed[, "Gene_name"]

  ## Identifiy genes with more than one transcript detected as expressed
  gene.table <- table(gene.names)
  multi.genes <- gene.table[gene.table > 1]
  if (verbose) print(paste(length(multi.genes), "genes detected with multiple peak sites expressed"))
  multi.gene.names <- names(multi.genes)

  peaks.use <- high.expressed.peaks[which(gene.names %in% multi.gene.names)]
  if (verbose) print(paste(length(peaks.use), "individual peak sites to test"))

  ## make pseudo-bulk profiles out of cells
  ## set a seed to allow replication of results
  set.seed(seed.use)
  if (length(population.1) == 1) {
    cells.1 <- names(Seurat::Idents(apa.seurat.object))[which(Seurat::Idents(apa.seurat.object) == population.1)]
  } else{
      cells.1 <- population.1
    }

  cells.1 = sample(cells.1)
  cell.sets1 <- split(cells.1, sort(1:length(cells.1)%%num.splits))

  ## create a profile set for first cluster
  profile.set1 = matrix(, nrow = length(peaks.use), ncol = length(cell.sets1))
  for (i in 1:length(cell.sets1)) {
    this.set <- cell.sets1[[i]]
    sub.matrix <- Seurat::GetAssayData(apa.seurat.object, slot = "counts", assay = "RNA")[peaks.use, this.set]
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
      cells.2 <- names(Seurat::Idents(apa.seurat.object))[which(Seurat::Idents(apa.seurat.object) == population.2)]
    } else {
      cells.2 <- population.2
    }
  }

  cells.2 = sample(cells.2)
  cell.sets2 <- split(cells.2, sort(1:length(cells.2)%%num.splits))

  profile.set2 = matrix(, nrow = length(peaks.use), ncol = length(cell.sets2))
  for (i in 1:length(cell.sets2)) {
    this.set <- cell.sets2[[i]]
    sub.matrix <- Seurat::GetAssayData(apa.seurat.object, slot = "counts", assay = "RNA")[peaks.use, this.set]
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

  dexseq.feature.table <- Tool(apa.seurat.object, "Sierra")[, c("Gene_name", "Gene_part", "Peak_number")]
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
  feature.type <- Tool(apa.seurat.object, "Sierra")[rownames(dxrSig_subset), c("FeaturesCollapsed")]
  dxrSig_subset$feature_type = feature.type
  
  if (include.annotations) {
    junction.annot <- Tool(apa.seurat.object, "Sierra")[rownames(dxrSig_subset), c("Junctions", "pA_motif", "pA_stretch")]
    dxrSig_subset <- cbind(dxrSig_subset, junction.annot)
    
    dxrSig_subset <- dxrSig_subset[, c("groupID", "feature_type", "Junctions",  "pA_motif", "pA_stretch",
                                       "population1_pct", "population2_pct",
                                       "pvalue", "padj", "log2fold_target_comparison")]
    colnames(dxrSig_subset) <- c("gene_name", "genomic_feature(s)", "Junctions",  "pA_motif", "pA_stretch",
                                 "population1_pct", "population2_pct", "pvalue", "padj",  "Log2_fold_change")
  } else {
    dxrSig_subset <- dxrSig_subset[, c("groupID", "feature_type", "population1_pct", "population2_pct",
                                       "pvalue", "padj", "log2fold_target_comparison")]
    colnames(dxrSig_subset) <- c("gene_name", "genomic_feature(s)", "population1_pct",
                                 "population2_pct", "pvalue", "padj",  "Log2_fold_change") 
  }

  dxrSig_subset <- dxrSig_subset[order(dxrSig_subset$padj, decreasing = FALSE), ]

  return(dxrSig_subset)
}


#######################################################################
#'
#' Apply DEXSeq to detect differential peak usage to a Single-Cell Experiment object
#'
#' Apply DEXSeq to detect differential peak usage been select populations. Works by building
#' a 'pseudo-bulk' profile of cell populations by aggregating counts from individual cells
#' into a smaller number of profiles, defined by num.splits.
#'
#' @param peaks.sce.object SCE object of peaks
#' @param population.1 a target population of cells (can be an ID/cluster label or a set of cell barcode IDs)
#' @param population.2 comparison population of cells. If NULL (default), uses all non-population.1 cells
#' @param exp.thresh minimum percent expression threshold (for a population of cells) to include a peak
#' @param fc.thresh threshold for log2 fold-change difference for returned results
#' @param adj.pval.thresh threshold for adjusted P-value for returned results
#' @param num.splits the number of pseudo-bulk profiles to create per identity class (default: 6)
#' @param seed.use seed use
#' @param feature.type genomic feature types to run analysis on (degault: all)
#' @param include.annotations whether to include junction, polyA motif and stretch annotations in output (default: FALSE)
#' @param filter.pA.stretch whether to filter out peaks annotated as proximal to an A-rich region (default: FALSE)
#' @param verbose whether to print outputs (TRUE by default)
#' @param do.MAPlot make an MA plot of results (FALSE by default)
#' @param return.dexseq.res return the raw and unfiltered DEXSeq results object (FALSE by default)
#' @param ncores Number of cores to use for multithreading
#' @return a data-frame of results.
#' @examples
#' 
#' \dontrun{
#' apply_DEXSeq_test_sce(apa.seurat.object, population.1 = "1", population.2 = "2")
#' }
#'
apply_DEXSeq_test_sce <- function(peaks.sce.object, 
                                  population.1, 
                                  population.2 = NULL, 
                                  exp.thresh = 0.1,
                                  fc.thresh=0.25, 
                                  adj.pval.thresh = 0.05, 
                                  num.splits = 6, 
                                  seed.use = 1,
                                  feature.type = c("UTR3", "UTR5", "exon", "intron"),
                                  include.annotations = FALSE,
                                  filter.pA.stretch = FALSE, 
                                  verbose = TRUE,
                                  do.MAPlot = FALSE, 
                                  return.dexseq.res = FALSE, 
                                  ncores = 1) {

  if (!'DEXSeq' %in% rownames(x = installed.packages())) {
    stop("Please install DEXSeq before using this function
         (http://bioconductor.org/packages/release/bioc/html/DEXSeq.html)")
  }

  ## reduce counts in a cluster to num.splits cells for genes with > 1 peak
  high.expressed.peaks <- GetExpressedPeaks(peaks.sce.object, population.1, population.2, threshold = exp.thresh)
  length(high.expressed.peaks)

  ## Filter peaks according to feature type
  annot.subset <- peaks.sce.object@metadata$Sierra[high.expressed.peaks, ]
  peaks.to.use <- apply(annot.subset, 1, function(x) {
    ifelse(sum(x[feature.type] == "YES") >= 1, TRUE, FALSE)
  })
  peaks.to.use <- names(peaks.to.use[which(peaks.to.use == TRUE)])
  high.expressed.peaks <- intersect(high.expressed.peaks, peaks.to.use)
  if (verbose) print(paste(length(high.expressed.peaks), "expressed peaks in feature types", toString(feature.type)))

  ## Check if A-rich peaks are to be filtered out
  if (filter.pA.stretch) {
    if (is.null(peaks.sce.object@metadata$Sierra$pA_stretch)) {
      stop("pA_stretch not in annotation data: please run nnotate_gr_from_gtf with
           an input genome to provide required annotation.")
    } else{
      annot.subset <- peaks.sce.object@metadata$Sierra[high.expressed.peaks, ]
      peaks.non.arich <- rownames(subset(annot.subset, pA_stretch == FALSE))
      high.expressed.peaks <- intersect(high.expressed.peaks, peaks.non.arich)
      if (verbose) print(paste(length(high.expressed.peaks), "peaks after filtering out A-rich annotations"))
    }
  }

  annotations.highly.expressed <- Tool(apa.seurat.object, "Sierra")[high.expressed.peaks, ]
  
  ## Make sure that gene annotations are not empty
  annotations.highly.expressed <- subset(annotations.highly.expressed, Gene_name != "")
  high.expressed.peaks <- rownames(annotations.highly.expressed)
  gene.names <- annotations.highly.expressed[, "Gene_name"]

  ## Identifiy genes with more than one transcript detected as expressed
  gene.table <- table(gene.names)
  multi.genes <- gene.table[gene.table > 1]
  if (verbose) print(paste(length(multi.genes), "genes detected with multiple peak sites expressed"))
  multi.gene.names <- names(multi.genes)

  peaks.use <- high.expressed.peaks[which(gene.names %in% multi.gene.names)]
  if (verbose) print(paste(length(peaks.use), "individual peak sites to test"))

  ## make pseudo-bulk profiles out of cells
  ## set a seed to allow replication of results
  set.seed(seed.use)
  if (length(population.1) == 1) {
    cells.1 <- names(colData(peaks.sce.object)$CellIdent)[which(colData(peaks.sce.object)$CellIdent == population.1)]
  } else{
    cells.1 <- population.1
  }

  cells.1 = sample(cells.1)
  cell.sets1 <- split(cells.1, sort(1:length(cells.1)%%num.splits))

  ## create a profile set for first cluster
  profile.set1 = matrix(, nrow = length(peaks.use), ncol = length(cell.sets1))
  for (i in 1:length(cell.sets1)) {
    this.set <- cell.sets1[[i]]
    sub.matrix <- SingleCellExperiment::counts(peaks.sce.object)[peaks.use, this.set]
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
    cells.2 <- setdiff(colnames(peaks.sce.object), cells.1)
  } else {
    if (length(population.2) == 1) {
      cells.2 <- names(colData(peaks.sce.object)$CellIdent)[which(colData(peaks.sce.object)$CellIdent == population.2)]
    } else {
      cells.2 <- population.2
    }
  }

  cells.2 = sample(cells.2)
  cell.sets2 <- split(cells.2, sort(1:length(cells.2)%%num.splits))

  profile.set2 = matrix(, nrow = length(peaks.use), ncol = length(cell.sets2))
  for (i in 1:length(cell.sets2)) {
    this.set <- cell.sets2[[i]]
    sub.matrix <- SingleCellExperiment::counts(peaks.sce.object)[peaks.use, this.set]
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

  dexseq.feature.table <- peaks.sce.object@metadata$Sierra[, c("Gene_name", "Gene_part", "Peak_number")]
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

  population.1.pct <- get_percent_expression(peaks.sce.object, population.1, remainder=FALSE, geneSet=rownames(dxrSig_subset))
  if (is.null(population.2)) {
    population.2.pct <- get_percent_expression(peaks.sce.object, population.1, remainder=TRUE, geneSet=rownames(dxrSig_subset))
    population.2 <- "Remainder"
  } else {
    population.2.pct <- get_percent_expression(peaks.sce.object, population.2, remainder=FALSE, geneSet=rownames(dxrSig_subset))
  }

  dxrSig_subset$population1_pct <- population.1.pct
  dxrSig_subset$population2_pct <- population.2.pct

  ## Add Genomic feature type
  feature.type <- peaks.sce.object@metadata$Sierra[rownames(dxrSig_subset), c("FeaturesCollapsed")]
  dxrSig_subset$feature_type = feature.type

  if (include.annotations) {
    junction.annot <- Tool(apa.seurat.object, "Sierra")[rownames(dxrSig_subset), c("Junctions", "pA_motif", "pA_stretch")]
    dxrSig_subset <- cbind(dxrSig_subset, junction.annot)
    
    dxrSig_subset <- dxrSig_subset[, c("groupID", "feature_type", "Junctions",  "pA_motif", "pA_stretch",
                                       "population1_pct", "population2_pct",
                                       "pvalue", "padj", "log2fold_target_comparison")]
    colnames(dxrSig_subset) <- c("gene_name", "genomic_feature(s)", "Junctions",  "pA_motif", "pA_stretch",
                                 "population1_pct", "population2_pct", "pvalue", "padj",  "Log2_fold_change")
  } else {
    dxrSig_subset <- dxrSig_subset[, c("groupID", "feature_type", "population1_pct", "population2_pct",
                                       "pvalue", "padj", "log2fold_target_comparison")]
    colnames(dxrSig_subset) <- c("gene_name", "genomic_feature(s)", "population1_pct",
                                 "population2_pct", "pvalue", "padj",  "Log2_fold_change") 
  }

  dxrSig_subset <- dxrSig_subset[order(dxrSig_subset$padj, decreasing = FALSE), ]

  return(dxrSig_subset)
}


### Given a results table from DUTest, order the peaks according to position
make_utr3_peak_location_table <- function(peaks.object, res.table, strand, peaks.expressed, utr3.mappings) {
  
  if (strand == "1" | strand == "+") {
    sort.decrease <- FALSE
  } else if (strand == "-1" | strand == "-") {
    sort.decrease <- TRUE
  } else{
    stop("Invalid strand specification - options are 1, +, -1 or -")
  }
  
  locations.res.table <- c()
  
  for (this.gene in unique(res.table$gene_name)) {
    all.sites <- SelectGenePeaks(peaks.object = peaks.object, gene = this.gene, feature.type = "UTR3")
    sites.expressed <- intersect(all.sites, peaks.expressed)
    
    ## Pull out the differentiall used sites/s (can be more than one)
    diff.site.set <- rownames(subset(res.table, gene_name == this.gene))
    
    for (diff.site in diff.site.set) {
      exons.use <- unique(subset(utr3.mappings, peakID %in% diff.site)$exon_name)
      sites.use <- unique(subset(utr3.mappings, exon_name %in% exons.use)$peakID)
      sites.expressed <- intersect(sites.expressed, sites.use)
      num.sites <- length(sites.expressed)
      
      if (length(sites.expressed) > 1) {
        ## Order the expressed sites and record the relative location of 
        ## the differentially used sites.
        start.sites <- as.numeric(sub(".*:.*:(.*)-.*:.*", "\\1", sites.expressed))
        names(start.sites) <- sites.expressed
        start.sites <- sort(start.sites, decreasing = sort.decrease)
        site.diff <- which(names(start.sites) %in% diff.site)
        
        if (length(site.diff) == 1 & length(num.sites) == 1 & length(diff.site) == 1) {
          this.res <- data.frame(SiteLocation = site.diff,
                                 NumSites = num.sites,
                                 row.names = diff.site,
                                 stringsAsFactors = FALSE)
          locations.res.table <- rbind(locations.res.table, this.res)
        }
        
      }
    }
  }
  
  return(locations.res.table)
}

############################################################
#'
#' Identify peaks expressed within a certain percentage of cells
#'
#' Selects peaks that are considered expressed above some provided criteria within a target or
#' background cluster. Considers peaks expressed in some x\% of cells to be highly expressed. Returns the
#' union of peaks identified from the target and background cluster
#'
#' @param peaks.object the peaks object either Seurat of SingleCellExperiment class.
#' @param population.1 target cluster
#' @param population.2 background cluster. If NULL (deafult) all non-target cells
#' @param threshold percentage threshold of detected (non-zero) expression for including a peak
#' @return an array of peak (or gene) names
#' @examples
#' 
#' \dontrun{
#'     get_highly_expressed_peaks(seurat.object, "1")
#'     get_highly_expressed_peaks(seurat.object, cluster1 = "1", cluster2 = "2")
#'  }
#' @export
#'
GetExpressedPeaks <- function(peaks.object, population.1, population.2=NULL, threshold=0.05) {

  if (class(peaks.object) == "Seurat") {
    expressed.peaks <- get_expressed_peaks_seurat(peaks.seurat.object = peaks.object,
                                                  population.1 = population.1,
                                                  population.2 = population.2,
                                                  threshold=threshold)
  } else if (class(peaks.object) == "SingleCellExperiment") {
    expressed.peaks <- get_expressed_peaks_sce(peaks.sce.object = peaks.object,
                                                  population.1 = population.1,
                                                  population.2 = population.2,
                                                  threshold=threshold)
  }

  return(expressed.peaks)
}

############################################################
#'
#' Identify highly expressed peaks
#'
#' Selects peaks that are considered expressed above some provided criteria within a target or
#' background cluster. Considers peaks expressed in some x\% of cells to be highly expressed. Returns the
#' union of peaks identified from the target and background cluster
#'
#' @param peaks.seurat.object the peak-count Seurat object
#' @param population.1 target cluster
#' @param population.2 background cluster. If NULL (deafult) all non-target cells
#' @param threshold percentage threshold of detected (non-zero) expression for including a peak
#' @return an array of peak (or gene) names
#' @examples
#' \dontrun{
#'   get_highly_expressed_peaks(seurat.object, "1")
#'   get_highly_expressed_peaks(seurat.object, population.1 = "1", population.2 = "2")
#' }
get_expressed_peaks_seurat <- function(peaks.seurat.object, population.1, population.2=NULL, threshold=0.05) {

  if (length(population.1) == 1){ # cluster identity used as input
    foreground.set = names(Seurat::Idents(peaks.seurat.object)[Seurat::Idents(peaks.seurat.object)==population.1])
  } else { # cell identity used as input
    foreground.set = population.1
  }
  if (is.null(population.2)) {
    remainder.set = names(Seurat::Idents(peaks.seurat.object)[Seurat::Idents(peaks.seurat.object)!=population.1])
  } else {
    if (length(population.2) == 1) { # cluster identity used as input
      remainder.set = names(Seurat::Idents(peaks.seurat.object)[Seurat::Idents(peaks.seurat.object)==population.2])
    } else { # cell identity used as input
      remainder.set = population.2
    }
  }

  peak.names = rownames(peaks.seurat.object)

  # Get the peaks/APA sites expressed in the foreground set based on proportion of non-zeros
  this.data <- Seurat::GetAssayData(peaks.seurat.object, slot = "data", assay="RNA")
  nz.row.foreground = tabulate(this.data[, foreground.set]@i + 1, nbins = nrow(peaks.seurat.object))
  nz.prop.foreground = nz.row.foreground/length(foreground.set)
  peaks.foreground = peak.names[which(nz.prop.foreground > threshold)]

  # Now identify the peaks/APA sites expressed in the background set
  nz.row.background = tabulate(this.data[, remainder.set]@i + 1, nbins = nrow(peaks.seurat.object))
  nz.prop.background = nz.row.background/length(remainder.set)
  peaks.background = peak.names[which(nz.prop.background > threshold)]

  return(union(peaks.foreground, peaks.background))
}


############################################################
#'
#' Identify highly expressed peaks
#'
#' Selects peaks that are considered expressed above some provided criteria within a target or
#' background cluster. Considers peaks expressed in some x\% of cells to be highly expressed. Returns the
#' union of peaks identified from the target and background cluster
#'
#' @param peaks.sce.object the peak-count SCE object
#' @param population.1 target population
#' @param population.2 background population If NULL (deafult) all non-population.1 cells
#' @param threshold percentage threshold of detected (non-zero) expression for including a peak
#' @return an array of peak (or gene) names
#' @examples
#' \dontrun{
#' get_expressed_peaks_sce(peak.sce, "1")
#' get_expressed_peaks_sce(peak.sce, population.1 = "1", population.2 = "2")
#' }
get_expressed_peaks_sce <- function(peaks.sce.object, population.1, population.2=NULL, threshold=0.05) {

  if (length(population.1) == 1){ # cluster identity used as input
    foreground.set = names(colData(peaks.sce.object)$CellIdent[colData(peaks.sce.object)$CellIdent==population.1])
  } else { # cell identity used as input
    foreground.set = population.1
  }
  if (is.null(population.2)) {
    remainder.set = names(colData(peaks.sce.object)$CellIdent[colData(peaks.sce.object)$CellIdent!=population.1])
  } else {
    if (length(population.2) == 1) { # cluster identity used as input
      remainder.set = names(colData(peaks.sce.object)$CellIdent[colData(peaks.sce.object)$CellIdent==population.2])
    } else { # cell identity used as input
      remainder.set = population.2
    }
  }

  peak.names = rownames(peaks.sce.object)

  # Get the peaks/APA sites expressed in the foreground set based on proportion of non-zeros
  this.data <- SingleCellExperiment::counts(peaks.sce.object)
  nz.row.foreground = tabulate(this.data[, foreground.set]@i + 1, nbins = nrow(peaks.sce.object))
  nz.prop.foreground = nz.row.foreground/length(foreground.set)
  peaks.foreground = peak.names[which(nz.prop.foreground > threshold)]

  # Now identify the peaks/APA sites expressed in the background set
  nz.row.background = tabulate(this.data[, remainder.set]@i + 1, nbins = nrow(peaks.sce.object))
  nz.prop.background = nz.row.background/length(remainder.set)
  peaks.background = peak.names[which(nz.prop.background > threshold)]

  return(union(peaks.foreground, peaks.background))
}

############################################################

get_percent_expression <- function(peaks.object, this.cluster, remainder=FALSE, geneSet = rownames(seurat.object)) {

  if (class(peaks.object) == "Seurat") {
    if (length(this.cluster) == 1){ # cluster identity used as input
      foreground.set = names(Seurat::Idents(peaks.object)[Seurat::Idents(peaks.object)==this.cluster])
    } else { # cell identity used as input
      foreground.set = this.cluster
    }

    if (remainder) {
      cell.set <- setdiff(colnames(peaks.object), foreground.set)
    } else{
      cell.set <- foreground.set
    }

    peak.names = rownames(peaks.object)

    # Get the peaks/APA sites expressed in the foreground set based on proportion of non-zeros
    this.data <- Seurat::GetAssayData(peaks.object, slot = "data", assay="RNA")
    nz.row.cells = tabulate(this.data[, cell.set]@i + 1, nbins = length(peak.names))
    nz.prop.cells = nz.row.cells/length(cell.set)
    names(nz.prop.cells) = peak.names
    nz.prop.cells = nz.prop.cells[geneSet]

    return(nz.prop.cells)
  } else if (class(peaks.object) == "SingleCellExperiment") {

    if (length(this.cluster) == 1){ # cluster identity used as input
      foreground.set = names(colData(peaks.object)$CellIdent[colData(peaks.object)$CellIdent==this.cluster])
    } else { # cell identity used as input
      foreground.set = this.cluster
    }

    if (remainder) {
      cell.set <- setdiff(colnames(peaks.object), foreground.set)
    } else{
      cell.set <- foreground.set
    }

    peak.names = rownames(peaks.object)

    # Get the peaks/APA sites expressed in the foreground set based on proportion of non-zeros
    this.data <- SingleCellExperiment::counts(peaks.object)
    nz.row.cells = tabulate(this.data[, cell.set]@i + 1, nbins = length(peak.names))
    nz.prop.cells = nz.row.cells/length(cell.set)
    names(nz.prop.cells) = peak.names
    nz.prop.cells = nz.prop.cells[geneSet]

    return(nz.prop.cells)
  }

}
