#########################################################################################################
## AnnotatePeaksFromGTF
##
#' Annotates a set of peak coordinates from a GTF
#'
#' Annotate a set of peak coordinates according to genomic features the coordinates fall on -
#' 3'UTR, exon, intron and 5'UTR, and annotate proximity to motifs. Motifs include the
#' canonical polyA motif, A-rich regions and T-rich regions.
#'
#' @param peak.sites.file a file of peak coordinates.
#' @param gtf.file GTF reference file.
#' @param output.file file to write the annotations to.
#' @param genome genome object. If NOT NULL then will perform pA motif analysis.
#' @param invert_strand Boolean to signifiy if strand of gr peaks should be inversed
#' @param annotationType can be assigned "any" or "within". Default is "any" which states that the peak with gr must overlap annotation feature (eg exon)
#' @param transcriptDetails Boolean. If false will only return gene name. If true will return internal transcript position feature (eg exon/intron)
#' @param annotation_correction Boolean. When multiple overlapping genes are identified will prioritise gene based on annotation. 3'UTR annotation trumps all other annotation.
#' @param pA_motif_max_position Any AAUAAA after this position are not considered (default 50nt)
#' @param AAA_motif_min_position Any polyA/polyT stretches before this postion are not considered (default 10)
#' @param polystretch_length : the length of A or T to search for (default 13)
#' @param max_mismatch number of allowed mismatches for motif matching (default 1)
#' @param append.chr.peaks : When TRUE (default) appends the character "chr" on chromosome entry in peaks file. 
#' @param check.chr if TRUE (default) and append.chr.peaks is also TRUE, check whether "chr" characters have already been added. 
#'
#' @return NULL. writes output to file
#'
#' @examples 
#' 
#' extdata_path <- system.file("extdata",package = "Sierra")
#' peak.merge.output.file <- paste0(extdata_path, "/TIP_merged_peaks.txt")
#' reference.file <- paste0(extdata_path,"/Vignette_cellranger_genes_subset.gtf")
#' 
#' 
#'  genome <- BSgenome.Mmusculus.UCSC.mm10::BSgenome.Mmusculus.UCSC.mm10
#' 
#' 
#'  AnnotatePeaksFromGTF(peak.sites.file = peak.merge.output.file, 
#'                     gtf.file = reference.file, 
#'                     output.file = "TIP_merged_peak_annotations.txt", 
#'                     genome = genome)
#'
#'
#' @export
#'
AnnotatePeaksFromGTF <- function(peak.sites.file,
                                 gtf.file,
                                 output.file,
                                 genome = NULL,
                                 invert_strand = FALSE,
                                 annotationType ="any",
                                 transcriptDetails = TRUE,
                                 annotation_correction = TRUE,
                                 pA_motif_max_position = 50,
                                 AAA_motif_min_position = 10,
                                 polystretch_length=13,
                                 max_mismatch=1,
                                 append.chr.peaks = TRUE,
                                 check.chr = TRUE) {

  ## Import the GTF reference
  gtf_gr <- rtracklayer::import(gtf.file)
  gtf_TxDb <- GenomicFeatures::makeTxDbFromGFF(gtf.file, format="gtf")

  ## Read in the peaks
  peak.table <- read.table(peak.sites.file, header = TRUE, sep="\t", stringsAsFactors = FALSE)
  print(paste("Annotating ", nrow(peak.table), " peak coordinates."))

  ## First need to convert strand labels for compatibility with GenomicRanges
  all.peaks <- peak.table$polyA_ID
  strand = sub(".*:.*:.*-.*:(.*)", "\\1", all.peaks)
  strand = plyr::mapvalues(x = strand, from = c("1", "-1"), to = c("+", "-"))
  peak.remainder = sub(".*:(.*:.*-.*):.*", "\\1", all.peaks)
  peaks.use = paste0(peak.remainder, ":", strand)

  ## Also need to ensure MT chromosomes are labelled 'M'
  chrs.all <- sub("(.*):.*-.*:.*", "\\1", peaks.use)
  if (length(which(chrs.all == "MT")) > 0)
  { chrs.all <- plyr::mapvalues(chrs.all, from="MT", to="M") }
  
  peaks.use.chr.update <- peaks.use
  # Append "chr" to peaks to match 
  if (append.chr.peaks == TRUE) { 
    
    if (check.chr == TRUE) {
      ## Check if 'chr' characters are already preceding chromosome names
      chr.counts <- sum(startsWith(unique(chrs.all), "chr") == TRUE)
      if (chr.counts > 0) {
        print(paste("\'chr\' character(s) already preceding chromosome name for",
                      chr.counts, "chromosomes. Will skip adding \'chr\' characters",
                      "for this chromosome set. Set check.chr=FALSE if the addition",
                      "of \'chr\' characters is required."))
      } else {
        ## If not add the chr prefix
        chrs.all <- paste0("chr", chrs.all)
      }
    } else {
      chrs.all <- paste0("chr", chrs.all)
    }
    peaks.use.chr.update <- paste0(chrs.all, sub(".*(:.*-.*:.*)", "\\1", peaks.use))
    
  }
  
  ## If genome object provided, check what chromosome names are matching
  if (!is.null(genome) & isS4(genome)) {  
    
    chr.set <- unique(chrs.all)
    available.chr <- BSgenome::seqnames(genome)
    peaks.diff <- setdiff(chr.set, available.chr)
    peaks.overlap <- intersect(chr.set, available.chr)
    if (length(peaks.diff) > 0) {
      print(paste("The following chromosome names are present in the peak file, but not in the genome file, and will not be annotated:", 
                    paste(peaks.diff, collapse = ", ")))
    }
    peaks.keep.idx <- which(chrs.all %in% peaks.overlap)
    peaks.use.chr.update <- peaks.use.chr.update[peaks.keep.idx]
    all.peaks <- all.peaks[peaks.keep.idx]
  }

  gr <- GenomicRanges::GRanges(peaks.use.chr.update)

  ## Everything should be in order - run annotation
  annot.df <- annotate_gr_from_gtf(gr = gr,
                                   gtf_gr = gtf_gr,
                                   gtf_TxDb = gtf_TxDb,
                                   genome = genome,
                                   invert_strand = invert_strand,
                                   annotationType = annotationType,
                                   transcriptDetails = transcriptDetails,
                                   annotation_correction = annotation_correction,
                                   pA_motif_max_position = pA_motif_max_position,
                                   AAA_motif_min_position = AAA_motif_min_position,
                                   polystretch_length = polystretch_length,
                                   max_mismatch = max_mismatch
  )
  rownames(annot.df) <- as.character(all.peaks)
  
  ## As a final step add the junctions to the output
  annot.df$Junctions <- peak.table[peaks.keep.idx, "exon.intron"]

  write.table(annot.df, file = output.file, quote = FALSE, sep = "\t")
}


################################################################################################
##
## This function has been designed to be called from annotate_gr_from_gtf
#'
#' This function has been designed to be called from annotate_gr_from_gtf
#'
#' @param gr a granges object of peaks to annotate
#' @param reference_gr a granges object of annotation info
#' @param annotationType a granges object of peaks to annotate
#'
gene_Labels<- function(gr, reference_gr, annotationType)
{
  all_hits <- GenomicAlignments::findOverlaps(gr , reference_gr, type= annotationType)
  if (length(reference_gr) ==0)
  { warning("No entries in reference to annotate. Cannot continue")
    return(NULL)
  }
  if (length(all_hits) == 0)
  { sanityCheck <- length(intersect(GenomeInfoDb::seqlevels(gr), GenomeInfoDb::seqlevels(reference_gr)))
    msg <- paste0("No peaks aligned to any entry within gtf reference.",
                      "\n Sanity check: ", sanityCheck,
                      " seqnames (i.e. chromosomes) match between peak and reference file")
    warning(msg)
    return(NULL)
  }

  identified_gene_symbols <- reference_gr[S4Vectors::subjectHits(all_hits)]$gene_name
  idx_to_annotate <- S4Vectors::queryHits(all_hits)

  multi_annotations <- which(table(idx_to_annotate) > 1)
  unique_annotations <- unique(idx_to_annotate)

  multi_gene_IDs <- unlist(sapply(names(multi_annotations),FUN=function(x) {
    newID <- paste(unique(identified_gene_symbols[which(x== idx_to_annotate)]),collapse=",")
    # rep(newID, length(which(x== idx_to_annotate)))
  }))
  multi_idx <- as.numeric(names(multi_gene_IDs))  # These are the indexes to annotate

  to_convert <- lapply(multi_idx,FUN = function(x) {which(idx_to_annotate == x)})

  if (length(to_convert) > 0)
  {   for(i in 1:length(to_convert))
  {
    identified_gene_symbols[to_convert[[i]]] <- multi_gene_IDs[[i]]
  }
  }

  return(list (identified_gene_symbols=identified_gene_symbols, idx_to_annotate=idx_to_annotate))
}

#########################################################################################################
## annotate_gr_from_gtf
##
#' Annotates a granges object with overlapping genes from gtf file.
#'
#'  gr is the genomic ranges that need to be annotation. Ideally original input should be in the format:
#'
#'       chr8:70331172-70331574:+   # chr:start-end:strand
#'
#'   This could already exist within an R object or you can copy it in via readClipboard.
#'
#'  gr <- GRanges(readClipboard())
#'
#' You need to run the following code:
#'         gtf_file <- "u:/Reference/hg38/hg38_gene.gtf.gz"
#'         gtf_file <- "u:/Reference/mm10/mm10_gene.gtf.gz"
#'         gtf_file <- "u:/Reference/mm10/cellranger_genes.gtf.gz"
#'        gtf_gr <- rtracklayer::import(gtf_file)
#'        gtf_TxDb <- GenomicFeatures::makeTxDbFromGFF(gtf_file, format="gtf")
#'
#' annotationType can be c("any", "start", "end", "within", "equal"),
#' @param gr a granges object of peaks to annotate
#' @param invert_strand Boolean to signifiy if strand of gr peaks should be inversed
#' @param gtf_gr granges gtf file that contains annotation information
#' @param annotationType can be assigned "any" or "within". Default is "any" which states that the peak with gr must overlap annotation feature (eg exon)
#' @param transcriptDetails Boolean. If false will only return gene name. If true will return internal transcript position feature (eg exon/intron)
#' @param gtf_TxDb  same as gtf_gr but as a TxDb object.
#' @param annotation_correction Boolean. When multiple overlapping genes are identified will
#' prioritise gene based on annotation. 3'UTR annotation trumps all other annotation.
#' @param genome genome object. If NOT NULL then will perform pA motif analysis.
#' @param pA_motif_max_position Any AAUAAA after this position are not considered (default 50nt)
#' @param AAA_motif_min_position Any polyA/polyT stretches before this postion are not considered (default 10)
#' @param polystretch_length : the length of A or T to search for (default 13)
#' @param max_mismatch : The number of mismatches tolerated in polystretch
#'
#' @examples 
#' library(Sierra)
#' 
#'  # Generate peaks for Cxcl12, Arhgap10, Mast4,  using mm10 coordinates:
#'  gr_peaks <- GenomicRanges::GRanges(c("chr6:117174600-117175065:+",
#'             "chr6:117180975-117181367:+",
#'             "chr8:77250366-77250686:-",
#'             "chr8:77426400-77517833:-",
#'             "chr13:102905701-102906230:-",
#'             "chr13:103139934-103171545:-"))
#'             
#'  # Load other files from vignette           
#'  extdata_path <- system.file("extdata",package = "Sierra")
#'  reference.file <- paste0(extdata_path,"/Vignette_cellranger_genes_subset.gtf")
#'  
#'  # convert gtf file to both granges and a TXDb object
#'  gtf_gr <- rtracklayer::import(reference.file)
#'  gtf_TxDb <- GenomicFeatures::makeTxDbFromGFF(reference.file, format="gtf")
#'  
#'  genome <- BSgenome.Mmusculus.UCSC.mm10::BSgenome.Mmusculus.UCSC.mm10          
#'   
#'  annotate_gr_from_gtf(gr = gr_peaks, gtf_gr = gtf_gr,
#'                       gtf_TxDb = gtf_TxDb, genome = genome)
#'
#'  annotate_gr_from_gtf(gr = gr_peaks, gtf_gr = gtf_gr,
#'                       gtf_TxDb = gtf_TxDb, genome = genome, transcriptDetails=TRUE)                                            
#' @return a dataframe with appended columns containing annotation
#' @export
##
## Written March 2019
annotate_gr_from_gtf <- function(gr, invert_strand = FALSE, gtf_gr = NULL,
                       annotationType ="any",
                       transcriptDetails = FALSE, gtf_TxDb,
                       annotation_correction = TRUE, genome = NULL,
                       pA_motif_max_position = 50,
                       AAA_motif_min_position = 10,
                       polystretch_length=13, max_mismatch=1
                       )
{
  # Checking values passed to function are meaning full.
  if (is.null(gtf_gr))
  {warning("No gtf file provided")
    return(NULL)
  }
  #
  if (invert_strand)
  { gr <- invertStrand(gr) }

  GenomicRanges::mcols(gtf_gr) <- GenomicRanges::mcols(gtf_gr)[c("type","gene_id","gene_name")]



  # check for compatibility between chromosome labels:
  if (length(intersect(GenomeInfoDb::seqlevels(gr), GenomeInfoDb::seqlevels(gtf_gr))) == 0)
  { # remove chr prefix from both data sets
    GenomeInfoDb::seqlevels(gr) <-   gsub(pattern = "chr",replacement = "",x = GenomeInfoDb::seqlevels(gr))
    GenomeInfoDb::seqlevels(gtf_gr) <- gsub(pattern = "chr",replacement = "",x = GenomeInfoDb::seqlevels(gtf_gr))
  }
  
  # Do some simple checks to ensure granges object has required fields
  if (length(grep(pattern="type",x=colnames(as.data.frame(gtf_gr)))) ==0)
  {
    message("Annotation reference does not have a type field. This field is used
            to extract gene, transcript, exon, UTR information. Cannot continue")
    return(NULL)
  }

  # Check that the type field has required annotations
  listed_annotations <- names(table(gtf_gr$type))
  if (length(grep(pattern = "gene",x = listed_annotations)))
  { genes_gr <- gtf_gr[gtf_gr$type == "gene"] }
  else if (length(grep(pattern = "transcript",x = listed_annotations)))
  { genes_gr <- gtf_gr[gtf_gr$type == "transcript"] }
  else if (length(grep(pattern = "exon",x = listed_annotations)))
  { genes_gr <- gtf_gr[gtf_gr$type == "exon"] }
  else # No way can continue
  {
    message("Reference has no recognised labels to annotate. Cannot continue")
    return(NULL)
  }

  annotate_info <- gene_Labels(gr, genes_gr,annotationType)
  if (is.null(annotate_info))
  { warning("")
    return(NULL)
  }

  df <- as.data.frame(gr)
  df$gene_id <- ""
#  df$gene_id[idx_to_annotate] <-  identified_gene_symbols
#  df$gene_id[multi_idx] <- multi_gene_IDs

  df$gene_id[annotate_info$idx_to_annotate] <- annotate_info$identified_gene_symbols # identified_gene_symbols
  df_with_gene_labels <- df    # Making a copy. Annotation tracks will be labelled with gene IDs rather than "YES"

  if (transcriptDetails)
  {
    cat("\nAnnotating 3' UTRs")
    df$UTR3 <- ""
    UTR_3_GR <- GenomicFeatures::threeUTRsByTranscript(gtf_TxDb, use.names=TRUE)
    all_UTR_3_hits <- GenomicAlignments::findOverlaps(gr , UTR_3_GR,type = annotationType)
    idx_to_annotate_3UTR <- S4Vectors::queryHits(all_UTR_3_hits)
    df$UTR3[idx_to_annotate_3UTR] <- "YES"
    ok_to_annotate <- 1: nrow(df) # To identify which candidates can be annotated. Initally all samples can be annotated.
    if (annotation_correction)
    {
      # Want to record gene name (symbol) for 3'UTRs. The genomicfeatures function removes this info.
      # Therefore subset UTRs and then overlap them to defined 3'UTRs.
      # First need to know what annotations exist, and if 3'UTRs are actually annotated
      listed_annotations <- names(table(gtf_gr$type))
      if (length(grep(pattern = "three_prime_utr",x = listed_annotations)))
      {
        UTR_3_GR <- gtf_gr[gtf_gr$type == "three_prime_utr"]
        UTR_annotate_info <- gene_Labels(gr, UTR_3_GR ,annotationType)
        df_with_gene_labels$UTR3 <- rep(NA, nrow(df_with_gene_labels))
        df_with_gene_labels$UTR3[UTR_annotate_info$idx_to_annotate] <- UTR_annotate_info$identified_gene_symbols
      }
      else if (length(grep(pattern = "UTR",x = listed_annotations)))
      { 
        UTR_GR <- gtf_gr[gtf_gr$type == "UTR"]   # This will be used to retrieve gene names from ALL UTRs
        real_3UTRs_idx <- GenomicAlignments::findOverlaps(UTR_GR , UTR_3_GR,type = annotationType)

        UTR_annotate_info <- gene_Labels(gr, UTR_GR[S4Vectors::queryHits(real_3UTRs_idx)] ,annotationType)
        df_with_gene_labels$UTR3 <- rep(NA, nrow(df_with_gene_labels))
        df_with_gene_labels$UTR3[UTR_annotate_info$idx_to_annotate] <- UTR_annotate_info$identified_gene_symbols
      }
      else
      {
        warning("Setting annotation_correction to FALSE as GTF does not have type metadata field containing 'UTR' or 'three_prime_utr'")
        annotation_correction <- FALSE
      }

      # Copy relevant updated annotations
      # Grab index of annotated entries and copy to main df.
      if (annotation_correction)
      { df$gene_id[UTR_annotate_info$idx_to_annotate] <- UTR_annotate_info$identified_gene_symbols
        ok_to_annotate <- setdiff(ok_to_annotate, UTR_annotate_info$idx_to_annotate)
      }
    }

    cat("\nAnnotating 5' UTRs")
    df$UTR5 <- ""
    UTR_5_GR <- GenomicFeatures::fiveUTRsByTranscript(gtf_TxDb)
    all_UTR_5_hits <- GenomicAlignments::findOverlaps(gr , UTR_5_GR,type = annotationType)
    identified_5UTRs <- UTR_5_GR[S4Vectors::subjectHits(all_UTR_5_hits)]$exon_id
    idx_to_annotate_5UTR <- S4Vectors::queryHits(all_UTR_5_hits)
    if (annotation_correction)
    { # We don't annotate entries that were previously annotated to derive from a 3'UTR of a single gene.
      idx_to_annotate_5UTR <- intersect(idx_to_annotate_5UTR, ok_to_annotate)
    }
    df$UTR5[idx_to_annotate_5UTR] <- "YES"


    cat("\nAnnotating introns")
    df$intron <- ""
    introns_GR <- GenomicFeatures::intronsByTranscript(gtf_TxDb)
    all_intron_hits <- GenomicAlignments::findOverlaps(gr , introns_GR, type = annotationType)
    identified_introns <- introns_GR[S4Vectors::subjectHits(all_intron_hits)]
    idx_to_annotate_introns <- S4Vectors::queryHits(all_intron_hits)
    if (annotation_correction)
    { # We don't annotate entries that were previously annotated to derive from a 3'UTR of a single gene.
      idx_to_annotate_introns <- intersect(idx_to_annotate_introns, ok_to_annotate)
    }
    df$intron[idx_to_annotate_introns] <- "YES"

    cat("\nAnnotating exons")
    my_exons <- gtf_gr[gtf_gr$type == "exon"]
    all_exon_hits <-  GenomicAlignments::findOverlaps(gr , my_exons, type = annotationType)
    idx_to_annotate_exons <- S4Vectors::queryHits(all_exon_hits)
    df$exon <- ""
    if (annotation_correction)
    { # We don't annotate entries that were previously annotated to derive from a 3'UTR of a single gene.
      idx_to_annotate_exons <- intersect(idx_to_annotate_exons, ok_to_annotate)
    }
    df$exon[idx_to_annotate_exons] <- "YES"

    cat("\nAnnotating CDS")
    my_CDS <- gtf_gr[gtf_gr$type == "CDS"]
    all_CDS_hits <-  GenomicAlignments::findOverlaps(gr , my_CDS, type = annotationType)
    idx_to_annotate_CDS <- S4Vectors::queryHits(all_CDS_hits)
    df$CDS <- ""
    if (annotation_correction)
    { # We don't annotate entries that were previously annotated to derive from a 3'UTR of a single gene.
      idx_to_annotate_CDS <- intersect(idx_to_annotate_CDS, ok_to_annotate)
    }
    df$CDS[idx_to_annotate_CDS] <- "YES"
  }

  if (! is.null(genome))
  { cat("\nAnalysing genomic motifs surrounding peaks (this can take some time)\n")
    if (isS4(genome))
    {
      ## Check if 'chr' character is/should be appended to chromosome number
      available.chr <- BSgenome::seqnames(genome)
      chr.counts <- sum(startsWith(available.chr, "chr") == TRUE)
      if (grepl("chr.", as.character(gr@seqnames)[1]) == FALSE & chr.counts > 0) 
      {  gr <- paste("chr",as.character(gr),sep='') }

      # Set up progress bar so user can "watch" progress
      pb <- txtProgressBar(min = 0, max = length(gr), style = 3)
      motif_details <- lapply(X = 1:length(gr),
                                FUN = function(i, gr, pb) {
                                  setTxtProgressBar(pb, i)
                                  nextgr <- gr[i]
                                  BaseComposition(genome,coord=nextgr,mismatch=max_mismatch, AT_length=polystretch_length)
                                  }, 
                                  gr,pb
                              )
      cat("\n")  # Keep output neat as text bar does not append return character
      

      df$pA_motif <- unlist( lapply(motif_details, FUN= function(x) {
            pA_motif_position <- FALSE
            motif_pos <- unlist(x$pA_motif_pos)
            if (length(motif_pos) > 0)
            {
              if (!is.na(motif_pos) & length(motif_pos) > 0)
                pA_motif_position <- (max(motif_pos) < pA_motif_max_position)
            }
            return (pA_motif_position) } ))
  
      df$pA_stretch <- unlist( lapply(motif_details, FUN= function(x) {
          pA_stretch_position <- FALSE
          motif_pos <- unlist(x$pA_stretch_pos)
          if (length(motif_pos) > 0)
          {
            if (!is.na(motif_pos) & length(motif_pos) > 0)
              pA_stretch_position <- (max(motif_pos) > AAA_motif_min_position)
          }
          return (pA_stretch_position) } ))
  
      df$pT_stretch <- unlist( lapply(motif_details, FUN= function(x) {
          pT_stretch_position <- FALSE
          motif_pos <- unlist(x$pT_stretch_pos)
          if (length(motif_pos) > 0)
          {
            if (!is.na(motif_pos) & length(motif_pos) > 0)
              pT_stretch_position <- (max(motif_pos) > AAA_motif_min_position)
          }
          return (pT_stretch_position) } ))
    }
    else
    {
      warning("Genome object is not a BSgenome S4 object.
              Cannot annotate for sequence motifs")
    } # if (isS4(genome))

    
     # BaseComposition(genome=genome, coord = as.character(gr))

  }
  return(df)
}


#############################################################################
##
##
#' Identify polyA motif and/or polyA stretches from provided genomic coordinates
#'
#'
#'       chrom <- chr8
#'
#'       start <- 70331172
#'
#'       stop <- 70331574
#'
#'       strand <- "+"
#'
#'
#' You need to run the following code:
#'      genome <-  GenomicFeatures::makeTxDbFromGFF(gtf_file, format="gtf")
#'
#' annotationType can be c("any", "start", "end", "within", "equal"),
#' @param genome genome object of organism.
#' @param chrom chromosome
#' @param start Upstream start position
#' @param stop downstream end position
#' @param strand '+' or '-'. This will define the applied direction of offset
#' @param coord coordinates
#' @param offset The
#' @param length How many nucleotides of DNA sequence to return
#' @param mismatch : The max number of mismatches allowed in poly A/T stretch (default 1)
#' @param AT_length : length of A/T to search for within input sequence (default 13)
#' @return a dataframe with appended columns containing annotation
#'
#' chrom <- 'chr16'
#' start <- 49896378
#' stop  <- 49911102
#' strand <- '+'
#' genome <- BSgenome.Mmusculus.UCSC.mm10::BSgenome.Mmusculus.UCSC.mm10
#'
#' # Tmem126a intronic peak that coincides with a poly(A) rich region.
#' chrom  <- 'chr7'
#' start  <- 90451180
#' stop   <- 90451380
#' strand <- '-'
#' output <-BaseComposition(genome=genome, chrom=chrom, start=start, stop=stop, strand=strand)
#'
#' # Dync1h1 intronic peak that coincides with a long poly(T) rich region
#' coord <- "chr12:110609400-110609800:1"
#' output <-BaseComposition(genome=genome, coord=coord)
#'
#' TO DO:
#' * If peak falls at end of exon then need to obtain sequence from next exon. This
#' would require passing exon junction information.
#'
#' @export
##
## Written August 2019
BaseComposition <- function(genome=NULL,  chrom=NULL, start=NULL, stop=NULL, strand=NULL,
                            coord = NULL, offset = -50, length = 250,
                            mismatch=1, AT_length=13)
{
  # Check inputs
  if (! isS4(genome))
  { warning("genome is not a BSgenome S4 object. Cannot continue.")
    return(NULL)
  }
  
  if (! is.null(coord))
  {
    if (! is.null(genome) & (! is.null(chrom) | ! is.null(start) |
          ! is.null(stop) | ! is.null(strand)))
      warning("Multiple coodinates passed, will be using what was passed to coord")
    # "Cd47:16:49896378-49911102:1"
    # "Cd47:16:49914609-49915010:1"
    # "Tmem126a:7:90451180-90451380:-1"


    if (length(gregexpr(":",coord)[[1]]) == 3) # coord should look like this "Cd47:16:49896378-49911102:1"
    {
      coord_detail <- strsplit(coord, split = ":")
      if (grepl("chr.", coord_detail[[1]][2]) == FALSE) {
        chrom <- paste("chr",coord_detail[[1]][2],sep='')
      }
      strand <- coord_detail[[1]][4]

      if (strand == 1)
      { strand <- '+'   }
      else if (strand == -1)
      { strand <- '-'   }

      coord_detail <- strsplit(coord_detail[[1]][3],"-")
      start <- as.numeric(coord_detail[[1]][1])
      stop <- as.numeric(coord_detail[[1]][2])
    }

    if (length(gregexpr(":",coord)[[1]]) == 2) # coord should look like this "chr16:49896378-49911102:+"
    {
      coord_detail <- strsplit(coord, split = ":")
      chrom <- coord_detail[[1]][1]
      strand <- coord_detail[[1]][3]

      coord_detail <- strsplit(coord_detail[[1]][2],"-")
      start <- as.numeric(coord_detail[[1]][1])
      stop <- as.numeric(coord_detail[[1]][2])
    }
  }
  
  if ((is.na(start)) | (is.na(stop)))
  { mesg <- paste("Could not define/identify start stop coordinates from coord:", coord)
    warning(mesg)
    return(NULL)
  }

  # Require start to be smaller number
  if ( start > stop)
  {
    tmp <- start
    start <- stop
    stop <- tmp
  }

  # Now to extract sequence up and downstream of peak (relative to +ve strand)
  # First extract downstream sequence
  ## Default parameters for '+' strand
  seq_start_position <- stop + offset
  seq_end_position <- stop + length + offset
  
  if(chrom == "chrMT")
  {
    chrom <- 'chrM'
  }

  sequ <- BSgenome::getSeq(genome, chrom, seq_start_position, seq_end_position)   # Always get +ve strand
  sequ_upstream <- {}

  if (strand == '-')  # sequ is actually upstream if we are thinking abour -ve strand
  {   sequ_upstream <- sequ   }

  # Now extract upstream sequence
  seq_start_position <- start - length -offset
  seq_end_position <- start - offset

  sequ_tmp <- BSgenome::getSeq(genome, chrom, seq_start_position, seq_end_position)   # Always get +ve strand

  if (strand == '-')
  {
    sequ <- Biostrings::reverseComplement(sequ_tmp)
    sequ_upstream <- Biostrings::reverseComplement(sequ_upstream)
  } else # strand '+'
  {
    sequ_upstream <- sequ_tmp
  }

  #<- longestConsecutive(seq, "A") # returns a integer
  A_pattern <- paste(rep("A",AT_length),collapse='')
  T_pattern <- paste(rep("T",AT_length),collapse='')
  pA_motif  <-  Biostrings::matchPattern(pattern = "AATAAA", subject = sequ)
  pA_stretch <- Biostrings::matchPattern(pattern=A_pattern, subject=sequ, max.mismatch=1)
  pT_stretch <- Biostrings::matchPattern(pattern=T_pattern, subject=sequ_upstream, max.mismatch=1)

  return(list(pA_motif_pos = start(pA_motif)[1], pA_stretch_pos = start(pA_stretch)[1],
              pT_stretch_pos = start(pT_stretch)[1],
              sequence = list(upstream=sequ_upstream, downstream=sequ) ))
}


###########################################################################
##
## convert_coord will attempt to identfy coordinates and convert to a universal
##    format
##
##
##  all_coords_formatted <- unlist(lapply(X = all_coords,convert_coord))
##
convert_coord <- function(coord)
{
  # Currently will convert the following:
  #
  #  "Prkar1a:11:109669227-109669656:1"     to chr11:109669227-109669656:+

  if (length(gregexpr(":",coord)[[1]]) == 3) # coord should look like this "Cd47:16:49896378-49911102:1"
  {
    coord_detail <- strsplit(coord, split = ":")
    if (grepl("chr.", coord_detail[[1]][2]) == FALSE) {
      chrom <- paste("chr",coord_detail[[1]][2],sep='')
    }
    strand <- coord_detail[[1]][4]

    if (strand == 1)
    { strand <- '+'   }
    else if (strand == -1)
    { strand <- '-'   }
    coord_to_return <- paste(chrom,":",coord_detail[[1]][3],":",strand,sep='')

  }
  return(coord_to_return)
}



