
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
  if (length(all_hits) == 0)
  {
    warning("No samples matched")
    return(-1)
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
#' @param AAA_motif_min_position Any polyA stretches before this postion are not considered (default 10)
#'
#' @return a dataframe with appended columns containing annotation
##
## Written March 2019
annotate_gr_from_gtf <- function(gr, invert_strand = FALSE, gtf_gr = NULL,
                       annotationType ="any",
                       transcriptDetails = FALSE, gtf_TxDb,
                       annotation_correction = TRUE, genome = NULL,
                       pA_motif_max_position = 50,
                       AAA_motif_min_position = 10)
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

  genes_gr <- gtf_gr[gtf_gr$type == "gene"]

  annotate_info <- gene_Labels(gr, genes_gr,annotationType)

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
      else
      {
        UTR_GR <- gtf_gr[gtf_gr$type == "UTR"]   # This will be used to retrieve gene names from ALL UTRs
        real_3UTRs_idx <- GenomicAlignments::findOverlaps(UTR_GR , UTR_3_GR,type = annotationType)

        UTR_annotate_info <- gene_Labels(gr, UTR_GR[S4Vectors::queryHits(real_3UTRs_idx)] ,annotationType)
        df_with_gene_labels$UTR3 <- rep(NA, nrow(df_with_gene_labels))
        df_with_gene_labels$UTR3[UTR_annotate_info$idx_to_annotate] <- UTR_annotate_info$identified_gene_symbols
      }

      # Copy relevant updated annotations
      # Grab index of annotated entries and copy to main df.
      df$gene_id[UTR_annotate_info$idx_to_annotate] <- UTR_annotate_info$identified_gene_symbols
      ok_to_annotate <- setdiff(ok_to_annotate, UTR_annotate_info$idx_to_annotate)
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

 # browser()


  if (! is.null(genome))
  {

    motif_details <- lapply(X = paste("chr",as.character(gr),sep=''),
                            FUN = function(x) {baseComposition(genome,coord=x)})

    df$pA_motif <- unlist( lapply(motif_details, FUN= function(x) {
          pA_motif_position <- (x[1] < pA_motif_max_position)
          if (is.na(pA_motif_position))
          {  pA_motif_position <- FALSE }
          return (pA_motif_position) } ))

    df$pA_stretch <- unlist( lapply(motif_details, FUN= function(x) {
          pA_stretch_position <- (x[2] > AAA_motif_min_position)
          if (is.na(pA_stretch_position))
          {  pA_stretch_position <- FALSE }
          return (pA_stretch_position) } ))

     # baseComposition(genome=genome, coord = as.character(gr))

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
#' @param end downstream end position
#' @param strand '+' or '-'. This will define the applied direction of offset
#' @param coord coordinates
#' @param offset The
#' @param length How many nucleotides of DNA sequence to return
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
#'
#' TO DO:
#' * If peak falls at end of exon then need to obtain sequence from next exon. This
#' would require passing exon junction information.
#'
#' @export
##
## Written August 2019
baseComposition <- function(genome=NULL,  chrom=NULL, start=NULL, stop=NULL, strand=NULL,
                            coord = NULL, offset = -50, length = 250)
{
  # Check inputs
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
      chrom <- paste("chr",coord_detail[[1]][2],sep='')
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


  # Require start to be smaller number
  if ( start > stop)
  {
    tmp <- start
    start <- stop
    stop <- tmp
  }

  # Default parameters for '+' strand
  seq_start_position <- stop
  seq_end_position <- stop + length



  # modify start / stop according to strand
  if (strand == '-')
  {     seq_start_position <- start - length
        seq_end_position <- start
        offset <- offset * -1
  }

  # apply offset
  seq_start_position <- seq_start_position + offset
  seq_end_position <- seq_end_position + offset


  sequ <- BSgenome::getSeq(genome, chrom, seq_start_position, seq_end_position)   # Always get +ve strand

  if (strand == '-')
    sequ <- reverseComplement(sequ)

#  browser()
  #<- longestConsecutive(seq, "A") # returns a integer
  pA_motif  <-  matchPattern(pattern = "AATAAA", subject = sequ)
  pA_stretch <- matchPattern(pattern="AAAAAAAAAAAAA", subject=sequ, max.mismatch=1)


  return(list(pA_motif_pos = start(pA_motif)[1], pA_stretch_pos = start(pA_stretch)[1], sequence = sequ ))
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
    chrom <- paste("chr",coord_detail[[1]][2],sep='')
    strand <- coord_detail[[1]][4]

    if (strand == 1)
    { strand <- '+'   }
    else if (strand == -1)
    { strand <- '-'   }
    coord_to_return <- paste(chrom,":",coord_detail[[1]][3],":",strand,sep='')

  }
  return(coord_to_return)
}



