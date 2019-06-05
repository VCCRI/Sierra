#########################################################################################################
## annotate_gr_from_gtf
##
#' Annotates a granges object with overlapping genes from gtf file.
#'
#'  gr is the genomic ranges that need to be annotation. Ideally original input should be in the format:
#'       chr8:70331172-70331574:+   # chr:start-end:strand
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
#' @param annotationType can be assigned "within" or "all". Default is "within" which states that the peak with gr must be within annotation feature (eg exon)
#' @param transcriptDetails Boolean. If false will only return gene name. If true will return internal transcript position feature (eg exon/intron)
#' @param gtf_TxDb  same as gtf_gr but as a TxDb object.
#' @param annotation_correction Boolean. When multiple overlapping genes are identified will 
#' prioritise gene based on annotation. 3'UTR annotation trumps all other annotation. 
#' @return a dataframe with appended columns containing annotation
##
## Written March 2019
annotate_gr_from_gtf <- function(gr, invert_strand = FALSE, gtf_gr = NULL,
                       annotationType ="within",
                       transcriptDetails = FALSE, gtf_TxDb,
                       annotation_correction = TRUE)
{
  if (is.null(gtf_gr))
  {warning("No gtf file provided")
    return(NULL)
  }

  if (invert_strand)
  { gr <- invertStrand(gr) }

#  mcols(gtf_gr) <- mcols(gtf_gr)[c("type","gene_id","gene_type","gene_name")]
  GenomicRanges::mcols(gtf_gr) <- GenomicRanges::mcols(gtf_gr)[c("type","gene_id","gene_name")]

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

  # check for compatibility between chromosome labels:
  if (length(intersect(seqlevels(gr),seqlevels(gtf_gr))) == 0)
  { # remove chr prefix from both data sets
    seqlevels(gr) <-   gsub(pattern = "chr",replacement = "",x = seqlevels(gr))
    seqlevels(gtf_gr) <- gsub(pattern = "chr",replacement = "",x = seqlevels(gtf_gr))
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
      {### working here
        
#  browser()      
        UTR_3_GR <- gtf_gr[gtf_gr$type == "three_prime_utr"]
        UTR_annotate_info <- gene_Labels(gr, UTR_3_GR ,annotationType)
        df_with_gene_labels$UTR[UTR_annotate_info$idx_to_annotate] <- UTR_annotate_info$identified_gene_symbols
        
      }
      else
      {
        
        UTR_GR <- gtf_gr[gtf_gr$type == "UTR"]   # This will be used to retrieve gene names from ALL UTRs
        real_3UTRs_idx <- GenomicAlignments::findOverlaps(UTR_GR , UTR_3_GR,type = annotationType)
        
        UTR_annotate_info <- gene_Labels(gr, UTR_GR[S4Vectors::queryHits(real_3UTRs_idx)] ,annotationType)
        df_with_gene_labels$UTR[UTR_annotate_info$idx_to_annotate] <- UTR_annotate_info$identified_gene_symbols
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
  return(df)
}
