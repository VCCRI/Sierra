#############################################################################
#' splitBam
#'
#' Utility to split a bam file into multiple bam files based on the barcode
#'
#' Given a bam file that was processed by CellRanger, splitBam splits the
#' bam into multiple bam files, one per cell barcode.
#' Bam file needs to have the barcode stored in the "CB" field.
#'
#' @param bam CellRanger outputted bam file with the CB field
#' @param cellbc.df data frame of the cell barcode, needs to have the column names: "celltype" and "cellbc"
#' @param outdir directory to output the bam files. The bam files will be called [celltype].bam. If NULL no BAM file created.
#' @param yieldSize number of lines of bam files to load. Default: 1000000
#' @param gtf_gr gene model genomic ranges. Only used if geneSymbol is defined.
#' @param geneSymbol Gene symbol. Used to identify the genomic coordinates to extract reads from.
#' @param gi_ext The number of nucleotides to extend the genomic interval in extracting reads from (default 50).
#' @param rle_output: If TRUE will generate and return rle_list object
#'
#' @return a rleList of coverage for each cell type
#'
#' @examples
#'
#' @export
#' extdata_path <- system.file("extdata",package = "scpolya")
#' load(paste(extdata_path,"TIP_vignette_gene_Seurat.RData",sep="/"))
#' cellbc.df <- data.frame(celltype=genes.seurat@active.ident, cellbc= names(genes.seurat@active.ident))
#' bam <- "c:/BAM/Harvey/scpolyA/one_percent.bam"
#' splitBam(bam, cellbc.df, "c:/TEMP/")
#'
#'
#' # Example 2 extract reads that overlap a gene
#'
#' gtf_file <- "u:/Reference/mm10/cellranger_genes.gtf.gz"
#' gtf_gr <- rtracklayer::import(gtf_file)
#' extdata_path <- system.file("extdata",package = "scpolya")
#' load(paste(extdata_path,"TIP_vignette_gene_Seurat.RData",sep="/"))
#' cellbc.df <- data.frame(celltype=genes.seurat@active.ident, cellbc= names(genes.seurat@active.ident))
#' bam <- "R:/scpolyA_BAM_link/possorted_genome_bam.bam"
#' splitBam(bam, cellbc.df, outdir="c:/TEMP/", gtf_gr=gtf_gr, geneSymbol="Dnajc19")
#'
#'
#' @export
SplitBam <- function(bam, cellbc.df, outdir=NULL, yieldSize = 1000000,
                     gtf_gr = NULL, geneSymbol=NULL, gi_ext = 50,
                     rle_output=FALSE) {
 # require(GenomicAlignments)
#  require(rtracklayer)

  message("splitting bam file: ", bam)


  if (! is.null(geneSymbol))
  {
    # Need check that gene_name field exists
    idx <-which(gtf_gr$gene_name == geneSymbol)
    if (length(idx) == 0)
    { warning("Could not find gene name. Please check spelling (and case)")
      return(NULL)
    }
    GenomeInfoDb::seqlevelsStyle(gtf_gr) <- "NCBI"

    # Work out the genomic range to extract from
    start <- min(start(ranges(gtf_gr[idx])))
    end <- max(end(ranges(gtf_gr[idx])))
    chrom <- as.character(GenomicRanges::seqnames(gtf_gr[idx]))[1]  # should I check that all returned chromosomes are the same?
    gene_strand <- as.character(strand(gtf_gr[idx]))[1]
    toExtract_gr <- GenomicRanges::GRanges(seqnames=chrom, ranges=IRanges::IRanges(start-gi_ext , width=end-start+gi_ext), strand=gene_strand)
    param <- Rsamtools::ScanBamParam(tag=c("CB", "UB"),which = toExtract_gr)
  }
  else
  {
    param <- Rsamtools::ScanBamParam(tag=c("CB", "UB"))
    geneSymbol <- "all"   # This will be incorporated into filename
    gene.provided <- NULL
  }

   cov_rle <- IRanges::RleList(compress=FALSE)     # Coverage list (i.e wig like). Populated for each cell type

   ctypes <- unique(cellbc.df$celltype)
   print(ctypes)
   for(eachtype in ctypes) {
      message("processing cell type ", eachtype)
      aln.per.type <- NULL
      cellbc <- subset(cellbc.df, celltype == eachtype)$cellbc
      bamfile <- Rsamtools::BamFile(bam, yieldSize=yieldSize)
      open(bamfile)
      while (length(chunk0 <- GenomicAlignments::readGAlignments(bamfile,param=param))) {
        if (! is.null(gene.provided))
        { # Only want to reads that are same strand as gene
          idx <- which(as.character(strand(chunk0)) == gene_strand)
          chunk0 <- chunk0[idx]
        }
        cat("chunk0:", length(chunk0), "length of aln: ", length(aln.per.type), "\n")

        if (length(aln.per.type) == 0)
          aln.per.type <- chunk0[which(S4Vectors::mcols(chunk0)$CB %in% cellbc)]
        else
          aln.per.type <- c(aln.per.type, chunk0[which(S4Vectors::mcols(chunk0)$CB %in% cellbc)])
      } # read in yieldSize number of records per iteration
      close(bamfile)

      outfile <- ''
      if (length(aln.per.type) == 0) {
        message("No data found for ", eachtype)
      }
      else if (is.null(outdir)){
        message(eachtype)
      }
      else {
        outfile <- paste0(outdir, eachtype, ".", geneSymbol,".bam")
        message("Writing to ", outfile)

       # as(aln.per.type, "GAlignments")
        rtracklayer::export(aln.per.type, Rsamtools::BamFile(outfile))
      }
      if(rle_output)
      {
        cov_rle  <- c(cov_rle, GenomicRanges::coverage(aln.per.type)[chrom])
        names(cov_rle)[length(cov_rle)] <- paste0(eachtype, ".", geneSymbol,".bam")
      }

   } ## Loop over cell types
  invisible(cov_rle)
}




#######################################################################
#' merge_bam_coverage
#'
#' @param bamfile : A list of BAM files that are to be merged
#'
#'
merge_bam_coverage <- function(bamfiles, to_extract)
{

  warning("Function is not finished .. expect errors?")
  for(i in bamfiles)
  {
    bf <-Rsamtools::BamFile(i)

    open(bf)
    chunk0 <- GenomicAlignments::readGAlignments(bf)
    # gr <-GenomicRanges::GRanges(chunk0)
    GenomeInfoDb::seqlevelsStyle(chunk0) <- "UCSC"
    close(bf)
    idx <- which(as.character(BiocGenerics::strand(chunk0)) == gene_strand)
    tmp <-GenomicRanges::coverage

    gr <- GenomicRanges::GRanges(seqnames=chrom, ranges=IRanges::IRanges(start:end, width=1), strand=gene_strand)
    S4Vectors::mcols(gr) <- as.numeric(tmp[[chrom]])[start:end]

  }
  return (bam_coverage)
}

##################################################################
#' geneToGR converts a gene sybol to genomic ranges coordinate
#'
#' @param geneID : Gene symbol
#' @param gtf_gr : Granges object of a gtf file
#'
#' gtf_file <- "u:/Reference/mm10/cellranger_genes.gtf.gz"
#' gtf_gr <- rtracklayer::import(gtf_file)
#' gtf_gr=gtf_gr, geneSymbol="Dnajc19"
#' geneGR  <- geneToGR(gtf_gr=gtf_gr, geneSymbol="Dnajc19")
#'
geneToGR <- function(geneID, gtf_gr)
{
  if (! is.null(geneSymbol))
  {
    # Need check that gene_name field exists
    idx <-which(gtf_gr$gene_name == geneSymbol)
    if (length(idx) == 0)
    { warning("Could not find gene name. Please check spelling (and case)")
      return(NULL)
    }
    GenomeInfoDb::seqlevelsStyle(gtf_gr) <- "NCBI"

    # Work out the genomic range to extract from
    start <- min(GenomicRanges::start(GenomicRanges::ranges(gtf_gr[idx])))
    end <- max(GenomicRanges::end(GenomicRanges::ranges(gtf_gr[idx])))
    chrom <- as.character(GenomicRanges::seqnames(gtf_gr[idx]))[1]  # should I check that all returned chromosomes are the same?
    gene_strand <- as.character(GenomicRanges::strand(gtf_gr[idx]))[1]
    gr <- GenomicRanges::GRanges(seqnames=chrom, ranges=IRanges::IRanges(start, width=end - start), strand=gene_strand)
  }
  return(gr)
}


##############################################################3
#'
#'
#' Example of how to use function:
#'
#' library("data.table")
#' endothelial_cov <- fread(file="c:/BAM/Harvey/scpolyA/Porrello_Support_Files/Porrello_Endothelial.F-CycCl_vs_F-Act.wig.txt.gz", sep = "\t", header = TRUE)
#' EC_coverage.rle <- seqmonk_to_rle(endothelial_cov, col_idx = 13:28)
#'
seqmonk_to_rle <- function(df, col_idx = 13:28)
{
  colnames(df)[2:5] <- c("chrom", "start","end", "strand")
  sampleIDs <- colnames(df)[col_idx]
  coverage.rle <- list()

  for(i in col_idx)
  {
    coverage.rle[[length(coverage.rle)+1]] <- rle(df[,..i] )
    print(length(coverage.rle))
  }
  names(coverage.rle) <- sampleIDs
  coverage.rle$gr <- GenomicRanges::makeGRangesFromDataFrame(df[,2:5])
}

##################################################################
#'
#'
seqmonk_file_to_rle <- function(fn)
{
  coverage.rle <- list()

  # Quick scan of file to identify column names
  temp <- data.table::fread(file=fn, sep = "\t",nrows = 2,header = TRUE)
  all_col_names <- colnames(temp)
  col_idx <- 13:length(all_col_names)

  # Read in genomic coordinates and generate a genomic range object
  col_to_keep <- c("Chromosome","Start","End","Probe Strand")
  df <- data.table::fread(file=fn, sep = "\t", header = TRUE, select = col_to_keep)
#  browser()

  colnames(df) <- c("chrom", "start","end", "strand")
  coverage.rle$gr <- GenomicRanges::makeGRangesFromDataFrame(df)


  # Now extract all columns as rle objects

  for(i in all_col_names[col_idx])
  {
    print(paste("Starting :",i))
    df <- data.table::fread(file=fn, sep = "\t", header = TRUE, select = i)
    if (! is.null(df))
    {  coverage.rle[[length(coverage.rle)+1]] <- rle(as.numeric(unlist(df)))
      print(paste(i , ": complete"))
    }
    else
    { print(paste("NO DATA for:",i))}
  }
  names(coverage.rle) <- c("gr",all_col_names[col_idx])
  return(coverage.rle)
}

################################################################3
#'
#' load(file="c:/BAM/scRNA_polyA/FC.RData")
#' gtf_file <- "u:/Reference/mm10/cellranger_genes.gtf.gz"
#' gtf_gr <- rtracklayer::import(gtf_file)
#'
rle_to_WIG <- function(rle_input, gtf_gr=gtf_gr, geneSymbol="Dnajc19")
{
  toExtract <-geneToGR(geneID=geneSymbol, gtf_gr)
  tmp <- GenomicRanges::findOverlaps(FC$gr,toExtract)
  idx <- S4Vectors::queryHits(tmp)
  min_idx <- min(idx)
  max_idx <- max(idx)
  # Next need to extract idx coordinates from rle_input

}