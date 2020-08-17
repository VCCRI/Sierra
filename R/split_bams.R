
#############################################################################
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
#' @param rle_output If TRUE will generate and return rle_list object
#' @param exportFastqHeader If TRUE will generate a txt output file that has same prefix as bam file containing fastq header IDs
#' @param genomicRegion Granges object of genomic region to extract. Only used if geneSymbol not defined.
#' @param bamTags BAM field tag identifiers to extract. Default is c("CB", "UB").
#'
#' @return a rleList of coverage for each cell type
#'
#' @examples
#' library('Sierra')
#' \dontrun{
#' extdata_path <- system.file("extdata",package = "scpolya")
#' load(paste(extdata_path,"TIP_vignette_gene_Seurat.RData",sep="/"))
#' cellbc.df <- data.frame(celltype=genes.seurat@active.ident, 
#'                         cellbc= names(genes.seurat@active.ident))
#' bamfile <- c(paste0(extdata_path,"/Vignette_example_TIP_sham.bam")
#' 
#' SplitBam(bam, cellbc.df)
#' }
#'
#' # Example 2 extract reads that overlap a gene
#' 
#' extdata_path <- system.file("extdata",package = "Sierra")
#' gtf.file <- paste0(extdata_path,"/Vignette_cellranger_genes_subset.gtf")
#' gtf.gr <- rtracklayer::import(gtf.file)
#' 
#' load(paste(extdata_path,"TIP_vignette_gene_Seurat.RData",sep="/"))
#' cellbc.df <- data.frame(celltype=genes.seurat@active.ident, 
#'                        cellbc= names(genes.seurat@active.ident))
#'   
#' # Modify cellbc.df so that the barcodes match what is in the BAM file                     
#' cellbc.df$cellbc <- sub("(.*)-.*", "\\1", cellbc.df$cellbc)
#' cellbc.df$cellbc <- paste0(cellbc.df$cellbc, "-1")
#'                        
#'                        
#' bam.file <- paste0(extdata_path,"/Vignette_example_TIP_mi.bam")
#' outdir <-  tempdir()  # change this to a meaningful location
#' SplitBam(bam.file, cellbc.df, outdir=outdir, gtf_gr=gtf.gr, geneSymbol="Dnajc19")
#' 
#'
#' @export
#' 
SplitBam <- function(bam, cellbc.df, outdir=NULL, yieldSize = 1000000,
                     gtf_gr = NULL, geneSymbol=NULL, gi_ext = 50,
                     rle_output=FALSE, exportFastqHeader=FALSE, genomicRegion=NULL,
                     bamTags=c("CB", "UB")) {
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
    param <- Rsamtools::ScanBamParam(tag=bamTags, which = toExtract_gr, what=c('qname', 'flag', 'rname', 'strand', 'pos'))
    gene.provided <- geneSymbol
  }
  else if (! is.null(genomicRegion))
  {
    tryCatch({
      param <- Rsamtools::ScanBamParam(tag=bamTags, which = toExtract_gr, what=c('qname', 'flag', 'rname', 'strand', 'pos'))
    }, error = function(err) {
      stop(paste0("Problem detected with provided genomic region. Please ensure genomicRegion is a Granges object"))
      })
  }
  else
  {

    param <- Rsamtools::ScanBamParam(tag=bamTags)
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
        if (exportFastqHeader)
        {   outfile.readID <- paste0(outdir, eachtype, ".", geneSymbol,".txt")
            write((as.data.frame(aln.per.type)$qname), file=outfile.readID)
        }
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
#' merge_bam_coverage
#'
#' @param bamfiles : A list of BAM files that are to be merged
#'
#'
merge_bam_coverage <- function(bamfiles)
{

  warning("Function is not finished .. expect errors?")
  for(i in bamfiles)
  {
    bf <-Rsamtools::BamFile(i)

#    open(bf)
#    chunk0 <- GenomicAlignments::readGAlignments(bf)
#    GenomeInfoDb::seqlevelsStyle(chunk0) <- "UCSC"
#    close(bf)
#    idx <- which(as.character(BiocGenerics::strand(chunk0)) == gene_strand)
#    tmp <-GenomicRanges::coverage

#    gr <- GenomicRanges::GRanges(seqnames=chrom, ranges=IRanges::IRanges(start:end, width=1), strand=gene_strand)
#    S4Vectors::mcols(gr) <- as.numeric(tmp[[chrom]])[start:end]

  }
  return (bam_coverage)
}

##################################################################
#' geneToGR converts a gene symbol to genomic ranges coordinate
#' 
#' geneToGR converts a gene symbol to genomic ranges coordinate
#'
#' @param geneSymbol : Gene symbol
#' @param gtf_gr : Granges object of a gtf file
#'
#' @examples
#'     library('Sierra')
#'     extdata_path <- system.file("extdata",package = "Sierra")
#'     gtf.file <- paste0(extdata_path,"/Vignette_cellranger_genes_subset.gtf")
#'     gtf.gr <- rtracklayer::import(gtf.file)
#'     
#'     geneGR  <- geneToGR(geneSymbol= "Dnajc19",gtf_gr=gtf.gr)
#' @export
geneToGR <- function(geneSymbol, gtf_gr)
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
#
#
# Example of how to use function:
#
# library("data.table")
# endothelial_cov <- fread(file="c:/BAM/Harvey/scpolyA/Porrello_Support_Files/Porrello_Endothelial.F-CycCl_vs_F-Act.wig.txt.gz", sep = "\t", header = TRUE)
# EC_coverage.rle <- seqmonk_to_rle(endothelial_cov, col_idx = 13:28)
#
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
#
#
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
#' @param rle_input   rle input object
#' @param gtf_gr   GTF file as genomics ranges pbject
#' @param geneSymbol name of gene to interrogate
#'
rle_to_WIG <- function(rle_input, gtf_gr=gtf_gr, geneSymbol="Dnajc19")
{
  toExtract <-geneToGR(geneSymbol=geneSymbol, gtf_gr)
  tmp <- GenomicRanges::findOverlaps(FC$gr,toExtract)
  idx <- S4Vectors::queryHits(tmp)
  min_idx <- min(idx)
  max_idx <- max(idx)
  # Next need to extract idx coordinates from rle_input

}
