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
#' @param outdir directory to output the bam files. The bam files will be called [celltype].bam 
#' @param yieldSize number of lines of bam files to load. Default: 1000000 
#' 
#' @return 
#' 
#' @examples
#' 
#' extdata_path <- system.file("extdata",package = "scpolya")
#' load(paste(extdata_path,"TIP_vignette_gene_Seurat.RData",sep="/"))
#' cellbc.df <- data.frame(celltype=genes.seurat@active.ident, cellbc= names(genes.seurat@active.ident))
#' bam <- "c:/BAM/Harvey/scpolyA/one_percent.bam"
#' splitBam(bam, cellbc.df, "c:/TEMP/")
#' 
#' @export 
#' 
splitBam <- function(bam, cellbc.df, outdir, yieldSize = 1000000) { 
 # require(GenomicAlignments)
#  require(rtracklayer)
  
  message("splitting bam file: ", bam) 

   param <- Rsamtools::ScanBamParam(tag=c("CB", "UB"))
  
   ctypes <- unique(cellbc.df$celltype) 
   print(ctypes)    
   for(eachtype in ctypes) {
      message("processing cell type ", eachtype)
      outfile <- paste0(outdir, eachtype, ".bam") 
      aln.per.type <- NULL
      cellbc <- subset(cellbc.df, celltype == eachtype)$cellbc
      bamfile <- Rsamtools::BamFile(bam, yieldSize=yieldSize)
      open(bamfile)
      while (length(chunk0 <- GenomicAlignments::readGAlignments(bamfile,param=param))) {
        cat("chunk0:", length(chunk0), "length of aln: ", length(aln.per.type), "\n")
        
        if (length(aln.per.type) == 0)
          aln.per.type <- chunk0[which(S4Vectors::mcols(chunk0)$CB %in% cellbc)]
        else
          aln.per.type <- c(aln.per.type, chunk0[which(S4Vectors::mcols(chunk0)$CB %in% cellbc)])
      } # read in yieldSize number of records per iteration
      close(bamfile)
     
      if (length(aln.per.type) == 0) {
        message("No data found for ", eachtype)
      }
      else {
        message("Writing to ", outfile) 

       # as(aln.per.type, "GAlignments")
        rtracklayer::export(aln.per.type, Rsamtools::BamFile(outfile))
      }
      
   } ## Loop over cell types 

} 
