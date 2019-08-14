#############################################################################
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
#' @export 
#' 
splitBam <- function(bam, cellbc.df, outdir, yieldSize = 1000000) { 
   message("splitting bam file: ", bam) 
   require(GenomicAlignments) 
   require(rtracklayer)
  
   param <- ScanBamParam(tag=c("CB", "UB"))
  
   ctypes <- unique(cellbc.df$celltype) 
   print(ctypes)    
   for(eachtype in ctypes) {
      message("processing cell type ", eachtype)
      outfile <- paste0(outdir, eachtype, ".bam") 
      message("outputing to ", outfile)
      aln.per.type <- NULL
      cellbc <- subset(cellbc.df, celltype == eachtype)$cellbc
      bamfile <- BamFile(bam, yieldSize=yieldSize)
      open(bamfile)
      while (length(chunk0 <- readGAlignments(bamfile,param))) {
         chunk0 <- readGAlignments(bamfile,param=param)
         cat("chunk0:", length(chunk0), "length of aln: ", length(aln.per.type), "\n")
      
         if(is.null(aln.per.type)) { 
           aln.per.type <- chunk0[which(mcols(chunk0)$CB %in% cellbc)]
         } else { 
           aln.per.type <- c(aln.per.type, chunk0[which(mcols(chunk0)$CB %in% cellbc)])
         } 
         #print(class(aln.per.type))
      }
      close(bamfile)
      message("Writing to bam file") 
      rtracklayer::export(aln.per.type, BamFile(outfile))
      
   } ## Loop over cell types 

} 
