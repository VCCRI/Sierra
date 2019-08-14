## Utility to split a bam file into multiple bam files based on the barcode 

library(GenomicAlignments) 
library(rtracklayer)

splitBam <- function(bam, cellbc.df, outdir, yieldSize = 1000000) { 
   message("splitting bam file: ", bam) 
  
   param <- ScanBamParam(tag=c("CB", "UB"))
  
   ctypes <- unique(cellbc.df$celltype) 
   print(ctypes)    
   for(eachtype in ctypes) {
      message("processing cell type ", eachtype)
      outfile <- paste0(outdir, eachtype, ".bam") 
      message("outputing to ", outfile)
      aln.per.type <- NULL
      cellbc <- subset(cellbc.df, celltype == eachtype)$cellbc
      bamfile <- BamFile(bam, yieldSize=10000000)
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
