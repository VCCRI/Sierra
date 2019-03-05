library(GenomicRanges)
library(GenomicAlignments)
library(reshape2)
library(dplyr)


# Files to read in 
#polyA.sites.file     <- "results/10x_neurons_9k/neurons_9k_polyA_sites_cleaned.csv"
#reference.file  <- "~/ReferenceData/Mouse/mm10_GRCm38.p6_withrefseqid.csv"
#bamfile         <- "~/Documents/temp/neuron_9k_possorted_genome_bam.bam"
#countUMI <- TRUE 
#whitelist.file <- "data/barcodes.tsv"
#output.file <- "results/10x_neurons_9k/counts_by_polyA_v2.tab"
###################################################################

count_polyA <- function(polyA.sites.file, reference.file, bamfile, whitelist.file, output.file, countUMI=TRUE) { 

whitelist.bc <- read.table(whitelist.file, stringsAsFactors = FALSE)
whitelist.bc <- whitelist.bc[,1]
n.bcs <- length(whitelist.bc) 
message("There are ", n.bcs, " whitelist barcodes.")

n.columns <- n.bcs + 1
colheadings <- c("polyAID", whitelist.bc)
write(colheadings, file = output.file,  sep = "\t", ncolumns = n.columns)
genes.ref <- read.table(reference.file, 
                        header = TRUE, sep = ",", stringsAsFactors = FALSE)
chr.names <- as.character(unique(genes.ref$chr)[2:22]) 
genes.ref <- subset(genes.ref, chr %in% chr.names)
n.genes <- nrow(genes.ref)

polyA.sites <- read.table(polyA.sites.file, header = T, sep = ",", 
                          stringsAsFactors = FALSE)

for(each.chr in chr.names) { 
#each.chr <- "19"
message("Processing chr: ", each.chr)
for(strand in c(1, -1) ) { 
message(" and strand ", strand)
#strand <- 1 
isMinusStrand <- if(strand==1) FALSE else TRUE
polyA.sites.chr <- filter(polyA.sites, Chr == each.chr & Strand == strand) %>% 
                     select(GeneID, Chr, Start, End, Strand)
polyA.sites.chr$Start <- as.integer(polyA.sites.chr$Start)
polyA.sites.chr$End <- as.integer(polyA.sites.chr$End)
polyA.sites.chr <- filter(polyA.sites.chr, Start < End)

isMinusStrand <- if(strand==1) FALSE else TRUE
which <- GRanges(seqnames = each.chr, ranges = IRanges(1, max(polyA.sites.chr$End) ))

param <- ScanBamParam(tag=c("CB", "UB"),
                      which = which,
                      flag=scanBamFlag(isMinusStrand=isMinusStrand))

aln <- readGAlignments(bamfile, param=param)

nobarcodes <- which(is.na(mcols(aln)$CB))
noUMI <- which(is.na(mcols(aln)$UB))
to.remove <- union(nobarcodes, noUMI)
aln <- aln[-to.remove]
whitelist.pos <- which(mcols(aln)$CB %in% whitelist.bc)
aln <- aln[whitelist.pos]

# For de-duplicating UMIs, let's just remove a random read 
# when there is a duplicate 
if(countUMI) { 
   mcols(aln)$CB_UB <- paste0(mcols(aln)$CB, "_", mcols(aln)$UB)
   uniqUMIs <- which(!duplicated(mcols(aln)$CB_UB))
   aln <- aln[uniqUMIs]
} 

aln <- split(aln, mcols(aln)$CB) 

polyA.GR <- GRanges(seqnames = polyA.sites.chr$Chr, 
                    IRanges(start = polyA.sites.chr$Start, 
                            end = as.integer(polyA.sites.chr$End))) 
n.polyA <- length(polyA.GR)
barcodes.gene <- names(aln)
res <- sapply(barcodes.gene, function(x) countOverlaps(polyA.GR, aln[[x]]))


mat.to.write <- matrix(0L, nrow = n.polyA, ncol = n.bcs)
mat.to.write[,match(barcodes.gene, whitelist.bc)] <- res
polyA.ids <- paste0(polyA.sites.chr$GeneID, ":", polyA.sites.chr$Chr, ":", polyA.sites.chr$Start, 
                    "-", polyA.sites.chr$End, ":", polyA.sites.chr$Strand )
rownames(mat.to.write) <- polyA.ids

write.table(mat.to.write, file = output.file, quote = F, col.names = F, row.names = T, sep = "\t", append = T)

} # Loop for strand 

} # Loop for chr 

} # End function 
