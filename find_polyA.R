library(GenomicRanges)
library(GenomicAlignments)
library(reshape2)
library(dplyr)
library(doParallel)

#setwd("~/Dropbox (Sydney Uni)/scPolyA/")

# Files to read in 
#output.file     <- "results/10x_neurons_9k/test.tab"
#reference.file  <- "~/ReferenceData/Mouse/mm10_GRCm38.p6_withrefseqid.csv"
#bamfile         <- "~/Documents/temp/neuron_9k_possorted_genome_bam.bam"
#junctions.file  <- "data/regtools_junctions_50.bed"

# Cutoffs to set 
#min.jcutoff      <- 50 
#min.jcutoff.prop <- 0.05 
#min.cov.cutoff   <- 500 
#min.cov.prop     <- 0.05 
#min.peak.cutoff  <- 200 
#min.peak.prop    <- 0.05 
#ncores <- 12

###################################################################
## Helper function 
makeExons <- function(x) { 
   to.write <- paste0("(", as.character(x[1]), ",")
   niters <- length(x) -1 
   for(i in 1:niters) { 
      if(x[i+1]-x[i] == 1) { 
        next 
      } else {  
        to.write <- paste0(to.write, as.character(x[i]), ")(", 
                           as.character(x[i+1]), ",")
      }
   }
   to.write <- paste0(to.write, x[niters+1], ")")
   return(to.write)
   
}


find.polyA <- function(output.file, reference.file, bamfile, junctions.file, 
                       min.jcutoff=50, min.jcutoff.prop = 0.05, min.cov.cutoff = 500, 
                       min.cov.prop = 0.05, min.peak.cutoff=200, min.peak.prop = 0.05, ncores = 1) { 

genes.ref <- read.table(reference.file, 
                        header = TRUE, sep = ",", stringsAsFactors = FALSE)
chr.names <- as.character(unique(genes.ref$chr)) 
genes.ref <- subset(genes.ref, chr %in% chr.names)
n.genes <- nrow(genes.ref)


# Initiate the output file 
write("Gene, Chr, Strand, MaxPosition, Fit.max.pos, Fit.start, Fit.end, mu, sigma, k, exon/intron, exon.pos", file = output.file)

# Read in the junction information 
junctions <- read.table(junctions.file, sep = "\t",header = FALSE)
junctions <- cbind(junctions,
                   colsplit(junctions$V11, ",", c("blocks1","blocks2")))
junctions$start <- junctions$V2+junctions$blocks1
junctions$end <- junctions$V3-junctions$blocks2
junctions.GR <- GRanges(seqnames = junctions$V1, 
                        IRanges(start = junctions$start, 
                                end = junctions$end), counts = junctions$V5)

registerDoParallel(cores=ncores)

foreach(i = 1:n.genes) %dopar% {
#for(i in 1:n.genes) { 
#for(i in c(7509)) { 
#i <- 7509
gene.name <- genes.ref[i, "Gene"]
seq.name <- genes.ref[i,"chr"]
gene.start <- genes.ref[i,"start"]
gene.end <- genes.ref[i, "end"]
strand <- genes.ref[i,"strand"]

message(i, " :", gene.name)
isMinusStrand <- if(strand==1) FALSE else TRUE
which <- GRanges(seqnames = seq.name, ranges = IRanges(gene.start, gene.end))
param <- ScanBamParam(which = which,
                      flag=scanBamFlag(isMinusStrand=isMinusStrand))

aln <- readGAlignments(bamfile, param=param)
aln_cov <- coverage(aln)[seq.name][[1]]

data <- data.frame(pos = seq(gene.start, gene.end), 
                   coverage = runValue(aln_cov)[findRun(gene.start:gene.end, aln_cov)]) 

# Find the junction which overlaps this gene 
j.cutoff <- max(min.jcutoff,min.jcutoff.prop*max(data$coverage)) 
hits <- findOverlaps(which, junctions.GR)
this.junctions.GR <- junctions.GR[hits@to]
this.junctions.GR <- subset(this.junctions.GR, counts > j.cutoff)
n.junctions <- length(this.junctions.GR)
data.no.juncs <- data

## This is pretty slow way to do this filtering, 
## can definitely improve computationally!
if(n.junctions > 0) { 
    for(i in 1:n.junctions) { 
      j.start <- start(this.junctions.GR[i])
      j.end <- end(this.junctions.GR[i])
      data.no.juncs <- data.no.juncs %>% 
      filter(pos < j.start | pos > j.end)
    }
}

## Find peaks 

totalcov <- sum(data$coverage)
cutoff <- max(min.cov.cutoff,min.cov.prop*totalcov) 
covsum <- totalcov
maxpeakval <- max(data$coverage)
maxpeakcutoff <- max(min.peak.cutoff,min.peak.prop*maxpeakval ) 

#message("Finding exonic sites...")
n.points <- nrow(data.no.juncs)
if(n.points > 0) { 
  while(covsum > cutoff) {
  maxpeak <- which.max(data.no.juncs$coverage)
  if(data.no.juncs[maxpeak, "coverage"] < maxpeakcutoff) { break }
  start <- maxpeak - 300 
  end <- maxpeak + 299
  
  if(start < 1 ) { start <- 1 }
  if(end > n.points) { end <- n.points }
  
  maxval <- max(data.no.juncs$coverage)
  #message(start, ",", end, ",", maxval)
  #message("length: ", data.no.juncs[end,"pos"] - data.no.juncs[start, "pos"])
  #message("max peak: ", data.no.juncs[maxpeak, "pos"])
    
  fit.data <- data.frame(x = seq(1,end-start+1), y = data.no.juncs[start:end,"coverage"])
  nls.res <- NULL 
  tryCatch({
  nls.res <- nls( y ~ k*exp(-1/2*(x-mu)^2/sigma^2), 
                  start=c(mu=300,sigma=100,k=maxval) , data = fit.data)
  }, error = function(err) { print(err) })
  
  if(!is.null(nls.res)) { 
     residuals <- sum(summary(nls.res)$residuals )
     v <- summary(nls.res)$parameters[,"Estimate"]
     fitted.peak <- maxpeak - 300 + floor(v[1])
     from <- fitted.peak - 3*floor(v[2])
     to <- fitted.peak + 3*floor(v[2])
     
     # Handle the cases where the peak is too close to eithr the start or the end
     if(from < 1) { from = 1 } 
     if(to > n.points) { to = n.points}
     
     if(fitted.peak <= 0) { 
       peak.pos <- "Negative" 
     } else {
       peak.pos <- data.no.juncs[fitted.peak, "pos"]
     } 
     
     isGapped <- FALSE 
     if(to <= 0) { 
       to.pos <- "Negative" 
     } else {
       to.pos <- data.no.juncs[to, "pos"]
       # Check if this is a spliced region 
       pos.gaps <- diff(data.no.juncs[from:to, "pos"])
       
       if(length(which(pos.gaps > 1) > 0)) { 
          isGapped <- TRUE 
       } 
       
     } 

     
     if(isGapped) { 
        exon.pos <- makeExons(data.no.juncs[from:to, "pos"]) 
     } else { 
        exon.pos <- "NA"   
     }
     
     line <- paste(gene.name, seq.name, strand, data.no.juncs[maxpeak, "pos"], 
                peak.pos, 
                data.no.juncs[from, "pos"], 
                to.pos, 
                v[1], v[2], v[3], "non-juncs", exon.pos, sep="\t")
     #print(line)
     write(line,file=output.file,append=TRUE)
  } else { 
    line <- paste(gene.name, seq.name, strand, data.no.juncs[maxpeak, "pos"], 
               "NA", "NA", "NA", "NA", "NA", "NA", "non-juncs", "NA", sep="\t")
    write(line,file=output.file,append=TRUE)    
  }
  
  data.no.juncs[start:end, "coverage"] <- 0
  covsum <- sum(data.no.juncs$coverage)
  #print(covsum)
}
}

## Now let's see if there are any peaks in the introns 
## test each intron separate 
#message("Finding intronic peaks...")

reduced.junctions <- reduce(this.junctions.GR) 
n.rjunctions <- length(reduced.junctions)

if(n.junctions > 0) { 

  for(i in 1:n.rjunctions) { 
    #message(i)
    j.start <- start(reduced.junctions[i])
    j.end <- end(reduced.junctions[i])
    intron.data <- data %>% 
      filter(pos > j.start & pos < j.end)
    
    if(nrow(intron.data) == 0) { next }
    maxpeak <- which.max(intron.data$coverage)
    maxval <- intron.data[maxpeak, "coverage"]
    
    #message(maxpeak, "    ", maxval)
    if(maxval < maxpeakcutoff) { next }
    fit.data <- data.frame(x = seq(1,nrow(intron.data)), 
                           y = intron.data[,"coverage"])
    nls.res <- NULL 
    tryCatch({
    nls.res <- nls( y ~ k*exp(-1/2*(x-mu)^2/sigma^2), 
                    start=c(mu=maxpeak,sigma=100,k=maxval) , data = fit.data)
    }, error = function(err) { print(err) })
    
    if(!is.null(nls.res)) { 
       residuals <- sum(summary(nls.res)$residuals )
       v <- summary(nls.res)$parameters[,"Estimate"]
       #print(v)
       fitted.peak <- floor(v[1])
       from <- fitted.peak - 3*floor(v[2])
       to <- fitted.peak + 3*floor(v[2])
       if(from < 1) { from = 1 } 
       this.n.points <- nrow(intron.data) 
       if(to > this.n.points) { to = this.n.points}
       
       if(fitted.peak <= 0) { 
         peak.pos <- "Negative" 
       } else {
         peak.pos <- intron.data[fitted.peak, "pos"]
       } 
       
       if(to <= 0) { 
         to.pos <- "Negative" 
       } else {
         to.pos <- intron.data[to, "pos"]
       } 
       
       line=paste(gene.name, seq.name, strand, intron.data[maxpeak, "pos"], 
                  peak.pos, 
                  intron.data[from, "pos"], 
                  to.pos, 
                  v[1], v[2], v[3], "junctions", "NA", sep="\t")
       
       #line=paste(gene.name, seq.name, maxpeak, v[1], v[2], v[3], "junctions", sep=",")
       #print(line)
       write(line,file=output.file,append=TRUE)
    } else { 
       line=paste(gene.name, seq.name, strand, intron.data[maxpeak, "pos"], 
                  "NA", "NA", "NA", "NA", "NA", "NA", "junction", "NA", sep="\t")
       write(line,file=output.file,append=TRUE)      
    }
  }
  
}

} # End loop for genes 

} # End function 
