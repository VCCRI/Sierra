

###################################################################
#'
#' count polyA sites in cells
#'
#' Generates a count matrix
#'
#' @param polyA.sites.file a file containing polyA sites
#' @param reference.file reference (GTF) file
#' @param bamfile scRNA-seq BAM file
#' @param whitelist.file white list file
#' @param output.file name of output
#' @param countUMI whether to count UMIs (default: TRUE)
#' @return NULL. Writes counts to file.
#' @examples
#' count_polyA(polyA.sites.file, reference.file, bamfile, whitelist.file, output.file)
count_polyA <- function(polyA.sites.file, reference.file, bamfile, whitelist.file, output.file, countUMI=TRUE,
			ncores = 1) {

  lock <- tempfile()
  whitelist.bc <- read.table(whitelist.file, stringsAsFactors = FALSE)
  whitelist.bc <- whitelist.bc[,1]
  n.bcs <- length(whitelist.bc)
  message("There are ", n.bcs, " whitelist barcodes.")

  n.columns <- n.bcs + 1
  colheadings <- c("polyAID", whitelist.bc)
  write(colheadings, file = output.file,  sep = "\t", ncolumns = n.columns)
  genes.ref <- read.table(reference.file,
                          header = TRUE, sep = ",", stringsAsFactors = FALSE)
  chr.names <- as.character(unique(genes.ref$chr))
  genes.ref <- subset(genes.ref, chr %in% chr.names)
  n.genes <- nrow(genes.ref)

  polyA.sites <- read.table(polyA.sites.file, header = T, sep = "\t",
                            stringsAsFactors = FALSE)

  # Filter the polyA sites
  n.total.sites <- nrow(polyA.sites)
  print(head(polyA.sites))
  to.filter <- which(polyA.sites$Fit.max.pos == "Negative")
  to.filter <- union(to.filter, which(polyA.sites$Fit.start == "Negative"))
  to.filter <- union(to.filter, which(polyA.sites$Fit.end == "Negative"))
  to.filter <- union(to.filter,  which(is.na(polyA.sites$Fit.max.pos)))

  polyA.sites <- polyA.sites[-to.filter,]
  print(head(polyA.sites))
  n.filt.sites <- nrow(polyA.sites)
  message("There are ", n.total.sites, " unfiltered sites and ", n.filt.sites, " filtered sites")
  message("Doing counting for each filtered site...")

  doParallel::registerDoParallel(cores=ncores)

  print(chr.names)
  foreach::foreach(each.chr = chr.names, .packages = c("GenomicRanges", "GenomicAlignments", "Rsamtools", "dplyr")) %dopar% {
    #for(each.chr in chr.names) {

      message("Processing chr: ", each.chr)
      for(strand in c(1, -1) ) {
      message(" and strand ", strand)
      isMinusStrand <- if(strand==1) FALSE else TRUE
      polyA.sites.chr <- dplyr::filter(polyA.sites, Chr == each.chr & Strand == strand) %>%
                           dplyr::select(Gene, Chr, Fit.start, Fit.end, Strand)
      polyA.sites.chr$Fit.start <- as.integer(polyA.sites.chr$Fit.start)
      polyA.sites.chr$Fit.end <- as.integer(polyA.sites.chr$Fit.end)
      polyA.sites.chr <- dplyr::filter(polyA.sites.chr, Fit.start < Fit.end)

      isMinusStrand <- if(strand==1) FALSE else TRUE
      which <- GenomicRanges::GRanges(seqnames = each.chr, ranges = IRanges::IRanges(1, max(polyA.sites.chr$Fit.end) ))

      param <- Rsamtools::ScanBamParam(tag=c("CB", "UB"),
                            which = which,
                            flag=scanBamFlag(isMinusStrand=isMinusStrand))

      aln <- GenomicAlignments::readGAlignments(bamfile, param=param)

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

      aln <- GenomicRanges::split(aln, mcols(aln)$CB)

      polyA.GR <- GenomicRanges::GRanges(seqnames = polyA.sites.chr$Chr,
                          IRanges::IRanges(start = polyA.sites.chr$Fit.start,
                                  end = as.integer(polyA.sites.chr$Fit.end)))
      n.polyA <- length(polyA.GR)
      barcodes.gene <- names(aln)
      res <- sapply(barcodes.gene, function(x) GenomicRanges::countOverlaps(polyA.GR, aln[[x]]))


      mat.to.write <- matrix(0L, nrow = n.polyA, ncol = n.bcs)
      mat.to.write[,match(barcodes.gene, whitelist.bc)] <- res
      polyA.ids <- paste0(polyA.sites.chr$Gene, ":", polyA.sites.chr$Chr, ":", polyA.sites.chr$Fit.start,
                          "-", polyA.sites.chr$Fit.end, ":", polyA.sites.chr$Strand )
      rownames(mat.to.write) <- polyA.ids

      locked <- flock::lock(lock)
      write.table(mat.to.write, file = output.file, quote = F, col.names = F, row.names = T, sep = "\t", append = T)
      flock::unlock(locked)

    } # Loop for strand

  } # Loop for chr

} # End function

###################################################################
#'
#' Helper function
#'
#' Helper function
#'
#' @param x x
#' @return to write
#' @examples
#' makeExons(x)
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

###################################################################
#'
#' Do peak calling on a scRNA-seq BAM file
#'
#' Do peak calling on a scRNA-seq BAM file...
#'
#' @param output.file a file containing polyA sites
#' @param reference.file reference (GTF) file
#' @param bamfile scRNA-seq BAM file
#' @param junctions.file white list file
#' @param min.jcutoff name of output
#' @param min.jcutoff.prop whether to count UMIs (default: TRUE)
#' @param min.cov.cutoff min.cov.cutoff
#' @param min.cov.prop min.cov.prop
#' @param min.peak.cutoff min.peak.cutoff
#' @param min.peak.prop min.peak.prop
#' @param ncores number of cores to use
#' @return NULL. Writes counts to file.
#' @examples
#' find_polyA(output.file, reference.file, bamfile, junctions.file)
#'
#' @importFrom magrittr "%>%"
#' @importFrom foreach "%dopar%"
#'
find_polyA <- function(output.file, reference.file, bamfile, junctions.file,
                       min.jcutoff=50, min.jcutoff.prop = 0.05, min.cov.cutoff = 500,
                       min.cov.prop = 0.05, min.peak.cutoff=200, min.peak.prop = 0.05, ncores = 1) {

  lock <- tempfile()
  genes.ref <- read.table(reference.file,
                          header = TRUE, sep = ",", stringsAsFactors = FALSE)
  chr.names <- as.character(unique(genes.ref$chr))
  genes.ref <- subset(genes.ref, chr %in% chr.names)
  n.genes <- nrow(genes.ref)


  # Initiate the output file
  write("Gene\tChr\tStrand\tMaxPosition\tFit.max.pos\tFit.start\tFit.end\tmu\tsigma\tk\texon/intron\texon.pos", file = output.file)

  # Read in the junction information
  junctions <- read.table(junctions.file, sep = "\t",header = FALSE)
  junctions <- cbind(junctions,
                     reshape2::colsplit(junctions$V11, ",", c("blocks1","blocks2")))
  junctions$start <- junctions$V2+junctions$blocks1
  junctions$end <- junctions$V3-junctions$blocks2
  junctions.GR <- GenomicRanges::GRanges(seqnames = junctions$V1,
                          IRanges::IRanges(start = junctions$start,
                                  end = junctions$end), counts = junctions$V5)

  doParallel::registerDoParallel(cores=ncores)

  foreach::foreach(i = 1:n.genes) %dopar% {
    gene.name <- genes.ref[i, "Gene"]
    seq.name <- genes.ref[i,"chr"]
    gene.start <- genes.ref[i,"start"]
    gene.end <- genes.ref[i, "end"]
    strand <- genes.ref[i,"strand"]

    #message(i, " :", gene.name)
    isMinusStrand <- if(strand==1) FALSE else TRUE
    which <- GenomicRanges::GRanges(seqnames = seq.name, ranges = IRanges::IRanges(gene.start, gene.end))
    param <- Rsamtools::ScanBamParam(which = which,
                          flag=Rsamtools::scanBamFlag(isMinusStrand=isMinusStrand))

    aln <- GenomicAlignments::readGAlignments(bamfile, param=param)
    aln_cov <- GenomicRanges::coverage(aln)[seq.name][[1]]

    data <- data.frame(pos = seq(gene.start, gene.end),
                       coverage = S4Vectors::runValue(aln_cov)[S4Vectors::findRun(gene.start:gene.end, aln_cov)])

    # Find the junction which overlaps this gene
    j.cutoff <- max(min.jcutoff,min.jcutoff.prop*max(data$coverage))
    hits <- GenomicRanges::findOverlaps(which, junctions.GR)
    this.junctions.GR <- junctions.GR[hits@to]
    this.junctions.GR <- IRanges::subset(this.junctions.GR, counts > j.cutoff)
    n.junctions <- length(this.junctions.GR)
    data.no.juncs <- data

    ## This is pretty slow way to do this filtering,
    ## can definitely improve computationally!
    if(n.junctions > 0) {
      for(i in 1:n.junctions) {
        j.start <- IRanges::start(this.junctions.GR[i])
        j.end <- IRanges::end(this.junctions.GR[i])
        data.no.juncs <- data.no.juncs %>%
          dplyr::filter(pos < j.start | pos > j.end)
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
        }, error = function(err) { })

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
	  locked <- flock::lock(lock)
          write(line,file=output.file,append=TRUE)
	  flock::unlock(locked)
        } else {
          line <- paste(gene.name, seq.name, strand, data.no.juncs[maxpeak, "pos"],
                        "NA", "NA", "NA", "NA", "NA", "NA", "non-juncs", "NA", sep="\t")
	  locked <- flock::lock(lock)
          write(line,file=output.file,append=TRUE)
	  flock::unlock(locked)
        }

        data.no.juncs[start:end, "coverage"] <- 0
        covsum <- sum(data.no.juncs$coverage)
        #print(covsum)
      }
    }

    ## Now let's see if there are any peaks in the introns
    ## test each intron separate
    #message("Finding intronic peaks...")

    reduced.junctions <- GenomicRanges::reduce(this.junctions.GR)
    n.rjunctions <- length(reduced.junctions)

    if(n.junctions > 0) {

      for(i in 1:n.rjunctions) {
        #message(i)
        j.start <- IRanges::start(reduced.junctions[i])
        j.end <- IRanges::end(reduced.junctions[i])
        intron.data <- data %>%
          dplyr::filter(pos > j.start & pos < j.end)

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
        }, error = function(err) { })

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
	  locked <- flock::lock(lock)
          write(line,file=output.file,append=TRUE)
	  flock::unlock(locked)
        } else {
          line=paste(gene.name, seq.name, strand, intron.data[maxpeak, "pos"],
                     "NA", "NA", "NA", "NA", "NA", "NA", "junction", "NA", sep="\t")
	  flocked <- flock::lock(lock)
          write(line,file=output.file,append=TRUE)
	  flock::unlock(locked)
        }
      }

    }

  } # End loop for genes

} # End function
