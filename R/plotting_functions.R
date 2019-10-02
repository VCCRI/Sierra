##########################################################
#'
#' Calculate relative expression between two or more peaks
#'
#' Calculate a relative expression between two or more peaks by dividing
#' the expression of each peak by the mean of the peak expression for that gene -
#' or set of provided peaks
#'
#' @param seurat.object Seurat object
#' @param peak.set set of peaks
#' @param gene_name gene name for retrieving a set of peaks
#' @param feature.type features to consider. 3'UTR and exon by default.
#'
#' @return a matrix of relative expression
#'
#' @examples
#' get_relative_expression(peaks.seurat, gene_name = "Cxcl12")
#'
get_relative_expression <- function(seurat.object, peak.set = NULL, gene_name = NULL,
                                    feature.type = c("UTR3", "exon")) {

  ## make sure either a gene or peak set has been provided
  if (is.null(gene_name) & is.null(peak.set)) {
    print("Please provide a gene or set of peaks")
  }

  ## if no peaks are provided, use the gene name to select peaks
  if (is.null(peak.set)) {
    peak.set <- select_gene_polyas(seurat.object, this.gene, feature.type = feature.type)
  }

  ## access expression data for this set of peaks
  expression.data <- GetAssayData(seurat.object)[peak.set, ]

  if (length(peak.set) == 1) {
    return(expression.data)
  }

  ## Calculate population-level gene-mean and relative peak expression values
  population.names <- names(table(Idents(seurat.object)))
  population.relative.usage <- c()
  gene.population.means <- c()
  for (cl in population.names) {
    cell.set <- colnames(seurat.object)[which(Idents(seurat.object) == cl)]
    expression.set <- expression.data[, cell.set]

    ## Calculate relative usage of each peak
    peak.means <- apply(expression.set, 1, function(x){mean(exp(x) - 1)})
    if (mean(peak.means) == 0) {
      relative.usage <- peak.means
    } else {
      relative.usage <- peak.means / mean(peak.means)
    }

    population.relative.usage <- cbind(population.relative.usage, relative.usage)

    ## Calculate mean expression across the gene
    gene.mean <- mean(exp(as.matrix(expression.set)) - 1)
    gene.population.means <- append(gene.population.means, gene.mean)
  }
  colnames(population.relative.usage) <- population.names
  names(gene.population.means) <- population.names

  ### Divide peak expression for each cell by cell-type expression average
  relative.expression.data <- c()
  population.names <- names(table(Idents(seurat.object)))
  for (cl in population.names) {
    cell.set <- colnames(seurat.object)[which(Idents(seurat.object) == cl)]
    expression.set <- expression.data[, cell.set]
    this.mean <- gene.population.means[cl]
    rel.values <- population.relative.usage[, cl]
    relative.exp.values <- ( (exp(expression.set) - 1) / (this.mean + 1) ) * rel.values
    relative.expression.data <- cbind(relative.expression.data, relative.exp.values)
  }
  relative.expression.data <- relative.expression.data[, colnames(seurat.object)]

  relative.expression.data <- log2(relative.expression.data + 1)
  relative.expression.data <- as(relative.expression.data, "sparseMatrix")

  return(relative.expression.data)
}


#########################################################
#'
#' Produce an arrow plot of peak expression
#'
#' Produce an arrow plot of peak expression, utlising the gggenes package.
#'
#' @param peaks.seurat.object a Seurat object containing t-SNE coordinates and cluster ID's in @ident slot
#' @param col.set a vector of colour codes corresponding to the number of clusters
#' @param gene_name optional plot title
#' @param peaks.use whether to print the plot to output (default: TRUE).
#' @param population.ids size of the point (default: 0.75)
#' @param return.plot whether to return the ggplot object (default: FALSE)
#' @return NULL by default. Returns a ggplot2 object if return.plot = TRUE
#' @examples
#' do_arrow_plot(peaks.seurat.object, gene_name = Favouritegene1)
#'
#' @import ggplot2
#'
do_arrow_plot <- function(peaks.seurat.object, gene_name, peaks.use = NULL, population.ids = NULL,
                          return.plot = FALSE) {

  if (!'gggenes' %in% rownames(x = installed.packages())) {
    stop("Please install the gggenes package (dev. version) before using this function
         (https://github.com/wilkox/gggenes)")
  }

  peak.data = subset(Tool(peaks.seurat.object, "GeneSLICER"), Gene_name == gene_name)
  if (!is.null(peaks.use)) peak.data = subset(peak.data, rownames(peak.data) %in% peaks.use)
  n.peaks = nrow(peak.data)

  if (is.null(population.ids)) population.ids = names(table(Idents(peaks.seurat.object)))

  ave.expression = Seurat::AverageExpression(peaks.seurat.object, features = rownames(peak.data), verbose = FALSE)
  ave.expression = t(as.matrix(ave.expression$RNA))
  ave.expression = ave.expression[population.ids, ]
  ave.expression = log2(ave.expression + 1)

  peak.info = c()
  for (this.peak in rownames(peak.data)) {
    this.info = data.frame(Peak = rep(this.peak, nrow(ave.expression)),
                           start = rep(peak.data[this.peak, "start"], nrow(ave.expression)),
                           end = rep(peak.data[this.peak, "end"], nrow(ave.expression)),
                           strand = rep(peak.data[this.peak, "strand"], nrow(ave.expression)),
                           direction = rep(peak.data[this.peak, "strand"], nrow(ave.expression)))
    peak.info = rbind(peak.info, this.info)
  }
  peak.info$strand = plyr::mapvalues(peak.info$strand, from = c("+", "-"), to = c("forward", "reverse"))
  peak.info$direction = plyr::mapvalues(peak.info$direction, from = c("+", "-"), to = c("1", "-1"))

  gggenesData = data.frame(Cluster = rep(rownames(ave.expression), n.peaks),
                           Expression = as.vector(ave.expression))
  gggenesData = cbind(gggenesData, peak.info)

  pl <- ggplot(gggenesData, aes(xmin = start, xmax = end, y = Cluster, fill = Expression)) +
    gggenes::geom_gene_arrow() + ggtitle(paste0(gene_name, " peak-specific expression")) +
    facet_wrap(~ Cluster, scales = "free", ncol = 1) +
    gggenes::theme_genes() + scale_fill_gradient2(low="#d9d9d9", mid="red", high="brown",
    midpoint=min(gggenesData$Expression) + (max(gggenesData$Expression)-min(gggenesData$Expression))/2,
    name="Expression (log2)") + theme(legend.position = "bottom", legend.box = "horizontal") +
    guides(fill = guide_colourbar(barwidth = 10))
  print(pl)

  if (return.plot) return(pl)
}

##########################################################
#'
#' Generate a t-SNE plot using relative expression
#'
#' Given two or more peaks to plot, a relative expression score and
#' plot on t-SNE coordinates
#'
#' @param seurat.object Seurat object
#' @param peaks.to.plot Set of peaks to plot
#' @param do.plot Whether to plot to output (TRUE by default)
#' @param figure.title Optional figure title
#' @param pt.size size of the points on the t-SNE plot
#'
#' @return a ggplot2 object
#'
#' @examples
#' PlotRelativeExpressionTSNE(peaks.seurat, this.peak.set)
#'
#' @import ggplot2
#'
PlotRelativeExpressionTSNE <- function(seurat.object, peaks.to.plot, do.plot=TRUE, figure.title=NULL,
                                     return.plot = TRUE, pt.size = 0.5, txt.size = 14) {

  ## Check multiple peaks have been provided
  if (length(peaks.to.plot) < 2) {
    stop("Please provide at least two peaks for plotting relative expression")
  }

  relative.exp.data <- get_relative_expression(seurat.object, peak.set = peaks.to.plot)

  ggData <- data.frame(Expression = as.vector(t(as.matrix(relative.exp.data))),
                       Peak = unlist(lapply(peaks.to.plot, function(x) rep(x, ncol(relative.exp.data)))))

  # Pull out the t-SNE coordinates
  seurat.object.tsne1 <- seurat.object@reductions$tsne@cell.embeddings[, 1]
  names(seurat.object.tsne1) <- colnames(seurat.object)

  seurat.object.tsne2 <- seurat.object@reductions$tsne@cell.embeddings[, 2]
  names(seurat.object.tsne2) <- colnames(seurat.object)

  ggData$tSNE_1 = rep(seurat.object.tsne1, length(peaks.to.plot))
  ggData$tSNE_2 = rep(seurat.object.tsne2, length(peaks.to.plot))

  ggData$Cell_ID = rep(colnames(seurat.object), length(peaks.to.plot))

  ggData$Peak <- factor(ggData$Peak, levels = peaks.to.plot)
  pl <- ggplot(ggData, aes(tSNE_1, tSNE_2, color=Expression)) + geom_point(size=pt.size) + xlab("t-SNE 1") + ylab("t-SNE 2") +
    scale_color_gradient2(low="#d9d9d9", mid="red", high="brown", midpoint=min(ggData$Expression) +
                            (max(ggData$Expression)-min(ggData$Expression))/2, name="") +
    theme_classic(base_size = txt.size) + theme(panel.grid = element_blank(), strip.background = element_blank()) +
    theme(strip.text.x = element_text(size = txt.size)) + facet_wrap(~ Peak)

  if (!is.null(figure.title)) {
    pl <- pl + ggtitle(figure.title)
  }

  if (do.plot) {
    plot(pl)
  }

  if (return.plot) {
    return(pl)
  }
}

##########################################################
#'
#' Generate a t-SNE plot using relative expression
#'
#' Given two or more peaks to plot, calculate a relative expression score and
#' plot on UMAP coordinates
#'
#' @param seurat.object Seurat object
#' @param peaks.to.plot Set of peaks to plot
#' @param do.plot Whether to plot to output (TRUE by default)
#' @param figure.title Optional figure title
#' @param pt.size size of the points on the t-SNE plot
#'
#' @return a ggplot2 object
#'
#' @examples
#' PlotRelativeExpressionUMAP(peaks.seurat, this.peak.set)
#'
#' @import ggplot2
#'
PlotRelativeExpressionUMAP <- function(seurat.object, peaks.to.plot, do.plot=TRUE, figure.title=NULL,
                                     return.plot = TRUE, pt.size = 0.5, txt.size = 14) {

  ## Check multiple peaks have been provided
  if (length(peaks.to.plot) < 2) {
    stop("Please provide at least two peaks for plotting relative expression")
  }

  relative.exp.data <- get_relative_expression(seurat.object, peak.set = peaks.to.plot)

  ggData <- data.frame(Expression = as.vector(t(as.matrix(relative.exp.data))),
                       Peak = unlist(lapply(peaks.to.plot, function(x) rep(x, ncol(relative.exp.data)))))

  # Pull out the t-SNE coordinates
  seurat.object.umap1 <- seurat.object@reductions$umap@cell.embeddings[, 1]
  names(seurat.object.umap1) <- colnames(seurat.object)

  seurat.object.umap2 <- seurat.object@reductions$umap@cell.embeddings[, 2]
  names(seurat.object.umap2) <- colnames(seurat.object)

  ggData$UMAP_1 = rep(seurat.object.umap1, length(peaks.to.plot))
  ggData$UMAP_2 = rep(seurat.object.umap2, length(peaks.to.plot))

  ggData$Cell_ID = rep(colnames(seurat.object), length(peaks.to.plot))

  ggData$Peak <- factor(ggData$Peak, levels = peaks.to.plot)
  pl <- ggplot(ggData, aes(UMAP_1, UMAP_2, color=Expression)) + geom_point(size=pt.size) + xlab("UMAP 1") + ylab("UMAP 2") +
    scale_color_gradient2(low="#d9d9d9", mid="red", high="brown", midpoint=min(ggData$Expression) +
                            (max(ggData$Expression)-min(ggData$Expression))/2, name="") +
    theme_classic(base_size = txt.size) + theme(panel.grid = element_blank(), strip.background = element_blank()) +
    theme(strip.text.x = element_text(size = txt.size)) + facet_wrap(~ Peak)

  if (!is.null(figure.title)) {
    pl <- pl + ggtitle(figure.title)
  }

  if (do.plot) {
    plot(pl)
  }

  if (return.plot) {
    return(pl)
  }
}

#########################################################
#'
#' t-SNE plot of cell populations
#'
#' Based on t-SNE coordinates and cluster identities stored in a Seurat object,
#' plot t-SNE colouring cells according to cluster ID
#'
#' @param seurat.obect a Seurat object containing t-SNE coordinates and cluster ID's in @ident slot
#' @param col.set a vector of colour codes corresponding to the number of clusters
#' @param title optional plot title
#' @param do.plot whether to print the plot to output (default: TRUE).
#' @param pt.size size of the point (default: 0.75)
#' @param show.labels whether to show the cluster labels (default: FALSE)
#' @param return.plot whether to return the ggplot object (default: FALSE)
#' @param simple.theme whether to remove axis (default: FALSE)
#' @return NULL by default. Returns a ggplot2 object if return.plot = TRUE
#' @examples
#' plot_tsne(apa.seurat, show.labels = TRUE)
#'
#' @import ggplot2
#'
plot_tsne <- function(seurat.object, col.set=NULL, title=NULL, do.plot=TRUE, pt.size = 0.75, show.labels = FALSE,
                     return.plot = FALSE, simple.theme=FALSE) {
  seurat.object.names <- names(seurat.object@ident)

  # get the tSNE coordinates
  seurat.object.tsne1 <- seurat.object@dr$tsne@cell.embeddings[, 1]
  names(seurat.object.tsne1) <- names(seurat.object@ident)

  seurat.object.tsne2 <- seurat.object@dr$tsne@cell.embeddings[, 2]
  names(seurat.object.tsne2) <- names(seurat.object@ident)

  # If colors not provided use the deafult ggplot2 color scheme
  if (is.null(col.set)){
    col.set = scales::hue_pal()(length(table(seurat.object@ident)))
  }

  # create the ggplot data-frame and generate a dot plot
  ggData <- data.frame(tSNE_1=seurat.object.tsne1, tSNE_2=seurat.object.tsne2, cluster=seurat.object@ident)

  if (simple.theme == FALSE) {
    pl <- ggplot(ggData, aes(tSNE_1, tSNE_2, color=cluster)) + geom_point(size=pt.size) + xlab("t-SNE 1") + ylab("t-SNE 2") +
      scale_color_manual(values=col.set, breaks=names(table(seurat.object@ident)),
                         labels=names(table(seurat.object@ident)), name="Population") + guides(color = guide_legend(override.aes = list(size=4)))
  } else {
    pl <- ggplot(ggData, aes(tSNE_1, tSNE_2, color=cluster)) + geom_point(size=pt.size) + theme_void() +
      scale_color_manual(values=col.set, breaks=names(table(seurat.object@ident)),
                         labels=names(table(seurat.object@ident)), name="Population") + guides(color = guide_legend(override.aes = list(size=4)))
  }

  if (show.labels == TRUE) {
    ggData %>%
      dplyr::group_by(cluster) %>%
      dplyr::summarize(tSNE_1 = median(x = tSNE_1), tSNE_2 = median(x = tSNE_2)) -> centers
    centers[12, "tSNE_1"] = centers[12, "tSNE_1"] - 6

    pl <- pl + geom_text(data = centers, mapping = aes(label = cluster), size = 4.5, colour="black", fontface="bold") +
      theme(text = element_text(size = 16, family="Helvetica"))
  }


  if (!is.null(title)) {
    pl <- pl + ggtitle(title)
  }

  if (do.plot) {
    plot(pl)
  }

  if (return.plot) {
    return(pl)
  }
}


#########################################################
### Plot expression for a set of genes in a tSNE plot ###
#########################################################
plot_expression_tsne <- function(seurat.object, geneSet, do.plot=TRUE, figure.title=NULL,
                               return.plot = FALSE, pt.size = 0.75) {

  # Get data-frame containing expression, gene name and cluster
  ggData = getMultiGeneExpressionData(seurat.object, geneSet)

  # Pull out the t-SNE coordinates
  seurat.object.tsne1 <- seurat.object@reductions$tsne@cell.embeddings[, 1]
  names(seurat.object.tsne1) <- colnames(seurat.object)

  seurat.object.tsne2 <- seurat.object@reductions$tsne@cell.embeddings[, 2]
  names(seurat.object.tsne2) <- colnames(seurat.object)

  ggData$tSNE_1 = rep(seurat.object.tsne1, length(geneSet))
  ggData$tSNE_2 = rep(seurat.object.tsne2, length(geneSet))

  ggData$Cell_ID = rep(colnames(seurat.object), length(geneSet))

  ggData$Gene <- factor(ggData$Gene, levels = geneSet)
  pl <- ggplot(ggData, aes(tSNE_1, tSNE_2, color=Expression)) + geom_point(size=pt.size) + xlab("t-SNE 1") + ylab("t-SNE 2") +
    scale_color_gradient2(low="#d9d9d9", mid="red", high="brown", midpoint=min(ggData$Expression) +
                            (max(ggData$Expression)-min(ggData$Expression))/2, name="") +
    theme_bw(base_size = 14) + theme(panel.grid = element_blank()) +
    theme(strip.text.x = element_text(size = 14))
  if (length(geneSet) > 1) {
    pl <- pl + facet_wrap(~Gene)
  }

  if (is.null(figure.title)) {
    if (length(geneSet) == 1) {
      pl <- pl + ggtitle(paste0(geneSet, " expression vizualised on t-SNE coordinates"))
    } else{
      pl <- pl + ggtitle("Gene expression vizualsed on t-SNE coordinates")
    }
  } else {
    pl <- pl + ggtitle(figure.title)
  }

  if (do.plot) {
    plot(pl)
  }

  if (return.plot) {
    return(pl)
  }
}


#########################################
### Do a boxplot for a panel of genes ###
#########################################
do_box_plot <- function(seurat.object, geneSet, figure.title = NULL, num_col = NULL, do.plot = TRUE, col.set = NULL,
                      return.plot = FALSE) {
  ggData = getMultiGeneExpressionData(seurat.object, geneSet)

  if (is.null(col.set)){
    col.set = scales::hue_pal()(length(table(Idents(seurat.object))))
  }

  ## Boxplot
  ggData$Gene = factor(ggData$Gene, levels = geneSet)
  ggData$Cluster = factor(ggData$Cluster, levels = names(table(Idents(seurat.object))))
  pl <- ggplot(ggData, aes(y=Expression, x=Cluster, fill=Cluster)) + geom_boxplot(colour = "black", outlier.size = 0.75) +
    ylab("Log2 (normalised expression + 1)") + scale_fill_manual(values=col.set) + theme_bw(base_size = 18) +
    theme(legend.position="none", text = element_text(size = 18), axis.text.x = element_text(angle=90, hjust=1, vjust=0.5)) +
    xlab("")

  if (!is.null(num_col)) {
    pl <- pl + facet_wrap(~ Gene, ncol = num_col)
  } else{
    pl <- pl + facet_wrap(~ Gene)
  }

  if (!is.null(figure.title)) {
    pl <- pl + ggtitle(figure.title)
  }

  if (do.plot) {
    plot(pl)
  }

  if (return.plot) {
    return(pl)
  }
}

################################################################
### Generate a data-frame of expression values for multiple  ###
### genes for input to ggplot2. Tracks cluster identities.   ###
################################################################
getMultiGeneExpressionData <- function(seurat.object, geneSet, use.log10 = FALSE) {
  expression <- c()
  geneName <- c()
  cluster <- c()
  if (use.log10){
    log.fun = log10
  } else {
    log.fun = log2
  }
  # Iterate through the genes and fill the vectors
  for (gene in geneSet) {
    expression <- append(expression, log.fun(exp(GetAssayData(seurat.object)[gene, ])))
    geneName <- append(geneName, rep(gene, length(GetAssayData(seurat.object)[gene, ])))
    cluster <- append(cluster, as.character(Idents(seurat.object)))
  }
  #create a data-frame for ggplot and return
  ggData <- data.frame(Expression=expression, Gene=geneName, Cluster=cluster)
  return(ggData)
}




#####################################################################
##
#' plotCoverage
#'
#'
#'
#' @param wig_data can be a data frame or a genomic ranges object. Must be stranded.
#' @param bamfiles : BAM filenames that are to be displayed as data tracks
#' @param wig_same_strand Display same strand or opposing strand of wig data (compared to reference gene)
#' @param pdf_output : If true will create output pdf files
#' @param output_file_name : Used if pdf_output is true. Location of where files will be placed.
#' @param zoom_UTR : If TRUE will create a second figure which will zoom in on 3'UTR.
#' @return NULL by default.
#' @examples
#'
#' gtf_file <- "u:/Reference/mm10/cellranger_genes.gtf.gz"
#' gtf_gr <- rtracklayer::import(gtf_file)
#'
#' endothelial_cov <- read.table(file="c:/BAM/Harvey/scpolyA/Porrello_Support_Files/Porrello_Endothelial.F-CycCl_vs_F-Act.wig.txt.gz", sep = "\t", header = TRUE)
#' df <- endothelial_cov[,c("Chromosome","Start","End","Probe.Strand", "MIP1_1.BAM")]
#' df <- cbind(endothelial_cov[,c("Chromosome","Start","End","Probe.Strand")],endothelial_cov[, 13:15])
#' colnames(df)[1:4] <- c("chrom", "start","end", "strand")
#' wig_data<- GenomicRanges::makeGRangesFromDataFrame(df,keep.extra.columns=TRUE)
#'  plotCoverage(genome_gr=gtf_gr, geneSymbol="Prkar1a", wig_data= wig_data)
#'
#' plotCoverage(genome_gr=gtf_gr, geneSymbol="Dnajc19", bamfiles = "c:/TEMP/Bams/F-SH.Dnajc19.bam",wig_data= wig_data)
#'
#' @import Gviz
#' @export
plotCoverage<-function(genome_gr, geneSymbol="", wig_data=NULL, bamfiles=NULL, wig_same_strand=TRUE, genome=NULL, pdf_output = FALSE, 

                       output_file_name='', zoom_3UTR=FALSE)
{
  # Need check that gene_name field exists
  GenomeInfoDb::seqlevelsStyle(genome_gr) <- "UCSC"
  idx <-which(genome_gr$gene_name == geneSymbol)
  if (length(idx) == 0)
  { warning("Could not find gene name. Please check spelling (and case)")
    return(NULL)
  }
  # Work out the genomic range to extract from
  genome_gr <- genome_gr[idx]
  start <- min(IRanges::start(IRanges::ranges(genome_gr)))
  end <- max(IRanges::end(IRanges::ranges(genome_gr)))
  # should I check that all returned chromosomes and strands are the same? They should be the same
  # Currently just grabbing first entry
  chrom <- as.character(GenomicRanges::seqnames(genome_gr))[1]
  gene_strand <- as.character(BiocGenerics::strand(genome_gr))[1]
  toExtract_gr <- GenomicRanges::GRanges(seqnames=chrom, ranges=IRanges::IRanges(start-1 , width=end-start+3), strand=gene_strand)

  # Assemble gene annotation track
  gene_gr <- IRanges::subsetByOverlaps(genome_gr, toExtract_gr)
  GenomeInfoDb::seqlevelsStyle(gene_gr) <- "UCSC"
  GenomeInfoDb::seqlevels(gene_gr) <- chrom

  # Following lines is a hack to get gene symbol to be name of transcripts.
  gene_name_idx <- which(names(GenomicRanges::elementMetadata(gene_gr)) == "gene_name")
  gene_id_idx <- which(names(GenomicRanges::elementMetadata(gene_gr)) == "gene_id")
  names(GenomicRanges::elementMetadata(gene_gr))[gene_id_idx] <- "ensemble_id"
  names(GenomicRanges::elementMetadata(gene_gr))[gene_name_idx] <- "gene_id"

  gene_txdb <- GenomicFeatures::makeTxDbFromGRanges(gene_gr)

  gtrack <- Gviz::GeneRegionTrack(gene_txdb, start = start, end = end, chromosome=chrom, name= geneSymbol)

  ##### Assemble data track(s)
  dtrack <- list()  # Add data tracks assembled on this list
  
  ## Assemble wig data tracks
  wig_tracks <- list()
  if (! is.null(wig_data))
  {
    if (typeof(wig_data) != "S4")  # assume a dataframe or list which we can create several df.
    { nc <- ncol(wig_data)  # 4 columns are chrom, start, end, strand. Thereafter are sample data
      wig_data <-  GenomicRanges::makeGRangesFromDataFrame(wig_data,keep.extra.columns=TRUE)
    }
    GenomeInfoDb::seqlevelsStyle(wig_data) <- "UCSC"
  
    if (! wig_same_strand)
    {  toExtract_gr <- GenomicRanges::invertStrand(toExtract_gr)  }
    dtrack_gr <- IRanges::subsetByOverlaps(wig_data, toExtract_gr)
  
  
    GenomeInfoDb::seqlevels(dtrack_gr) <- chrom
    sample_col_idx <- 1: ncol(S4Vectors::mcols(wig_data))
    
    # Now assemble coverage plots
    for(i in sample_col_idx) 
    {
      tmp_gr <- dtrack_gr
      S4Vectors::mcols(tmp_gr) <- S4Vectors::mcols(tmp_gr)[i]
      dtrack_name <- names(S4Vectors::mcols(tmp_gr))
      wig_tracks[[length(wig_tracks)+1]] <- Gviz::DataTrack(tmp_gr, name=dtrack_name, type = "histogram", genome=genome) 
    }
    
  }

  
  ## First load any BAM files onto dtrack
  if (length(bamfiles) > 0)
  {
    toExtract_gr <- GenomicRanges::GRanges(seqnames=chrom, ranges=IRanges::IRanges(start-50 , width=end-start+50), strand=gene_strand)


    for(i in bamfiles)
    {
      bamHeader <- Rsamtools::scanBamHeader(i)
      if (length(grep(pattern = chrom,x = names(bamHeader[[i]]$targets))) == 0)
      {   GenomeInfoDb::seqlevelsStyle(toExtract_gr) <- "NCBI" }
      else
      {   GenomeInfoDb::seqlevelsStyle(toExtract_gr) <- "UCSC" }

      param <- Rsamtools::ScanBamParam(which = toExtract_gr)
      bf <-Rsamtools::BamFile(i)

      open(bf)
      chunk0 <- GenomicAlignments::readGAlignments(bf,param=param)
      GenomeInfoDb::seqlevelsStyle(chunk0) <- "UCSC"
      close(bf)
      idx <- which(as.character(BiocGenerics::strand(chunk0)) == gene_strand)
      if (length(idx) == 0)
      {        next; }
      tmp <-GenomicRanges::coverage(chunk0[idx])

      gr <- GenomicRanges::GRanges(seqnames=chrom, ranges=IRanges::IRanges(start:end, width=1), strand=gene_strand)
      S4Vectors::mcols(gr) <- as.numeric(tmp[[chrom]])[start:end]
      dtrack[[length(dtrack)+1]] <- Gviz::DataTrack(gr, name=i, type = "histogram", genome=genome)
    }
  }

  
  if (length(wig_tracks) > 0)
  {
    dtrack <- c(wig_tracks, dtrack)
  }

  toPlot <- c(gtrack, dtrack)

  if (pdf_output)
  {
    if (output_file_name == '')
    { warning("No file name provided")
      pdf_output = FALSE
    }
    else
    {  pdf(file=output_file_name,width = 24,height = 18)
    }
  }
  Gviz::plotTracks(toPlot, from = start, to = end, chromosome= chrom, transcriptAnnotation = "gene")

  if (zoom_3UTR)
  {
    idx <- which(genome_gr$type == 'three_prime_utr')
    # Work out the genomic range of UTR
    start <- min(IRanges::start(IRanges::ranges(genome_gr[idx])))
    end <- max(IRanges::end(IRanges::ranges(genome_gr[idx])))
    Gviz::plotTracks(toPlot, from = start, to = end, chromosome= chrom, transcriptAnnotation = "gene")
  }



  if (pdf_output)
  {  dev.off() }

}


