
##########################################################
#'
#' Calculate relative expression between two or more peaks
#'
#' Calculate a relative expression between two or more peaks by dividing
#' the expression of each peak by the mean of the peak expression for that gene -
#' or set of provided peaks
#'
#' @param peaks.object Seurat object
#' @param peak.set set of peaks
#' @param gene.name gene name for retrieving a set of peaks
#' @param feature.type features to consider. 3'UTR and exon by default.
#'
#' @return a matrix of relative expression
#'
#' @examples
#' \dontrun{
#'     GetRelativeExpression(peaks.seurat.object, gene.name = "Cxcl12")
#' }
#'
#' @export
#'
GetRelativeExpression <- function(peaks.object, peak.set = NULL, gene.name = NULL,
                                          feature.type = c("UTR3", "exon")) {

  if (class(peaks.object) == "Seurat") {
    relative.expression.data <- get_relative_expression_seurat(peaks.seurat.object = peaks.object,
                                                               peak.set = peak.set, gene.name = gene.name,
                                                               feature.type = feature.type)
    return(relative.expression.data)
  } else if (class(peaks.object) == "SingleCellExperiment") {
    relative.expression.data <- get_relative_expression_sce(peaks.sce.object = peaks.object,
                                                               peak.set = peak.set, gene.name = gene.name,
                                                               feature.type = feature.type)
    return(relative.expression.data)
  }

}


##########################################################
#'
#' Calculate relative expression between two or more peaks
#'
#' Calculate a relative expression between two or more peaks by dividing
#' the expression of each peak by the mean of the peak expression for that gene -
#' or set of provided peaks
#'
#' @param peaks.seurat.object Seurat object
#' @param peak.set set of peaks
#' @param gene.name gene name for retrieving a set of peaks
#' @param feature.type features to consider. 3'UTR and exon by default.
#'
#' @return a matrix of relative expression
#'
#' @examples
#' \dontrun{
#' get_relative_expression_seurat(peaks.seurat.object, gene.name = "Cxcl12")
#' }
get_relative_expression_seurat <- function(peaks.seurat.object, peak.set = NULL, gene.name = NULL,
                                    feature.type = c("UTR3", "exon")) {

  ## make sure either a gene or peak set has been provided
  if (is.null(gene.name) & is.null(peak.set)) {
    stop("Please provide a gene or set of peaks")
  }

  ## if no peaks are provided, use the gene name to select peaks
  if (is.null(peak.set)) {
    peak.set <- SelectGenePeaks(peaks.seurat.object, gene.name, feature.type = feature.type)
  }

  ## Check that peaks correspond to the same gene
  peak.gene.names = sub("(.*).*:.*:.*-.*:.*", "\\1", peak.set)
  if (length(unique(peak.gene.names)) > 1) {
    stop("Multiple genes detected in peak set - please ensure input peaks are from one gene")
  }

  ## access expression data for this set of peaks
  expression.data <- GetAssayData(peaks.seurat.object)[peak.set, ]

  if (length(peak.set) == 1) {
    return(expression.data)
  }

  cell.names <- colnames(peaks.seurat.object)

  population.names <- Idents(peaks.seurat.object)

  ## Calculate population-level gene-mean and relative peak expression values
  population.names <- names(table(Idents(peaks.seurat.object)))
  population.relative.usage <- c()
  gene.population.means <- c()
  for (cl in population.names) {
    cell.set <- colnames(peaks.seurat.object)[which(Idents(peaks.seurat.object) == cl)]
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
  population.names <- names(table(Idents(peaks.seurat.object)))
  for (cl in population.names) {
    cell.set <- colnames(peaks.seurat.object)[which(Idents(peaks.seurat.object) == cl)]
    expression.set <- expression.data[, cell.set]
    this.mean <- gene.population.means[cl]
    rel.values <- population.relative.usage[, cl]
    relative.exp.values <- ( (exp(expression.set) - 1) / (this.mean + 1) ) * rel.values
    relative.expression.data <- cbind(relative.expression.data, relative.exp.values)
  }
  relative.expression.data <- relative.expression.data[, colnames(peaks.seurat.object)]

  relative.expression.data <- log2(relative.expression.data + 1)
  relative.expression.data <- as(relative.expression.data, "sparseMatrix")

  return(relative.expression.data)
}

##########################################################
#'
#' Calculate relative expression between two or more peaks
#'
#' Calculate a relative expression between two or more peaks by dividing
#' the expression of each peak by the mean of the peak expression for that gene -
#' or set of provided peaks
#'
#' @param peaks.sce.object Seurat object
#' @param peak.set set of peaks
#' @param gene.name gene name for retrieving a set of peaks
#' @param feature.type features to consider. 3'UTR and exon by default.
#'
#' @return a matrix of relative expression
#'
#' @examples
#' \dontrun{
#' get_relative_expression(peaks.seurat, gene.name = "Cxcl12")
#' }
get_relative_expression_sce <- function(peaks.sce.object, peak.set = NULL, gene.name = NULL,
                                    feature.type = c("UTR3", "exon")) {

  ## make sure either a gene or peak set has been provided
  if (is.null(gene.name) & is.null(peak.set)) {
    print("Please provide a gene or set of peaks")
  }

  ## if no peaks are provided, use the gene name to select peaks
  if (is.null(peak.set)) {
    peak.set <- SelectGenePeaks(peaks.sce.object, gene.name, feature.type = feature.type)
  }

  ## Check that peaks correspond to the same gene
  peak.gene.names = sub("(.*).*:.*:.*-.*:.*", "\\1", peak.set)
  if (length(unique(peak.gene.names)) > 1) {
    stop("Multiple genes detected in peak set - please ensure input peaks are from one gene")
  }

  ## access expression data for this set of peaks
  expression.data <- peaks.sce.object@assays$data$lnorm_counts[peak.set, ]

  if (length(peak.set) == 1) {
    return(expression.data)
  }

  cell.names <- colnames(peaks.sce.object)

  ## Calculate population-level gene-mean and relative peak expression values
  population.names <- names(table(colData(peaks.sce.object)$CellIdent))
  population.relative.usage <- c()
  gene.population.means <- c()
  for (cl in population.names) {
    cell.set <- colnames(peaks.sce.object)[which(colData(peaks.sce.object)$CellIdent == cl)]
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
  population.names <- names(table(colData(peaks.sce.object)$CellIdent))
  for (cl in population.names) {
    cell.set <- colnames(peaks.sce.object)[which(colData(peaks.sce.object)$CellIdent == cl)]
    expression.set <- expression.data[, cell.set]
    this.mean <- gene.population.means[cl]
    rel.values <- population.relative.usage[, cl]
    relative.exp.values <- ( (exp(expression.set) - 1) / (this.mean + 1) ) * rel.values
    relative.expression.data <- cbind(relative.expression.data, relative.exp.values)
  }
  relative.expression.data <- relative.expression.data[, colnames(peaks.sce.object)]

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
#' @param gene_name optional plot title
#' @param peaks.use whether to print the plot to output (default: TRUE).
#' @param population.ids size of the point (default: 0.75)
#' @param return.plot whether to return the ggplot object (default: FALSE)
#' @return NULL by default. Returns a ggplot2 object if return.plot = TRUE
#' @examples
#' \dontrun{
#' do_arrow_plot(peaks.seurat.object, gene_name = Favouritegene1)
#' }
#' @import ggplot2
#'
do_arrow_plot <- function(peaks.seurat.object, gene_name, peaks.use = NULL, population.ids = NULL,
                          return.plot = FALSE) {

  if (!'gggenes' %in% rownames(x = installed.packages())) {
    stop("Please install the gggenes package (dev. version) before using this function
         (https://github.com/wilkox/gggenes)")
  }

  peak.data = subset(Tool(peaks.seurat.object, "Sierra"), Gene_name == gene_name)
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
#' @param peaks.object peak object either Seurat or SingleCellExperiment class
#' @param peaks.to.plot Set of peaks to plot
#' @param do.plot Whether to plot to output (TRUE by default)
#' @param figure.title Optional figure title
#' @param return.plot boolean of whether to return plot. Default is TRUE.
#' @param pt.size size of the points on the t-SNE plot. Default 0.5
#' @param txt.size size of text. Default 14
#'
#' @return a ggplot2 object
#'
#' @examples
#' \dontrun{
#'     PlotRelativeExpressionTSNE(peaks.seurat, this.peak.set)
#'  }
#' @import ggplot2
#'
#' @export
#'
PlotRelativeExpressionTSNE <- function(peaks.object, peaks.to.plot, do.plot=FALSE, figure.title=NULL,
                                     return.plot = TRUE, pt.size = 0.5, txt.size = 14) {

  ## Check multiple peaks have been provided
  if (length(peaks.to.plot) < 2) {
    stop("Please provide at least two peaks for plotting relative expression")
  }

  relative.exp.data <- GetRelativeExpression(peaks.object, peak.set = peaks.to.plot)

  ggData <- data.frame(Expression = as.vector(t(as.matrix(relative.exp.data))),
                       Peak = unlist(lapply(peaks.to.plot, function(x) rep(x, ncol(relative.exp.data)))))

  # Pull out the t-SNE coordinates
  if (class(peaks.object) == "Seurat") {
    peaks.object.tsne1 <- peaks.object@reductions$tsne@cell.embeddings[, 1]
    peaks.object.tsne2 <- peaks.object@reductions$tsne@cell.embeddings[, 2]
  } else if (class(peaks.object) == "SingleCellExperiment") {
    peaks.object.tsne1 <- peaks.object@reducedDims$tsne[, 1]
    peaks.object.tsne2 <- peaks.object@reducedDims$tsne[, 2]
  }

  names(peaks.object.tsne1) <- colnames(peaks.object)
  names(peaks.object.tsne2) <- colnames(peaks.object)

  ggData$tSNE_1 = rep(peaks.object.tsne1, length(peaks.to.plot))
  ggData$tSNE_2 = rep(peaks.object.tsne2, length(peaks.to.plot))

  ggData$Cell_ID = rep(colnames(peaks.object), length(peaks.to.plot))

  ggData$Peak <- factor(ggData$Peak, levels = peaks.to.plot)
  pl <- ggplot(ggData, aes(tSNE_1, tSNE_2, color=Expression)) + geom_point(size=pt.size) + xlab("t-SNE 1") + ylab("t-SNE 2") +
    scale_color_gradient2(low="#d9d9d9", mid="red", high="brown", midpoint=min(ggData$Expression) +
                            (max(ggData$Expression)-min(ggData$Expression))/2, name="") +
    theme_classic(base_size = txt.size) + theme(strip.background = element_blank()) +
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
#' @param peaks.object Seurat object
#' @param peaks.to.plot Set of peaks to plot
#' @param do.plot Whether to plot to output (TRUE by default)
#' @param figure.title Optional figure title
#' @param return.plot Boolean of whether to return plot (default TRUE)
#' @param pt.size size of the points on the t-SNE plot. Default 0.5
#' @param txt.size size of text. Default 14
#'
#' @return a ggplot2 object
#'
#' @examples
#' \dontrun{
#'    PlotRelativeExpressionUMAP(peaks.seurat, this.peak.set)
#'  }
#' @import ggplot2
#'
#' @export
#'
PlotRelativeExpressionUMAP <- function(peaks.object, peaks.to.plot, do.plot=FALSE, figure.title=NULL,
                                     return.plot = TRUE, pt.size = 0.5, txt.size = 14) {

  ## Check multiple peaks have been provided
  if (length(peaks.to.plot) < 2) {
    stop("Please provide at least two peaks for plotting relative expression")
  }

  relative.exp.data <- GetRelativeExpression(peaks.object, peak.set = peaks.to.plot)

  ggData <- data.frame(Expression = as.vector(t(as.matrix(relative.exp.data))),
                       Peak = unlist(lapply(peaks.to.plot, function(x) rep(x, ncol(relative.exp.data)))))

  # Pull out the UMAP coordinates
  if (class(peaks.object) == "Seurat") {
    peaks.object.umap1 <- peaks.object@reductions$umap@cell.embeddings[, 1]
    peaks.object.umap2 <- peaks.object@reductions$umap@cell.embeddings[, 2]
  } else if (class(peaks.object) == "SingleCellExperiment") {
    peaks.object.umap1 <- peaks.object@reducedDims$umap[, 1]
    peaks.object.umap2 <- peaks.object@reducedDims$umap[, 2]
  }

  names(peaks.object.umap1) <- colnames(peaks.object)
  names(peaks.object.umap2) <- colnames(peaks.object)

  ggData$UMAP_1 = rep(peaks.object.umap1, length(peaks.to.plot))
  ggData$UMAP_2 = rep(peaks.object.umap2, length(peaks.to.plot))

  ggData$Cell_ID = rep(colnames(peaks.object), length(peaks.to.plot))

  ggData$Peak <- factor(ggData$Peak, levels = peaks.to.plot)
  pl <- ggplot(ggData, aes(UMAP_1, UMAP_2, color=Expression)) + geom_point(size=pt.size) + xlab("UMAP 1") + ylab("UMAP 2") +
    scale_color_gradient2(low="#d9d9d9", mid="red", high="brown", midpoint=min(ggData$Expression) +
                            (max(ggData$Expression)-min(ggData$Expression))/2, name="") +
    theme_classic(base_size = txt.size) + theme(strip.background = element_blank()) +
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
#' Generate a box plot plot using relative expression
#'
#' Given two or more peaks to plot, a relative expression score and
#' generate a box plot according to cell identities
#'
#' @param peaks.object Peak object of either Seurat or SCE class
#' @param peaks.to.plot Set of peaks to plot
#' @param do.plot Whether to plot to output (TRUE by default)
#' @param figure.title Optional figure title
#' @param return.plot Boolean (default True) identifying if plot should be returned.
#' @param pt.size Size of the points on the t-SNE plot (default 0.5)
#' @param col.set col set (default NULL)
#' @param txt.size sie of text (default 14)
#'
#' @return a ggplot2 object
#'
#' @examples
#' \dontrun{
#'    PlotRelativeExpressionBox(peaks.object, this.peak.set)
#'  }
#' @import ggplot2
#'
#' @export
#'
PlotRelativeExpressionBox <- function(peaks.object, peaks.to.plot, do.plot=FALSE, figure.title=NULL,
                                      return.plot = TRUE, pt.size = 0.5, col.set = NULL, txt.size = 14) {

  ## Check multiple peaks have been provided
  if (length(peaks.to.plot) < 2) {
    stop("Please provide at least two peaks for plotting relative expression")
  }

  relative.exp.data <- GetRelativeExpression(peaks.object, peak.set = peaks.to.plot)

  ggData <- data.frame(Expression = as.vector(t(as.matrix(relative.exp.data))),
                       Peak = unlist(lapply(peaks.to.plot, function(x) rep(x, ncol(relative.exp.data)))))

  # Pull out the t-SNE coordinates
  if (class(peaks.object) == "Seurat") {
    peaks.object.tsne1 <- peaks.object@reductions$tsne@cell.embeddings[, 1]
    peaks.object.tsne2 <- peaks.object@reductions$tsne@cell.embeddings[, 2]
    cell.idents <- Idents(peaks.object)
    if (is.null(col.set)){
      col.set = scales::hue_pal()(length(table(Idents(peaks.object))))
    }
  } else if (class(peaks.object) == "SingleCellExperiment") {
    peaks.object.tsne1 <- peaks.object@reducedDims$tsne[, 1]
    peaks.object.tsne2 <- peaks.object@reducedDims$tsne[, 2]
    cell.idents <- colData(peaks.object)$CellIdent
    if (is.null(col.set)){
      col.set = scales::hue_pal()(length(table(colData(peaks.object)$CellIdent)))
    }
  }

  names(peaks.object.tsne1) <- colnames(peaks.object)
  names(peaks.object.tsne2) <- colnames(peaks.object)

  ggData$tSNE_1 = rep(peaks.object.tsne1, length(peaks.to.plot))
  ggData$tSNE_2 = rep(peaks.object.tsne2, length(peaks.to.plot))

  ggData$Cell_ID = rep(colnames(peaks.object), length(peaks.to.plot))

  ggData$Peak <- factor(ggData$Peak, levels = peaks.to.plot)

  ## Add cell population identities. Order according to order of input
  ggData$Identity <- rep(cell.idents, length(peaks.to.plot))
  ggData$Identity = factor(ggData$Identity, levels = names(table(cell.idents)))

  pl <- ggplot(ggData, aes(y=Expression, x=Identity, fill=Identity)) + geom_boxplot(colour = "black", outlier.size = 0.75) +
    ylab("Relative expression") + scale_fill_manual(values=col.set) + theme_classic(base_size = 18) +
    theme(legend.position="none", text = element_text(size = txt.size), axis.text.x = element_text(angle=90, hjust=1, vjust=0.5)) +
    theme(strip.background = element_blank()) + theme(strip.text.x = element_text(size = txt.size)) +
    xlab("") + facet_wrap(~ Peak)

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
#' Generate a violin plot plot using relative expression
#'
#' Given two or more peaks to plot, a relative expression score and
#' generate a violin plot according to cell identities
#'
#' @param peaks.object Peak object of either Seurat or SCE class
#' @param peaks.to.plot Set of peaks to plot
#' @param do.plot Whether to plot to output (TRUE by default)
#' @param figure.title Optional figure title
#' @param return.plot Boolean of whether to return plot (default TRUE)
#' @param pt.size size of the points on the t-SNE plot
#' @param col.set default NULL
#' @param txt.size size of text. Default 14
#' @param add.jitter whether to add a geom_jitter to the plot (default: TRUE)
#' @param jitter.pt.size size of point for geom_jitter (default = 0.25)
#'
#' @return a ggplot2 object
#'
#' @examples
#' \dontrun{
#'    PlotRelativeExpressionViolin(peaks.object, this.peak.set)
#' }
#' @import ggplot2
#'
#' @export
#'
PlotRelativeExpressionViolin <- function(peaks.object, peaks.to.plot, do.plot=FALSE, figure.title=NULL,
                                      return.plot = TRUE, pt.size = 0.5, col.set = NULL, txt.size = 14,
                                      add.jitter = TRUE, jitter.pt.size = 0.25) {

  ## Check multiple peaks have been provided
  if (length(peaks.to.plot) < 2) {
    stop("Please provide at least two peaks for plotting relative expression")
  }

  relative.exp.data <- GetRelativeExpression(peaks.object, peak.set = peaks.to.plot)

  ggData <- data.frame(Expression = as.vector(t(as.matrix(relative.exp.data))),
                       Peak = unlist(lapply(peaks.to.plot, function(x) rep(x, ncol(relative.exp.data)))))

  # Pull out the t-SNE coordinates
  if (class(peaks.object) == "Seurat") {
    peaks.object.tsne1 <- peaks.object@reductions$tsne@cell.embeddings[, 1]
    peaks.object.tsne2 <- peaks.object@reductions$tsne@cell.embeddings[, 2]
    cell.idents <- Idents(peaks.object)
    if (is.null(col.set)){
      col.set = scales::hue_pal()(length(table(Idents(peaks.object))))
    }
  } else if (class(peaks.object) == "SingleCellExperiment") {
    peaks.object.tsne1 <- peaks.object@reducedDims$tsne[, 1]
    peaks.object.tsne2 <- peaks.object@reducedDims$tsne[, 2]
    cell.idents <- colData(peaks.object)$CellIdent
    if (is.null(col.set)){
      col.set = scales::hue_pal()(length(table(colData(peaks.object)$CellIdent)))
    }
  }

  names(peaks.object.tsne1) <- colnames(peaks.object)
  names(peaks.object.tsne2) <- colnames(peaks.object)

  ggData$tSNE_1 = rep(peaks.object.tsne1, length(peaks.to.plot))
  ggData$tSNE_2 = rep(peaks.object.tsne2, length(peaks.to.plot))

  ggData$Cell_ID = rep(colnames(peaks.object), length(peaks.to.plot))

  ggData$Peak <- factor(ggData$Peak, levels = peaks.to.plot)

  ## Add cell population identities. Order according to order of input
  ggData$Identity <- rep(cell.idents, length(peaks.to.plot))
  ggData$Identity = factor(ggData$Identity, levels = names(table(cell.idents)))

  pl <- ggplot(ggData, aes(y=Expression, x=Identity, fill=Identity)) + geom_violin(colour = "black") +
    ylab("Relative expression") + scale_fill_manual(values=col.set) + theme_classic(base_size = 18) +
    theme(legend.position="none", text = element_text(size = txt.size), axis.text.x = element_text(angle=90, hjust=1, vjust=0.5)) +
    theme(strip.background = element_blank()) + theme(strip.text.x = element_text(size = txt.size)) +
    xlab("") + facet_wrap(~ Peak)

  if (add.jitter) {
    pl <- pl + geom_jitter(size = jitter.pt.size)
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

##########################################################
#'
#' Plot global shifts in 3'UTR length
#'
#' Plot global shifts in 3'UTR lengths between cell populations.
#' Input is a table of results from the Detect3UTRLengthShift functions.
#' By default evaluates whether there is a significant shift in 3'UTR length
#' between upregulated and downregulated peaks using the Wilcoxon Rank-sum test.
#' 
#' @param results.table table produced by the DetectUTRLengthShift function
#' @param plot.title optional title
#' @param do.ranksum.test whether to perform a ranksum test on the shift in UTR usage
#' @param return.plot whether to return the ggplot2 object
#' @param do.plot whether to print the figure to output
#' 
#' 
#' @examples
#' 
#' extdata_path <- system.file("extdata",package = "Sierra")
#' results.file <- paste0(extdata_path,"/Cycling_vs_resting_fibro_UTR_length_res.RData")
#' load(results.file)
#' 
#' PlotUTRLengthShift(res.table)
#'
#' @import ggplot2
#' @import dplyr
#' @importFrom genefilter plot
#' 
#' @export
#'
PlotUTRLengthShift <- function(results.table,
                                plot.title = "Global shift in 3'UTR length",
                                do.ranksum.test = TRUE,
                                return.plot = TRUE,
                                do.plot = FALSE) {

  locations.res.table.up <- subset(results.table, FC_direction == "Up")
  pos.upreg <- apply(as.matrix(locations.res.table.up[, c("SiteLocation","NumSites")]), 1, 
                function(x) {relative_location(x[1], x[2])})
  
  locations.res.table.down <- subset(results.table, FC_direction == "Down")
  pos.downreg <- apply(as.matrix(locations.res.table.down[, c("SiteLocation","NumSites")]), 1, 
                function(x) {relative_location(x[1], x[2])})
  
  if (do.ranksum.test) {
    this.test <- wilcox.test(pos.upreg, pos.downreg)
    print("Wilcoxon Rank-sum test comparing relative peak locations for up- vs down-regulated peaks:")
    print(paste0("P-value = ", this.test$p.value)) 
  }
  
  ggData <- data.frame(Peak_location = c(pos.upreg, pos.downreg),
                       FC_direction = c(rep("Up", length(pos.upreg)), rep("Down", length(pos.downreg))))
  ggData$FC_direction <- factor(ggData$FC_direction, levels = c("Up", "Down"))
  
  pl.density <- ggplot(ggData, aes(Peak_location, stat(count), fill = FC_direction)) + 
    geom_density(alpha = 0.8) + ylab("") + theme_void(base_size = 18) + 
    theme(axis.text = element_blank(), axis.title=element_blank(), axis.ticks = element_blank()) +
    ggtitle(plot.title) + theme(legend.position = "none", plot.title = element_text(hjust = 0.5)) +
    scale_fill_brewer(palette="Set1")
  
  pl.histogram <- ggplot(ggData, aes(Peak_location, fill = FC_direction)) + 
    geom_histogram(position = position_dodge(), colour = "black", binwidth = 0.1, alpha = 0.8) +
    theme_classic(base_size = 18) + xlab("Relative peak location") + ylab("Peak count") +
    scale_fill_brewer(palette="Set1") + theme(legend.position = "right") +
    guides(fill = guide_legend(title = "Fold-change\ndirection", title.position = "top"))
    
  
  ### Combine the plots together
  pl.combined <- cowplot::plot_grid(pl.density, pl.histogram, ncol=1, 
                                    rel_heights = c(0.3, 0.8), axis = "lr", align = "v")
  
  if (do.plot) {
    plot(pl.combined)
  }
  
  if (return.plot) {
    return(pl.combined)
  }
}

####################################################
#'
#' Given a peak position in a 3'UTR out of some n number of peaks,
#' relative to the terminating exon, calculate the relative position
#' of the query peak location on a scale of 0 to 1, where 0 indicates
#' the most proximal location and 1 indicates most distal.
#'
#' @param location location
#' @param n number of locations
#'
#'
relative_location <- function(location, n) {
  make_range <- function(x){(x-min(x))/(max(x)-min(x))}
  relative.locations <- make_range( (1:n) / n )
  this.location <- relative.locations[location]
  return(this.location)
}

#####################################################################
##
#' PlotCoverage
#'
#' @description
#' Plots read coverage across a gene for a set of BAM files and/or wig data.
#'
#' @param genome_gr : genome granges object
#' @param geneSymbol : Name of gene symbol
#' @param wig_data can be a data frame or a genomic ranges object. Must be stranded.
#' @param bamfiles : BAM filenames that are to be displayed as data tracks
#' @param wig_same_strand Display same strand or opposing strand of wig data (compared to reference gene)
#' @param genome  : genome object
#' @param bamfile.tracknames : BAM track display names. Assumed to be in same order as bamfiles.
#' @param wig_data.tracknames : WIG track display names. Assumed to be in same order as wig_data.
#' @param pdf_output : If true will create output pdf files
#' @param output_file_name : Used if pdf_output is true. Location of where files will be placed.
#' @param zoom_3UTR : If TRUE will create a second figure which will zoom in on 3'UTR.
#' @return NULL by default.
#' @examples
#' 
#' extdata_path <- system.file("extdata",package = "Sierra")
#' reference.file <- paste0(extdata_path,"/Vignette_cellranger_genes_subset.gtf")
#' gtf_gr <- rtracklayer::import(reference.file)
#' bam.files <- c(paste0(extdata_path,"/Vignette_example_TIP_mi.bam"),
#'                  paste0(extdata_path,"/Vignette_example_TIP_sham.bam"))
#' 
#' 
#' PlotCoverage(genome_gr = gtf_gr, geneSymbol = "Lrrc58", genome = "mm10", 
#'            bamfiles = bam.files, bamfile.tracknames=c("MI", "sham"))
#'
#' @import Gviz
#' @export
PlotCoverage<-function(genome_gr, geneSymbol="", wig_data=NULL, bamfiles=NULL, wig_same_strand=TRUE, 
                       genome=NULL, pdf_output = FALSE, wig_data.tracknames=NULL, bamfile.tracknames=NULL,
                       output_file_name='', zoom_3UTR=FALSE)
{
  # Check that gene_name field exists
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


  ## Load BAM files onto dtrack
  if (length(bamfiles) > 0)
  { 
    # Set track naming
    if (length(bamfile.tracknames) > 0)
    {  # Defined track names has been passed to function.
      if (length(bamfile.tracknames) == length(bamfiles))
      { names(bamfile.tracknames) <- bamfiles }
      else
      { warning("BAM track names does not match number of bam files passed. 
                Replacing with filenames.") 
        bamfile.tracknames <- bamfiles
        names(bamfile.tracknames) <- bamfiles
      }
    }
    else
    { # Default is to use bamfile names.
      bamfile.tracknames <- bamfiles
      names(bamfile.tracknames) <- bamfiles
    }

    # Extend gene window 50nt in both directions    
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
      dtrack[[length(dtrack)+1]] <- Gviz::DataTrack(gr, name=bamfile.tracknames[i], type = "histogram", genome=genome)
      
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


