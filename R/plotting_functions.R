#########################################################
#'
#' Produce an arrow plot of peak expression
#'
#' Produce an arrow plot of peak expression, utlising the gggenes package. 
#'
#' @param peaks.seurat.object a Seurat object containing t-SNE coordinates and cluster ID's in @ident slot
#' @param col.set a vector of colour codes corresponding to the number of clusters
#' @param gene_name optional plot title
#' @param peaks.use 
#' @param population.ids names of identified cell populations
#' @param return.plot whether to return the ggplot object (default: FALSE)
#' @return NULL by default. Returns a ggplot2 object if return.plot = TRUE
#' @examples
#' do_arrow_plot(peaks.seurat.object, gene_name = Favouritegene1)
#'
#' @import ggplot2
#' @export
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
#  ave.expression = ave.expression[cl.use, ]   # Ralph to check if the following line is an OK replacement
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
  
  pl <- ggplot2::ggplot(gggenesData, ggplot2::aes(xmin = start, xmax = end, y = Cluster, fill = Expression)) +
    gggenes::geom_gene_arrow() + ggplot2::ggtitle(paste0(gene_name, " peak-specific expression")) +
    ggplot2::facet_wrap(~ Cluster, scales = "free", ncol = 1) +
    gggenes::theme_genes() + ggplot2::scale_fill_gradient2(low="#d9d9d9", mid="red", high="brown", 
    midpoint=min(gggenesData$Expression) + (max(gggenesData$Expression)-min(gggenesData$Expression))/2, 
    name="Expression (log2)") + ggplot2::theme(legend.position = "bottom", legend.box = "horizontal") +
    ggplot2::guides(fill = ggplot2::guide_colourbar(barwidth = 10))
  print(pl)
  
  if (return.plot) return(pl)
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
                            (max(ggData$Expression)-min(ggData$Expression))/2, name="") + theme_bw() +
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
    col.set = scales::hue_pal()(length(table(seurat.object@ident)))
  }

  ## Boxplot
  ggData$Gene = factor(ggData$Gene, levels = geneSet)
  ggData$Cluster = factor(ggData$Cluster, levels = names(table(seurat.object@ident)))
  pl <- ggplot(ggData, aes(y=Expression, x=Cluster, fill=Cluster)) + geom_boxplot(colour = "black", outlier.size = 0.75) +
    ylab("Log2 (normalised expression + 1)") + scale_fill_manual(values=col.set) + theme_bw(base_size = 16) +
    theme(legend.position="bottom", axis.ticks.x = element_blank(), axis.text.x = element_blank(), text = element_text(size = 16)) + xlab("")

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

