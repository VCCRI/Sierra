## ----setup,echo=FALSE------------------------------------------------------
library(knitr)
library(BiocStyle)

#Color Format
colFmt = function(x,color){
  outputFormat = knitr::opts_knit$get("rmarkdown.pandoc.to")
  if(outputFormat == 'latex')
    paste("\\textcolor{",color,"}{",x,"}",sep="")
  else if(outputFormat == 'html')
    paste("<font color='",color,"'>",x,"</font>",sep="")
  else
    x
}


## ----dimplot, out.width='100%', fig.cap = 'Dimplot output',echo=FALSE------
knitr::include_graphics('DimPlot.png')

## ----Seurat.FeaturePlot, out.width='100%', fig.cap = 'Seurat FeaturePlot  output',echo=FALSE----
knitr::include_graphics('Seurat.FeaturePlot.png')

## ----PlotRelativeExpressionTSNE, out.width='100%', fig.cap = 'PlotRelativeExpressionTSNE output',echo=FALSE----
knitr::include_graphics('PlotRelativeExpressionTSNE.png')

## ----PlotRelativeExpressionBox, out.width='100%', fig.cap = 'PlotRelativeExpressionBox output',echo=FALSE----
knitr::include_graphics('PlotRelativeExpressionBox.png')

## ----SplitBam, out.width='100%', fig.cap = 'SplitBam output',echo=FALSE----
knitr::include_graphics('PlotCoverage_CXCL12.png')

## ----sessionInfo-----------------------------------------------------------
sessionInfo()

