
## sierra

SIngle cEll diffeRential gene-paRt Analysis

An R package designed to analyse of read pileups (i.e. peaks) from single cell RNASeq data. Please read vignette for a demonstration on how to use this software.

The following commands will install Sierra on your local machine:

```
install.packages("devtools")
library(devtools)
devtools::install_github("VCCRI/Sierra", build = TRUE, build_opts = c("--no-resave-data", "--no-manual"))
```