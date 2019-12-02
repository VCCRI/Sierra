
## Sierra

Sierra is an R package designed for detecting differential transcript usage analysis from polyA-enriched single cell RNASeq data. Sierra identifies coordinates of read pileups (i.e. peaks) and performs UMI counting, followed by differential usage analysis between defined cell populations. Please read vignette for a demonstration on how to use this software.

The following commands will install Sierra on your local machine:

```
install.packages("devtools")
devtools::install_github("VCCRI/Sierra", build = TRUE, build_vignettes = TRUE, build_opts = c("--no-resave-data", "--no-manual"))
```
To access the vignette:

```
library(Sierra)
browseVignettes("Sierra")
```