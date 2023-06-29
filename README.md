
# Sierra

Sierra is an R package designed for detecting differential transcript usage from polyA-captured single cell (sc)RNA-seq data. Sierra identifies coordinates of read pileups (i.e. peaks) and performs UMI counting, followed by differential usage analysis between defined cell populations. Please read the vignette for a demonstration on how to use this software.

## Installation

To install Sierra on your local machine, open an R session and use the following commands:

```
install.packages("devtools")
devtools::install_github("VCCRI/Sierra", build = TRUE, build_vignettes = TRUE, build_opts = c("--no-resave-data", "--no-manual"))
```
The vignette contains a detailed walk-through of how to use Sierra. The vignette can be access from R using:

```
library(Sierra)
browseVignettes("Sierra")
```

## Alternative installation

If you have trouble building the vignette, an alternative is to view it on the [wiki](https://github.com/VCCRI/Sierra/wiki/Sierra-Vignette) and install the R package without it:

```
devtools::install_github("VCCRI/Sierra", build = TRUE)
```

## Method overview

The manuscript describing Sierra, including validation and example applications, is published in [Genome Biology](https://genomebiology.biomedcentral.com/articles/10.1186/s13059-020-02071-7).

Briefly, the Sierra pipeline requires as input a BAM file, such as produced by the 10x Genomics CellRanger software, the reference GTF file used for alignment and a BED file of junctions derived from the BAM file as produced by [RegTools](https://regtools.readthedocs.io/en/latest/).

1. Splice-aware peak calling is used to identify read pileups corresponding to potential polyA sites in the dataset.

2. UMI counting is performed against a set of peak coordinates, first unified in the case of multiple experiments. 

3. Peak coordinates are annotated with various features, including the genomic features they fall on.

4. Differential transcript usage (DTU) is evaluated between defined cell populations by applying the differential exon method DEXSeq to test for differences in the relative usage of peaks, with pseudo-bulk profiles of cells used to define replicates.

5. DTU genes can be visualised using read coverage plots, or by plotting the relative expression of peaks. 

## Citation

If you find Sierra useful in your research, please cite the following paper:

Patrick, R., Humphreys, D.T., Janbandhu, V., Oshlack A., Ho J.W.K., Harvey, R.P. and Lo K.K. Sierra: discovery of differential transcript usage from polyA-captured single-cell RNA-seq data. Genome Biol 21, 167 (2020). https://doi.org/10.1186/s13059-020-02071-7

## Contact

Sierra is maintained by Ralph Patrick, David Humphreys and Kitty Lo. For questions or feedback you can contact them at:

* Ralph Patrick: ralph.patrick at imb dot uq dot edu dot au

* David Humphreys: d.humphreys at victorchang dot edu dot au

* Kitty Lo: kitty.lo at gmail dot com

## Additional implementation

For a command line version of Sierra that can be used, for example, in a high performance computing environment, see [this implementation](https://github.com/GeertvanGeest/Sierra-commands) by Geert van Geest. 



