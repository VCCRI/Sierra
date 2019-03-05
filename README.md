# scpolya

The scPolyA package allows you to analyse any single cell RNAseq dataset that uses a 3' protocol. 

## Required inputs 

1. Bam file 
We require as input a bam file, ideally one that was created from CellRanger. scPolyA uses the CellRanger formatted tags that contains the UMI and cell barcode information in the bam file for doing the counting. 

2. Junctions file 
scPolya also requires a junctions file as a mandatory input into the peak caller. This file can be created using the `regtools` program. 

3. Reference annotation file 

4. Whitelist of cell barcodes 

## Required packages 

The following R packages are required to run scPolyA: 
* GenomicRanges
* GenomicAlignments 
* reshape2 
* dplyr 
* doParallel 


## Peak calling 

* output.file: name of the output file which will contain the location of the peaks 
* reference.file: file with the gene annotations 
* bamfile: bam file of your data, it is essential this bam file contains the UMI and barcode information in the tags 
* junctions.file: file made with `regtools` using the same bamfile as your input 

```{r}
find.polyA(output.file, reference.file, bamfile, junctions.file) 

```

## Counting per peak 

After you finish peak calling, you will have a file with the peak location information and we can now recount the data to create a per peak counts table. 

* polyA.sites.file: name of the output file from the `find.poly()` function
* reference.file: file with the gene annotations 
* bamfile: bam file of your data, it is essential this bam file contains the UMI and barcode information in the tags 
* whitelist.file: file with a whitelist of the barcodes used 
* output.file: name of the file with the resulting counts data 

```{r} 
count_polyA(polyA.sites.file, reference.file, bamfile, whitelist.file, output.file) 
```
