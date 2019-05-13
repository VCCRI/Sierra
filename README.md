
## scPolyA

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
find_polyA(output.file, reference.file, bamfile, junctions.file) 

```

## Counting per peak 

After you finish peak calling, you will have a file with the peak location information and we can now recount the data to create a per peak counts table. 

* polyA.sites.file: name of the output file from the `find_polyA()` function
* reference.file: file with the gene annotations 
* bamfile: bam file of your data, it is essential this bam file contains the UMI and barcode information in the tags 
* whitelist.file: file with a whitelist of the barcodes used 
* output.file: name of the file with the resulting counts data 

```{r} 
count_polyA(polyA.sites.file, reference.file, bamfile, whitelist.file, output.file) 
```

## Integrating multiple data-sets

If you have multiple data-sets, you will need to run an extra step to merge together the peaks called from the different sequencing runs so that counting is run on a unified set of peaks. To run peak merging, create a table with two columns: 'Peak_file', containing the files of peaks generated with find_polyA, and 'Identifier', containing labels for the data-sets. Specify an output file name (output.file) where the merged peaks will be written. Number of cores can be set to speed up the process. 

```{r}
peak.dataset.table = data.frame(Peak_file = c("Condition1_peaks.txt", "Condition2_peaks.txt"),
                                Identifier = c("Condition1", "Condition2"), stringsAsFactors = FALSE)

output.file = "conditions_merged_peaks.txt"

merged.peak.table = do_peak_merging(peak.dataset.table, output.file = output.file, ncores = 4)
```





