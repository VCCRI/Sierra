library(Sierra)

tmp_lib_check <- function() {
  stopifnot(identical(Sierra:::validate_strand_mode(0), 0L))
  stopifnot(identical(Sierra:::validate_strand_mode(1L), 1L))
  stopifnot(identical(Sierra:::validate_strand_mode(2), 2L))

  for (bad in list(-1, 3, 1.5, NA, NA_integer_, "reverse", c(1, 2), NULL)) {
    err <- tryCatch(
      {
        Sierra:::validate_strand_mode(bad)
        NULL
      },
      error = function(e) e
    )
    stopifnot(!is.null(err))
  }
}

test_alignment_helpers <- function() {
  gr <- GenomicRanges::GRanges(
    seqnames = c("chr1", "chr1", "chr1"),
    ranges = IRanges::IRanges(c(10, 20, 30), width = 10),
    strand = c("+", "-", "+")
  )
  aln <- methods::as(gr, "GAlignments")

  stopifnot(length(Sierra:::sierra_filter_alignments_by_strand(aln, "+", strandMode = 0L)) == 3L)
  stopifnot(length(Sierra:::sierra_filter_alignments_by_strand(aln, "+", strandMode = 1L)) == 2L)
  stopifnot(length(Sierra:::sierra_filter_alignments_by_strand(aln, "+", strandMode = 2L)) == 1L)
  stopifnot(length(Sierra:::sierra_filter_alignments_by_strand(aln, "-", strandMode = 1L)) == 1L)
  stopifnot(length(Sierra:::sierra_filter_alignments_by_strand(aln, "-", strandMode = 2L)) == 2L)

  f_unstranded <- Sierra:::sierra_scan_bam_flag(1, strandMode = 0L)
  f_fwd_plus <- Sierra:::sierra_scan_bam_flag(1, strandMode = 1L)
  f_fwd_minus <- Sierra:::sierra_scan_bam_flag(-1, strandMode = 1L)
  f_rev_plus <- Sierra:::sierra_scan_bam_flag(1, strandMode = 2L)
  f_rev_minus <- Sierra:::sierra_scan_bam_flag(-1, strandMode = 2L)
  f_factor_minus <- Sierra:::sierra_scan_bam_flag(factor("-1", levels = c("1", "-1")), strandMode = 1L)
  f_factor_rev_minus <- Sierra:::sierra_scan_bam_flag(factor("-1", levels = c("1", "-1")), strandMode = 2L)

  stopifnot(identical(f_unstranded, Rsamtools::scanBamFlag()))
  stopifnot(identical(f_fwd_plus, Rsamtools::scanBamFlag(isMinusStrand = FALSE)))
  stopifnot(identical(f_fwd_minus, Rsamtools::scanBamFlag(isMinusStrand = TRUE)))
  stopifnot(identical(f_rev_plus, Rsamtools::scanBamFlag(isMinusStrand = TRUE)))
  stopifnot(identical(f_rev_minus, Rsamtools::scanBamFlag(isMinusStrand = FALSE)))
  stopifnot(identical(f_factor_minus, Rsamtools::scanBamFlag(isMinusStrand = TRUE)))
  stopifnot(identical(f_factor_rev_minus, Rsamtools::scanBamFlag(isMinusStrand = FALSE)))
}

flip_gtf_strand <- function(gtf_in, gtf_out) {
  lines <- readLines(gtf_in)
  flipped <- vapply(lines, function(line) {
    if (startsWith(line, "#") || !nzchar(line)) return(line)
    fields <- strsplit(line, "\t", fixed = TRUE)[[1]]
    if (length(fields) >= 7) {
      if (fields[7] == "+") {
        fields[7] <- "-"
      } else if (fields[7] == "-") {
        fields[7] <- "+"
      }
    }
    paste(fields, collapse = "\t")
  }, character(1), USE.NAMES = FALSE)
  writeLines(flipped, gtf_out)
}

test_extdata_count_find_split <- function() {
  extdata_path <- system.file("extdata", package = "Sierra")
  gtf_file <- file.path(extdata_path, "Vignette_cellranger_genes_subset.gtf")
  bam_file <- file.path(extdata_path, "Vignette_example_TIP_sham.bam")
  junc_file <- file.path(extdata_path, "Vignette_example_TIP_sham_junctions.bed")
  wl_file <- file.path(extdata_path, "example_TIP_sham_whitelist_barcodes.tsv")
  peaks_file <- file.path(extdata_path, "TIP_merged_peaks.txt")

  work_dir <- tempfile("sierra_strand_test_")
  dir.create(work_dir, recursive = TRUE)
  on.exit(unlink(work_dir, recursive = TRUE), add = TRUE)

  # Subset peaks to a + strand gene (Ly6e on chr15) and a - strand gene (Lamp2 on chrX)
  peaks_df <- read.table(peaks_file, header = TRUE, sep = "\t", quote = "", stringsAsFactors = FALSE)
  sub_peaks <- subset(peaks_df, Gene %in% c("Ly6e", "Lamp2"))
  stopifnot(any(sub_peaks$Strand == 1), any(sub_peaks$Strand == -1))

  peaks_sub_path <- file.path(work_dir, "peaks_sub.txt")
  write.table(sub_peaks, peaks_sub_path, sep = "\t", quote = FALSE, row.names = FALSE)

  peaks_flipped <- sub_peaks
  peaks_flipped$Strand <- -peaks_flipped$Strand
  peaks_flipped_path <- file.path(work_dir, "peaks_flipped.txt")
  write.table(peaks_flipped, peaks_flipped_path, sep = "\t", quote = FALSE, row.names = FALSE)

  gtf_flipped_path <- file.path(work_dir, "genes_flipped.gtf")
  flip_gtf_strand(gtf_file, gtf_flipped_path)

  # CountPeaks: strandMode=1 on original == strandMode=2 on strand-flipped peaks/GTF
  dir_fwd <- file.path(work_dir, "count_fwd")
  dir_rev <- file.path(work_dir, "count_rev_flipped")
  dir_rev_orig <- file.path(work_dir, "count_rev_orig")
  dir_unstranded <- file.path(work_dir, "count_unstranded")

  CountPeaks(peaks_sub_path, gtf_file, bam_file, wl_file, dir_fwd, ncores = 1, strandMode = 1L)
  CountPeaks(peaks_flipped_path, gtf_flipped_path, bam_file, wl_file, dir_rev, ncores = 1, strandMode = 2L)
  CountPeaks(peaks_sub_path, gtf_file, bam_file, wl_file, dir_rev_orig, ncores = 1, strandMode = 2L)
  CountPeaks(peaks_sub_path, gtf_file, bam_file, wl_file, dir_unstranded, ncores = 1, strandMode = 0L)

  read_peak_totals <- function(d) {
    mtx_path <- if (file.exists(file.path(d, "matrix.mtx.gz"))) file.path(d, "matrix.mtx.gz") else file.path(d, "matrix.mtx")
    sites_path <- if (file.exists(file.path(d, "sitenames.tsv.gz"))) file.path(d, "sitenames.tsv.gz") else file.path(d, "sitenames.tsv")
    mtx <- if (grepl("\\.gz$", mtx_path)) Matrix::readMM(gzfile(mtx_path)) else Matrix::readMM(mtx_path)
    sites <- readLines(sites_path)
    if (length(sites) == 0L || nrow(mtx) == 0L) {
      return(numeric(0))
    }
    setNames(Matrix::rowSums(mtx), sub(":-?1$", "", sites))
  }

  tot_fwd <- read_peak_totals(dir_fwd)
  tot_rev_flipped <- read_peak_totals(dir_rev)
  tot_rev_orig <- read_peak_totals(dir_rev_orig)
  tot_unstranded <- read_peak_totals(dir_unstranded)

  stopifnot(sum(tot_fwd) > 0)
  stopifnot(identical(sort(names(tot_fwd)), sort(names(tot_rev_flipped))))
  stopifnot(isTRUE(all.equal(tot_fwd[names(tot_fwd)], tot_rev_flipped[names(tot_fwd)])))
  stopifnot(sum(tot_unstranded) >= sum(tot_fwd))
  stopifnot(sum(tot_rev_orig) < sum(tot_fwd))

  # FindPeaks: strandMode=1 on original finds peaks; strandMode=2 on flipped GTF finds matching coordinates
  out_find_fwd <- file.path(work_dir, "find_fwd.txt")
  out_find_rev_flip <- file.path(work_dir, "find_rev_flip.txt")
  FindPeaks(
    output.file = out_find_fwd,
    gtf.file = gtf_file,
    bamfile = bam_file,
    junctions.file = junc_file,
    chr.names = c("15", "X"),
    ncores = 1,
    strandMode = 1L
  )
  FindPeaks(
    output.file = out_find_rev_flip,
    gtf.file = gtf_flipped_path,
    bamfile = bam_file,
    junctions.file = junc_file,
    chr.names = c("15", "X"),
    ncores = 1,
    strandMode = 2L
  )
  df_find_fwd <- read.table(out_find_fwd, header = TRUE, sep = "\t", stringsAsFactors = FALSE)
  df_find_rev <- read.table(out_find_rev_flip, header = TRUE, sep = "\t", stringsAsFactors = FALSE)
  stopifnot(nrow(df_find_fwd) > 0)
  stopifnot(nrow(df_find_fwd) == nrow(df_find_rev))
  stopifnot(identical(df_find_fwd$Fit.max.pos, df_find_rev$Fit.max.pos))

  # SplitBam: strandMode=2 on flipped GTF matches strandMode=1 on original GTF
  gtf_gr <- rtracklayer::import(gtf_file)
  gtf_flip_gr <- rtracklayer::import(gtf_flipped_path)
  bcs <- readLines(wl_file)
  cellbc_df <- data.frame(
    celltype = rep("TestType", min(50L, length(bcs))),
    cellbc = head(bcs, 50L),
    stringsAsFactors = FALSE
  )
  rle_fwd <- SplitBam(bam_file, cellbc_df, outdir = NULL, gtf_gr = gtf_gr, geneSymbol = "Lamp2", rle_output = TRUE, strandMode = 1L)
  rle_rev <- SplitBam(bam_file, cellbc_df, outdir = NULL, gtf_gr = gtf_flip_gr, geneSymbol = "Lamp2", rle_output = TRUE, strandMode = 2L)
  stopifnot(isTRUE(all.equal(rle_fwd, rle_rev)))
}

tmp_lib_check()
test_alignment_helpers()
test_extdata_count_find_split()
cat("ALL_SIERRA_STRANDMODE_TESTS_PASSED\n")
