


#######################################################################################
#'
#' Generates a table of similarity measures between two sets of peaks
#'
#' Goes through the set of genes contained in peaks.1. For each gene-specific peak,
#' calculate the amount of overlapping nucleotides to the nearest peak in peaks.2. If the gene
#' is not available in peaks.2, ditance is set to -1e7.
#'
#' @param peaks.1 first set of peaks - used as a reference point
#' @param peaks.2 second set of peaks being compared
#' @param plot.distribution whether to plot the distribution of similarites between the peaks (default: FALSE).
#' @return a data-frame with peaks from peaks.1 mapped to the closest corresponding peak in peaks.2.
#' @examples
#' generate_similarity_table(peaks.1, peaks.2)
#'
#' @importFrom magrittr "%>%"
#'
generate_similarity_table <- function(peaks.1, peaks.2, ncores = 1) {

  ## Pull out gene names
  gene.names1 = sub("(*):.*", "\\1", peaks.1)
  gene.names2 = sub("(*):.*", "\\1", peaks.2)

  gene.names1.unique = unique(gene.names1)
  gene.names1.unique = gene.names1.unique[which(gene.names1.unique != "")]
  complete.apa.table = c()

  
  # Set up multiple workers
  system.name <- Sys.info()['sysname']
  new_cl <- FALSE
  if (system.name == "Windows") {
    new_cl <- TRUE
    cluster <- parallel::makePSOCKcluster(rep("localhost", ncores))
    doParallel::registerDoParallel(cluster)
  } else {
    doParallel::registerDoParallel(cores=ncores)
  }
  
  ## Go through the peaks by pulling out peaks associates with each gene
  complete.apa.table <- foreach::foreach(gene.name = gene.names1.unique, .combine = 'rbind') %dopar% {

    ## Get positions for first data-set
    gene.apa1 = peaks.1[which(gene.names1 == gene.name)]
    start.pos = as.numeric(sub(".*:(.*)-.*:.*", "\\1", gene.apa1))
    end.pos = as.numeric(sub(".*:.*-(.*):.*", "\\1", gene.apa1))
    peaks1.table = data.frame(Peak = gene.apa1, Start = start.pos, End = end.pos)
    rownames(peaks1.table) = gene.apa1

    ## Get positions for second data-set
    gene.apa2 = peaks.2[which(gene.names2 == gene.name)]
    start.pos = as.numeric(sub(".*:(.*)-.*:.*", "\\1", gene.apa2))
    end.pos = as.numeric(sub(".*:.*-(.*):.*", "\\1", gene.apa2))
    peaks2.table = data.frame(Peak = gene.apa2, Start = start.pos, End = end.pos)
    rownames(peaks2.table) = gene.apa2

    ## Check if the gene is actually in the next data-set
    if (nrow(peaks2.table) > 0) {
      gene.apa.table = c()
      for (this.apa in seq(1, nrow(peaks1.table))) {
        this.start = peaks1.table[this.apa, 2]
        this.end = peaks1.table[this.apa, 3]

        ## Find the closest peak
        start.difs = this.start - peaks2.table[, 2]
        end.difs = this.end - peaks2.table[, 3]
        total.difs = abs(start.difs) + abs(end.difs)
        closest.peak.index = which(total.difs == min(total.difs))
        similarity = 1 - total.difs[closest.peak.index]/((this.end + 1) - this.start)
        closest.peak = rownames(peaks2.table)[closest.peak.index]

        ## Check if there are more than one matches
        if (length(similarity) > 1) {
          if (similarity[1] > 0.1) {
            print(paste0("Tied peak similarities with peak ", rownames(peaks1.table)[this.apa],
                         " with proportion = ", similarity[1]))
          }
          similarity = similarity[1]
          closest.peak = closest.peak[1]
          closest.peak.index = closest.peak.index[1]
        }

        ## If closest peak far enough away similarity will be negative
        ## In that case, set it to 0
        if (similarity < 0) {
          similarity = 0
        }

        ## Calculate similarity from the other direction
        dataset2.start = peaks2.table[closest.peak.index, 2]
        dataset2.end = peaks2.table[closest.peak.index, 3]

        start.dif = dataset2.start - this.start
        end.dif = dataset2.end - this.end
        total.dif = abs(start.dif) + abs(end.dif)
        dataset2.similarity = 1 - total.dif/((dataset2.end + 1) - dataset2.start)
        if (dataset2.similarity < 0) {
          dataset2.similarity = 0
        }

        ## Add information to table
        peaks1.table[this.apa, ] %>%
          dplyr::mutate(Data2_Closest_Peak = closest.peak, Data2_Start_Dif = start.difs[closest.peak.index],
                 Data2_End_Dif = end.difs[closest.peak.index], Data2_Similarity = similarity,
                 Data2_Data1_Similarity = dataset2.similarity) %>% as.data.frame() ->
          newLine
        gene.apa.table = rbind(gene.apa.table, newLine)
      }

      #complete.apa.table = rbind(complete.apa.table, gene.apa.table)
      return(gene.apa.table)

    } else {
      ## Gene isn't in the comparison data-set - set closest peak to -1e7
      peaks1.table %>%
        dplyr::mutate(Data2_Closest_Peak = -1e7, Data2_Start_Dif = -1e7,
               Data2_End_Dif = -1e7, Data2_Similarity = 0, Data2_Data1_Similarity = 0) %>% as.data.frame() ->
        newLine

      return(newLine)
    }
  }
  
  if (new_cl) { ## Shut down cluster if on Windows
    ## stop cluster
    parallel::stopCluster(cluster)
  }

  rownames(complete.apa.table) = complete.apa.table$Peak
  return(complete.apa.table)
}



####################################################################################
#'
#' Generates a table of similarity measures within a set of peaks
#'
#' In some rare cases, called peaks will show a high degree of overlap, and before
#' merging two different sets of peaks, the similar peaks within a set first need
#' to be merged. This function looks for the most similar peak (non-self) within a set of peaks
#' and calculates the level of overlap.
#'
#' @param peaks.1 the set of peaks to merge
#' @param plot.distribution whether to plot the distribution of similarites between the peaks (default: FALSE).
#' @return a data-frame with peaks from peaks.1 mapped to the closest peak within itself
#' @examples
#' generate_similarity_table(peaks.1)
#'
#' @importFrom magrittr "%>%"
#' @importFrom foreach "%dopar%"
#'
generate_self_similarity_table <- function(peaks.1, ncores = 1) {

  peaks.2 = peaks.1

  ## Pull out gene names
  gene.names1 = sub("(*):.*", "\\1", peaks.1)
  gene.names2 = sub("(*):.*", "\\1", peaks.2)

  ## remove any empty string gene names
  peaks.1 = peaks.1[which(gene.names1 != "")]
  gene.names1 = gene.names1[which(gene.names1 != "")]

  gene.names1.unique = unique(gene.names1)

  # Set up multiple workers
  system.name <- Sys.info()['sysname']
  new_cl <- FALSE
  if (system.name == "Windows") {
    new_cl <- TRUE
    cluster <- parallel::makePSOCKcluster(rep("localhost", ncores))
    doParallel::registerDoParallel(cluster)
  } else {
    doParallel::registerDoParallel(cores=ncores)
  }
  
  ## Go through the peaks by pulling out peaks associates with each gene
  complete.apa.table <- foreach::foreach(gene.name = gene.names1.unique, .combine = 'rbind') %dopar% {
    ## Get positions for first data-set
    gene.apa1 = peaks.1[which(gene.names1 == gene.name)]
    start.pos = as.numeric(sub(".*:.*:(.*)-.*:.*", "\\1", gene.apa1))
    end.pos = as.numeric(sub(".*:.*:.*-(.*):.*", "\\1", gene.apa1))
    peaks1.table = data.frame(Peak = gene.apa1, Start = start.pos, End = end.pos, stringsAsFactors = FALSE)
    rownames(peaks1.table) = gene.apa1

    ## Get positions for second data-set
    gene.apa2 = peaks.2[which(gene.names2 == gene.name)]
    start.pos = as.numeric(sub(".*:.*:(.*)-.*:.*", "\\1", gene.apa2))
    end.pos = as.numeric(sub(".*:.*:.*-(.*):.*", "\\1", gene.apa2))
    peaks2.table = data.frame(Peak = gene.apa2, Start = start.pos, End = end.pos, stringsAsFactors = FALSE)
    rownames(peaks2.table) = gene.apa2

    ## Check if the gene is actually in the next data-set
    if (nrow(peaks2.table) > 1) {
      gene.apa.table = data.frame(matrix(NA, nrow = nrow(peaks1.table), ncol = 8))
      colnames(gene.apa.table) = c("Peak", "Start", "End", "Data2_Closest_Peak",
                                   "Data2_Start_Dif", "Data2_End_Dif", "Data2_Similarity", "Data2_Data1_Similarity")
      for (this.apa in seq(1, nrow(peaks1.table))) {
        this.start = peaks1.table[this.apa, 2]
        this.end = peaks1.table[this.apa, 3]

        ## Remove the test peak from peaks2.table
        peaks2.table = peaks1.table[-this.apa, ]

        ## Find the closest peak
        start.difs = this.start - peaks2.table[, 2]
        end.difs = this.end - peaks2.table[, 3]
        total.difs = abs(start.difs) + abs(end.difs)
        closest.peak.index = which(total.difs == min(total.difs))
        similarity = 1 - total.difs[closest.peak.index]/((this.end + 1) - this.start)
        closest.peak = rownames(peaks2.table)[closest.peak.index]

        ## Check if there are more than one matches
        if (length(similarity) > 1) {
          if (similarity[1] > 0.1) {
            print(paste0("Tied peak similarities with peak ", rownames(peaks1.table)[this.apa],
                         " with proportion = ", similarity[1]))
          }
          similarity = similarity[1]
          closest.peak = closest.peak[1]
          closest.peak.index = closest.peak.index[1]
        }

        ## If closest peak far enough away similarity will be negative
        ## In that case, set it to 0
        if (similarity < 0) {
          similarity = 0
        }

        ## Calculate similarity from the other direction
        dataset2.start = peaks2.table[closest.peak.index, 2]
        dataset2.end = peaks2.table[closest.peak.index, 3]

        start.dif = dataset2.start - this.start
        end.dif = dataset2.end - this.end
        total.dif = abs(start.dif) + abs(end.dif)
        dataset2.similarity = 1 - total.dif/((dataset2.end + 1) - dataset2.start)
        if (dataset2.similarity < 0) {
          dataset2.similarity = 0
        }

        ## Add information to table
        peaks1.table[this.apa, ] %>%
          dplyr::mutate(Data2_Closest_Peak = closest.peak, Data2_Start_Dif = start.difs[closest.peak.index],
                 Data2_End_Dif = end.difs[closest.peak.index], Data2_Similarity = similarity,
                 Data2_Data1_Similarity = dataset2.similarity) %>% as.data.frame() ->
          newLine
        gene.apa.table[this.apa, ] = newLine
      }

      ## Update the complete table
      return(gene.apa.table)

    } else {
      ## Gene isn't in the comparison data-set - set closest peak to -1e7
      peaks1.table %>%
        dplyr::mutate(Data2_Closest_Peak = -1e7, Data2_Start_Dif = -1e7,
               Data2_End_Dif = -1e7, Data2_Similarity = 0, Data2_Data1_Similarity = 0) %>% as.data.frame() ->
        newLine

      return(newLine)
    }

  }
  
  if (new_cl) { ## Shut down cluster if on Windows
    ## stop cluster
    parallel::stopCluster(cluster)
  }

  ## Add column names to the output table
  colnames(complete.apa.table) = c("Peak", "Start", "End", "Data2_Closest_Peak",
                                   "Data2_Start_Dif", "Data2_End_Dif", "Data2_Similarity", "Data2_Data1_Similarity")
  return(complete.apa.table)
}


#############################################################
#'
#' Merge a set of peaks
#'
#' Given a self-similarity table of peaks, identify peaks that should be merged. Merged peaks are
#' taken as the union of the two peaks. For two given peaks, A and B, they will be merged if at least one has
#' some x\% (75\% by default) or more overlap with the other, and the other has at least x-(y*x)\% overlap where
#' y is a percentage of allowed variance (25\% by default)
#'
#' @param apa.similarity.table the set of peaks to merge
#' @param sim.thresh The required similarity threshold for merging (default: 0.75)
#' @param allow.match.var The allowance for deviation from the sim.thresh for comparison peaks (default: 0.25)
#' @param return.peaks Whether to return a full table of results or simply a vector of merged peaks
#' @return a table of merged peaks with original merged peaks or a vector of merged peaks
#' @examples
#' generate_self_merged_peaks(apa.similarity.table)
#'
#' @importFrom magrittr "%>%"
#'
generate_self_merged_peaks <- function(apa.similarity.table, sim.thresh = 0.75, allow.match.var = 0.25,
                                       return.type = "peak_list") {

  apa.similarity.table %>%
    dplyr::mutate(Matched_Similarity =
             ((Data2_Similarity > sim.thresh & Data2_Data1_Similarity > (sim.thresh - allow.match.var * sim.thresh))) |
             (Data2_Data1_Similarity > sim.thresh & Data2_Similarity > (sim.thresh - allow.match.var * sim.thresh))) ->
    apa.similarity.table
  apa.similarity.table %>% dplyr::mutate(GeneName = sub("(*):.*", "\\1", apa.similarity.table$Peak)) -> apa.similarity.table
  rownames(apa.similarity.table) = apa.similarity.table$Peak

  table1.merge1.subset = subset(apa.similarity.table, Matched_Similarity == TRUE)
  rownames(table1.merge1.subset) = table1.merge1.subset$Peak
  peaks.merge1 = rownames(table1.merge1.subset)

  ### Create a mapping table of peaks to merge and new peak name
  data2.start.pos = as.numeric(sub(".*:.*:(.*)-.*:.*", "\\1", table1.merge1.subset$Data2_Closest_Peak))
  data2.end.pos = as.numeric(sub(".*:.*:.*-(.*):.*", "\\1", table1.merge1.subset$Data2_Closest_Peak))
  strand = sub(".*:.*:.*-.*:(.*)", "\\1", table1.merge1.subset$Peak)
  chr = sub(".*:(.*):.*-.*:.*", "\\1", table1.merge1.subset$Peak)

  peak.mapping.table = data.frame(Gene = table1.merge1.subset$GeneName, Chromosome = chr, Strand = strand,
                                  Data1_Start = table1.merge1.subset$Start, Data2_Start = data2.start.pos,
                                  Data1_End = table1.merge1.subset$End, Data2_End = data2.end.pos)

  ## Calculate the minimum start point and maximum end point to create the merged peak
  min.start = apply(peak.mapping.table, 1, function(x) min(as.numeric(x["Data1_Start"]), as.numeric(x["Data2_Start"])))
  peak.mapping.table %>% dplyr::mutate(Start_Min = min.start) -> peak.mapping.table
  max.end = apply(peak.mapping.table, 1, function(x) max(as.numeric(x["Data1_End"]), as.numeric(x["Data2_End"])))
  peak.mapping.table %>% dplyr::mutate(End_Max = max.end) -> peak.mapping.table

  ## Add updated peak name to the table
  peak.mapping.table %>% dplyr::mutate(MergedPeak = sprintf("%s:%s:%d-%d:%s", Gene, Chromosome, Start_Min, End_Max, Strand),
                                Data1_Peak = table1.merge1.subset$Peak, Data2_Peak = table1.merge1.subset$Data2_Closest_Peak) ->
    peak.mapping.table

  ## Remove duplicated merged peaks
  peak.mapping.table %>% dplyr::distinct(MergedPeak, .keep_all = TRUE) -> peak.mapping.table.unique

  ## Add all the non-overlapping peaks
  table1.distinct.subset = subset(apa.similarity.table, Matched_Similarity == FALSE)
  distinct.peaks = table1.distinct.subset$Peak

  ## Combined the merged and unique peaks
  combined.peaks = c(peak.mapping.table.unique$MergedPeak, table1.distinct.subset$Peak)
  if (return.type == "peak_list") {
        return(combined.peaks)
    } else if (return.type == "merged_table") {
        return(peak.mapping.table) ## return a table of peaks that should be merged
    } else if (return.type == "both") {
        res.output <- list(combined.peaks, peak.mapping.table)
        names(res.output) <- c("peak_list", "merged_table")
        return(res.output)
    } else {
    warning("invalid output option provided")
  }
}


#############################################################################
#'
#' Merge peaks across data-sets based on a reference
#'
#' Given a reference data-set, a list of data-sets for merging and set of merged peaks from the referece,
#' identify peaks that should be merged. Merged peaks are taken as the union of the peaks to be merged.
#' For two given peaks, A and B, they will be merged if at least one has some x\% (75\% by default) or more
#' overlap with the other, and the other has at least x-(y*x)\% overlap where y is a percentage of allowed
#' variance (25\% by default).
#'
#' @param dataset.1 the reference peak data-set
#' @param peak.dataset.list a list of peak data-sets
#' @param self.merged.peaks.list the set of self-merged peaks from the reference data-set
#' @param sim.thresh The required similarity threshold for merging (default: 0.75)
#' @param allow.match.var The allowance for deviation from the sim.thresh for comparison peaks (default: 0.25)
#' @return a data-frame containing peaks, their class (merged or unique) and the original peak from the reference
#' @examples
#' generate_merged_peak_table(dataset.1, peak.dataset.table, self.merged.peaks.list)
#'
#' @importFrom magrittr "%>%"
#'
generate_merged_peak_table <- function(dataset.1, peak.dataset.list, self.merged.peaks.list,
                                       sim.thresh = 0.75, allow.match.var = 0.25, ncores = 1) {
  peak.set1 = self.merged.peaks.list[[dataset.1]]

  gene.names1 = sub("(*):.*", "\\1", peak.set1)

  ## remove any empty string gene names
  combined.peak.similarity.table = data.frame(Peak = peak.set1, stringsAsFactors = FALSE)
  rownames(combined.peak.similarity.table) = combined.peak.similarity.table$Peak

  similarity.labels = c()
  rev.similarity.labels = c()
  closest.peak.labels = c()
  match.labels = c()
  for (i in 1:length(peak.dataset.list)) {
    comparison.dataset = names(peak.dataset.list)[i]
    if (dataset.1 != comparison.dataset) {
      comparison.peak.set = self.merged.peaks.list[[comparison.dataset]]
      this.peak.similarity.table = generate_similarity_table(peak.set1, comparison.peak.set, ncores = ncores)
      table.subset = this.peak.similarity.table[combined.peak.similarity.table$Peak,
                                                c("Data2_Closest_Peak", "Data2_Similarity", "Data2_Data1_Similarity")]
      table.subset %>%
        dplyr::mutate(Matched_Similarity =
                 ((Data2_Similarity > sim.thresh & Data2_Data1_Similarity > (sim.thresh - allow.match.var * sim.thresh))) |
                 (Data2_Data1_Similarity > sim.thresh & Data2_Similarity > (sim.thresh - allow.match.var * sim.thresh))) ->
        table.subset
      similarity.label = paste0(comparison.dataset, "_Similarity")
      similarity.labels = append(similarity.labels, similarity.label)
      rev.similarity.label = paste0(comparison.dataset, "_", dataset.1, "_Similarity")
      rev.similarity.labels = append(rev.similarity.labels, rev.similarity.label)
      closest.peak.label = paste0(comparison.dataset, "_Closest_Peak")
      closest.peak.labels = append(closest.peak.labels, closest.peak.label)
      match.label = paste0(comparison.dataset, "_Matched_Similarity")
      match.labels = append(match.labels, match.label)
      colnames(table.subset) = c(closest.peak.label, similarity.label, rev.similarity.label, match.label)
      combined.peak.similarity.table = cbind(combined.peak.similarity.table, table.subset)
    }
  }

  ## Identify the peaks that will be merged
  if (length(match.labels) > 1) {
    match.counts = apply(combined.peak.similarity.table[, match.labels], 1, function(x) sum(x==TRUE))
  } else {
    match.counts = unlist(lapply(combined.peak.similarity.table[, match.labels], function(x) sum(x==TRUE)))
    names(match.counts) = rownames(combined.peak.similarity.table)
  }

  peaks.to.merge = names(match.counts[which(match.counts > 0)])

  combined.peak.similarity.table.subset = combined.peak.similarity.table[peaks.to.merge, ]

  ### --- Create a mapping table of peaks to merge and new peak name for first table --- ###

  ## begin by identifying the miniminum start positions for the matched peaks
  start.positions = t(apply(combined.peak.similarity.table.subset[, c("Peak", closest.peak.labels)], 1, function(x)
    as.numeric(sub(".*:.*:(.*)-.*:.*", "\\1", x))))
  colnames(start.positions) = c("Peak", closest.peak.labels)

  ## Go through all closest start sites and where matched similarities are true, take the minimum
  min.start.positions = vector("numeric", nrow(start.positions))
  for (i in 1:nrow(start.positions)) {
    start.set = start.positions[i, ][c(TRUE, as.logical(combined.peak.similarity.table.subset[i, match.labels]))]
    min.start.positions[i] = min(start.set)
  }

  ## Now calculate the maximum end positions for matched peaks
  end.positions = t(apply(combined.peak.similarity.table.subset[, c("Peak", closest.peak.labels)], 1, function(x)
    as.numeric(sub(".*:.*:.*-(.*):.*", "\\1", x))))
  colnames(end.positions) = c("Peak", closest.peak.labels)

  ## Go through all closest start sites and where matched similarities are true, take the minimum
  max.end.positions = vector("numeric", nrow(end.positions))
  for (i in 1:nrow(end.positions)) {
    end.set = end.positions[i, ][c(TRUE, as.logical(combined.peak.similarity.table.subset[i, match.labels]))]
    max.end.positions[i] = max(end.set)
  }

  ## Now build a peak mapping table
  strand = sub(".*:.*:.*-.*:(.*)", "\\1", combined.peak.similarity.table.subset$Peak)
  gene = sub("(.*):.*:.*-.*:.*", "\\1", combined.peak.similarity.table.subset$Peak)
  chr = sub(".*:(.*):.*-.*:.*", "\\1", combined.peak.similarity.table.subset$Peak)

  peak.mapping.table = data.frame(OriginalPeak = combined.peak.similarity.table.subset$Peak, Gene = gene,
                                  Chromosome = chr, Strand = strand, Start_Min = min.start.positions,
                                  End_Max = max.end.positions, stringsAsFactors = FALSE)


  ## Add updated peak name to the table
  peak.mapping.table %>% dplyr::mutate(MergedPeak = sprintf("%s:%s:%d-%d:%s", Gene, Chromosome, Start_Min, End_Max, Strand)) -> peak.mapping.table

  ## Remove duplicated merged peaks
  peak.mapping.table %>% dplyr::distinct(MergedPeak, .keep_all = TRUE) -> peak.mapping.table.unique

  ## Now pull out the unique peaks
  peaks.unique = names(match.counts[which(match.counts == 0)])

  combined.peak.output.table = data.frame(Peak = c(peak.mapping.table.unique$MergedPeak, peaks.unique),
                                          Class = c(rep("Merged", nrow(peak.mapping.table.unique)),
                                                    rep(paste0("Unique_", dataset.1), length(peaks.unique))),
                                          OriginalPeak = c(peak.mapping.table.unique$OriginalPeak, peaks.unique),
                                          DataOrigin = rep(dataset.1, times = (nrow(peak.mapping.table.unique) + length(peaks.unique))),
                                          stringsAsFactors = FALSE)
  return(combined.peak.output.table)
}



####################################################################################
#'
#' Merge peaks across a list of data-sets
#'
#' Takes as input a list of named peaks obtained from running peak calling on multiple data-sets.
#' First goes through each peak set and check what peaks within each set should be merged (self-merging).
#' Merging is based on similarity criteria set by sim.thresh and allow.match.var.
#' Then compares each peak set as a reference to the remaining sets to identify peaks that should be merged.
#' Returns a list of peaks that have been merged, as well as the unique peaks from each data-set.
#'
#' @param peak.dataset.table a dataframe with two required columnss: one called "Peak_file" , which
#' contains file names of the peak data-sets to be merged and labels ("Identifier") for each file.
#' @param output.file file to write the set of merged peaks to
#' @param sim.thresh The required similarity threshold for merging (default: 0.75)
#' @param allow.match.var The allowance for deviation from the sim.thresh for comparison peaks (default: 0.25)
#' @param ncores number of cores to use (default 1)
#' @return NULL. writes out a set of merged peaks to output.file
#' @examples
#' MergePeakCoordinates(peak.dataset.table, output.file, ncores = 4)
#'
#' @importFrom magrittr "%>%"
#'
#' @export
#'
MergePeakCoordinates <- function(peak.dataset.table, output.file, sim.thresh = 0.75,
                            allow.match.var = 0.25, ncores = 1) {

  ## Create a named list from the peaks
  ## While reading in the peak files retain information on junctions - will use this later
  peak.junctions.list = c()
  peak.dataset.list = c()
  for (i in 1:nrow(peak.dataset.table)) {
    peak.table <- read.table(as.character(peak.dataset.table[i, "Peak_file"]),
                             header = TRUE, stringsAsFactors = FALSE)
    ## list of unique peak IDs
    this.peak.list = unique(peak.table$polyA_ID)
    peak.dataset.list = c(peak.dataset.list, list(this.peak.list))

    ## table mapping peak ID to junction status
    junction.table = peak.table[, c("polyA_ID", "exon.intron")]
    junction.table %>% dplyr::distinct(polyA_ID, .keep_all = TRUE) -> junction.table
    rownames(junction.table) <- junction.table$polyA_ID
    peak.junctions.list = c(peak.junctions.list, list(junction.table))
  }
  names(peak.dataset.list) = as.character(peak.dataset.table$Identifier)
  names(peak.junctions.list) = as.character(peak.dataset.table$Identifier)

  ### Generate the in-dataset similarity table
  self.merged.peaks.list = vector("list", length(peak.dataset.list))

  ## Go through the list and first generate self similarity tables to merge
  ## Return a list of peaks where similar peaks within each data-set have been merged
  for (i in 1:length(peak.dataset.list)) {
    dataset.name = names(peak.dataset.list)[i]
    this.peak.set = peak.dataset.list[[dataset.name]]
    print(paste0("Performing internal peak merging for ", dataset.name))
    this.peak.similarity.table = generate_self_similarity_table(this.peak.set, ncores = ncores)
    merging.res = generate_self_merged_peaks(this.peak.similarity.table, sim.thresh = sim.thresh,
                                             allow.match.var = allow.match.var, return.type = "both")
    ## List of all merged and original peaks
    peak.set.merged <- merging.res[["peak_list"]]
    self.merged.peaks.list[[i]] = peak.set.merged

    ## table just of the peaks that were merged
    merged.table  <- merging.res[["merged_table"]]

    ## use the table to add junction status for the merged peaks
    this.junction.table <- peak.junctions.list[[dataset.name]]
    junction.subset <- this.junction.table[merged.table$Data1_Peak, ]
    junction.subset$polyA_ID <- merged.table$MergedPeak
    junction.subset %>% dplyr::distinct(polyA_ID, .keep_all = TRUE) -> junction.subset
    rownames(junction.subset) <- junction.subset$polyA_ID
    this.junction.table <- rbind(this.junction.table, junction.subset)
    peak.junctions.list[[dataset.name]] <- this.junction.table
  }
  names(self.merged.peaks.list) = names(peak.dataset.list)

  ### Build a table of similarity between the reference data-set and remainder data-sets
  combined.merged.peaks.list = vector("list", length(peak.dataset.list))
  for (i in 1:length(peak.dataset.list)) {
    ref.dataset.name = names(peak.dataset.list)[i]
    print(paste0("Comparing peaks from ", ref.dataset.name, " to remaining data-sets"))
    this.merged.peak.table = generate_merged_peak_table(dataset.1 = ref.dataset.name, peak.dataset.list = peak.dataset.list,
                                                        self.merged.peaks.list = self.merged.peaks.list,
                                                        sim.thresh = sim.thresh, allow.match.var = allow.match.var,
                                                        ncores = ncores)
    combined.merged.peaks.list[[i]] = this.merged.peak.table
  }
  names(combined.merged.peaks.list) = names(peak.dataset.list)

  ### Identif the peaks that are among all the merged and the ones that are unique to each data-set
  ## start with the merged peaks
  all.peaks = do.call(rbind, combined.merged.peaks.list)

  merged.peaks = subset(all.peaks, Class == "Merged")
  unique.peaks = subset(all.peaks, Class != "Merged")
  unique.peaks %>% dplyr::distinct(Peak, .keep_all = TRUE) -> unique.peaks

  ## Remove duplicated merged peaks
  merged.peaks %>% dplyr::distinct(Peak, .keep_all = TRUE) -> merged.peaks.unique

  ## Add the set of unique peaks that aren't already in the merged list
  unique.peaks = unique.peaks[which(!(unique.peaks$Peak %in% merged.peaks.unique$Peak)), ]

  all.merged.peaks = rbind(merged.peaks.unique, unique.peaks)

  ### Run a final check on the peaks to clean up any that should be further merged
  final.peak.similarity.table = generate_self_similarity_table(all.merged.peaks$Peak, ncores = ncores)
  peaks.to.merge = generate_self_merged_peaks(final.peak.similarity.table, sim.thresh = sim.thresh,
                                              allow.match.var = allow.match.var, return.type = "merged_table")

  peaks.to.replace <- peaks.to.merge$Data1_Peak
  peaks.merge.index <- which(all.merged.peaks$Peak %in% peaks.to.replace)
  peaks.retain.index <- setdiff(1:nrow(all.merged.peaks), peaks.merge.index)
  all.merged.peaks.retained <- all.merged.peaks[peaks.retain.index, ] ## Keep these ones for later

  peaks.to.remerge <- all.merged.peaks[peaks.merge.index, ]
  final.peaks.combined <- all.merged.peaks.retained

  ## Iteratively check for peak similarity and merge until all peaks are considered unique
  while (nrow(peaks.to.merge) > 0) {
    ## First update the peak names to the merged peak names
    peaks.to.remerge$Peak <- plyr::mapvalues(peaks.to.remerge$Peak,
                                             from = peaks.to.merge$Data1_Peak,
                                             to = peaks.to.merge$MergedPeak)
    ## reduce to the set of unique merged peaks
    peaks.to.remerge %>% dplyr::distinct(Peak, .keep_all = TRUE) -> peaks.to.remerge

    ## Generate an updated similarity table and check for peaks to merge
    new.peak.similarity.table = generate_self_similarity_table(peaks.to.remerge$Peak, ncores = ncores)
    peaks.to.merge = generate_self_merged_peaks(new.peak.similarity.table, sim.thresh = sim.thresh,
                                                allow.match.var = allow.match.var, return.type = "merged_table")

    ## Separate peaks by those that can be retained now and those that
    ## will go through another round of mering
    peaks.to.replace <- peaks.to.merge$Data1_Peak
    peaks.merge.index <- which(peaks.to.remerge$Peak %in% peaks.to.replace)
    peaks.retain.index <- setdiff(1:nrow(peaks.to.remerge), peaks.merge.index)

    ## Append the retained peaks to the output table.
    merged.peaks.retained <- peaks.to.remerge[peaks.retain.index, ]
    final.peaks.combined <- rbind(final.peaks.combined, merged.peaks.retained)

    ## Peaks to remerge (if any)
    peaks.to.remerge <- peaks.to.remerge[peaks.merge.index, ]
  }

  final.merged.peaks <- final.peaks.combined

  final.merged.peaks %>% dplyr::mutate(Gene = sub("(.*):.*:.*-.*:.*", "\\1", Peak),
                                       Chr = sub(".*:(.*):.*-.*:.*", "\\1", Peak),
                                       Strand = sub(".*:.*:.*-.*:(.*)", "\\1", Peak),
                                       Fit.start = sub(".*:.*:(.*)-.*:.*", "\\1", Peak),
                                       Fit.end = sub(".*:.*:.*-(.*):.*", "\\1", Peak)) -> output.table
  output.table = output.table[, c("Gene", "Chr", "Strand", "Fit.start", "Fit.end", "Peak", "Class",
                                  "OriginalPeak", "DataOrigin")]
  colnames(output.table) = c("Gene", "Chr", "Strand", "Fit.start", "Fit.end",
                             "polyA_ID", "PeakClass", "OriginalPeak", "DataOrigin")

  ## Filter peaks

  ## Make sure the positions are numeric
  output.table$Fit.end <- as.numeric(output.table$Fit.end)
  output.table$Fit.start <- as.numeric(output.table$Fit.start)

  sites.diffs = output.table$Fit.end - output.table$Fit.start
  sites.keep = which(sites.diffs > 0)
  output.table = output.table[sites.keep, ]

  ## As a final step, map back status on junctions from the original input files
  data.classes <- unique(output.table$DataOrigin)


  junction.status <- rep("Unknown", nrow(output.table))
  output.table$'exon.intron' <- junction.status

  for (dataset.name in peak.dataset.table$Identifier) {
    dataset.index <- which(output.table$DataOrigin == dataset.name)
    original.peaks <- output.table[dataset.index, "OriginalPeak"]
    this.junction.table <- peak.junctions.list[[dataset.name]]
    original.peaks.overlap <- intersect(original.peaks, rownames(this.junction.table))
    this.junction.labels <- this.junction.table[original.peaks.overlap, "exon.intron"]

    ## add the junction labels to the output table
    peak.overlap.index <- which(output.table$OriginalPeak %in% original.peaks.overlap)
    output.table[peak.overlap.index, "exon.intron"] <- this.junction.labels
  }

  write.table(output.table, file = output.file, sep="\t", quote = FALSE, row.names = FALSE)
}
