


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
generate_similarity_table <- function(peaks.1, peaks.2, plot.distributions = FALSE) {

  ## Pull out gene names
  gene.names1 = sub("(*):.*", "\\1", peaks.1)
  gene.names2 = sub("(*):.*", "\\1", peaks.2)

  gene.names1.unique = unique(gene.names1)
  gene.names1.unique = gene.names1.unique[which(gene.names1.unique != "")]
  complete.apa.table = c()

  pb <- progress_bar$new(format = "Processing [:bar] :percent eta: :eta",
                         total = length(gene.names1.unique), clear=FALSE)
  pb$tick(0)
  for (gene.name in gene.names1.unique) {
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
          mutate(Data2_Closest_Peak = closest.peak, Data2_Start_Dif = start.difs[closest.peak.index],
                 Data2_End_Dif = end.difs[closest.peak.index], Data2_Similarity = similarity,
                 Data2_Data1_Similarity = dataset2.similarity) %>% as.data.frame() ->
          newLine
        gene.apa.table = rbind(gene.apa.table, newLine)
      }

      complete.apa.table = rbind(complete.apa.table, gene.apa.table)
    } else {
      ## Gene isn't in the comparison data-set - set closest peak to -1e7
      peaks1.table %>%
        mutate(Data2_Closest_Peak = -1e7, Data2_Start_Dif = -1e7,
               Data2_End_Dif = -1e7, Data2_Similarity = 0, Data2_Data1_Similarity = 0) %>% as.data.frame() ->
        newLine
      complete.apa.table = rbind(complete.apa.table, newLine)
    }
    pb$tick()
  }

  if (plot.distributions) {
    ## Plot the distribution of best similarities
    sub.apa.table = complete.apa.table[which(complete.apa.table$Data2_Similarity >= 0), ]

    ## plot distribution of peak similarity between data-sets
    ggData = data.frame(Similarity = complete.apa.table$Data2_Similarity, Data = rep("1", nrow(complete.apa.table)),
                        Start = complete.apa.table$Data2_Start_Dif, End = complete.apa.table$Data2_End_Dif)

    ## Plot Discordance differences using a box plot
    ggplot(ggData, aes(Similarity)) + geom_density() + theme_bw(base_size = 15) + ylab("Density") +
      ggtitle("T-regulatory to T-helper closest peak similarities")
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
generate_self_similarity_table <- function(peaks.1, plot.distributions = TRUE) {

  peaks.2 = peaks.1

  ## Pull out gene names
  gene.names1 = sub("(*):.*", "\\1", peaks.1)
  gene.names2 = sub("(*):.*", "\\1", peaks.2)

  ## remove any empty string gene names
  peaks.1 = peaks.1[which(gene.names1 != "")]
  gene.names1 = gene.names1[which(gene.names1 != "")]

  gene.names1.unique = unique(gene.names1)

  ## Initialise the results table. Number of final entries will equal APA with non empty gene names
  complete.apa.table = data.frame(matrix(NA, nrow = length(gene.names1), ncol = 8))
  colnames(complete.apa.table) = c("Peak", "Start", "End", "Data2_Closest_Peak",
                                   "Data2_Start_Dif", "Data2_End_Dif", "Data2_Similarity", "Data2_Data1_Similarity")
  current.start = 1
  current.end = 0

  pb <- progress_bar$new(format = "Processing [:bar] :percent eta: :eta",
                         total = length(gene.names1.unique), clear=FALSE)
  pb$tick(0)
  for (gene.name in gene.names1.unique) {
    ## Get positions for first data-set
    gene.apa1 = peaks.1[which(gene.names1 == gene.name)]
    start.pos = as.numeric(sub(".*:(.*)-.*:.*", "\\1", gene.apa1))
    end.pos = as.numeric(sub(".*:.*-(.*):.*", "\\1", gene.apa1))
    peaks1.table = data.frame(Peak = gene.apa1, Start = start.pos, End = end.pos, stringsAsFactors = FALSE)
    rownames(peaks1.table) = gene.apa1

    ## Get positions for second data-set
    gene.apa2 = peaks.2[which(gene.names2 == gene.name)]
    start.pos = as.numeric(sub(".*:(.*)-.*:.*", "\\1", gene.apa2))
    end.pos = as.numeric(sub(".*:.*-(.*):.*", "\\1", gene.apa2))
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
          mutate(Data2_Closest_Peak = closest.peak, Data2_Start_Dif = start.difs[closest.peak.index],
                 Data2_End_Dif = end.difs[closest.peak.index], Data2_Similarity = similarity,
                 Data2_Data1_Similarity = dataset2.similarity) %>% as.data.frame() ->
          newLine
        gene.apa.table[this.apa, ] = newLine
      }

      ## Update the complete table
      current.end = current.end + nrow(gene.apa.table)
      complete.apa.table[current.start:current.end, ] = gene.apa.table
      current.start = current.end + 1
      #complete.apa.table = rbind(complete.apa.table, gene.apa.table)
    } else {
      ## Gene isn't in the comparison data-set - set closest peak to -1e7
      peaks1.table %>%
        mutate(Data2_Closest_Peak = -1e7, Data2_Start_Dif = -1e7,
               Data2_End_Dif = -1e7, Data2_Similarity = 0, Data2_Data1_Similarity = 0) %>% as.data.frame() ->
        newLine
      current.end = current.end + 1
      complete.apa.table[current.start, ] = newLine
      current.start = current.start + 1
      #complete.apa.table = rbind(complete.apa.table, newLine)
    }
    pb$tick()
  }

  if (plot.distributions) {
    ## Plot the distribution of best similarities
    sub.apa.table = complete.apa.table[which(complete.apa.table$Data2_Similarity >= 0), ]

    ## plot distribution of peak similarity between data-sets
    ggData = data.frame(Similarity = complete.apa.table$Data2_Similarity, Data = rep("1", nrow(complete.apa.table)),
                        Start = complete.apa.table$Data2_Start_Dif, End = complete.apa.table$Data2_End_Dif)

    ## Plot Discordance differences using a box plot
    ggplot(ggData, aes(Similarity)) + geom_density() + theme_bw(base_size = 15) + ylab("Density") +
      ggtitle("T-regulatory to T-helper closest peak similarities")
  }

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
#' @return a vector of merged peaks
#' @examples
#' generate_similarity_table(peaks.1)
generate_self_merged_peaks <- function(apa.similarity.table, sim.thresh = 0.75, allow.match.var = 0.25) {

  apa.similarity.table %>%
    mutate(Matched_Similarity =
             ((Data2_Similarity > sim.thresh & Data2_Data1_Similarity > (sim.thresh - allow.match.var * sim.thresh))) |
             (Data2_Data1_Similarity > sim.thresh & Data2_Similarity > (sim.thresh - allow.match.var * sim.thresh))) ->
    apa.similarity.table
  apa.similarity.table %>% mutate(GeneName = sub("(*):.*", "\\1", apa.similarity.table$Peak)) -> apa.similarity.table
  rownames(apa.similarity.table) = apa.similarity.table$Peak

  table1.merge1.subset = subset(apa.similarity.table, Matched_Similarity == TRUE)
  rownames(table1.merge1.subset) = table1.merge1.subset$Peak
  peaks.merge1 = rownames(table1.merge1.subset)

  ### Create a mapping table of peaks to merge and new peak name
  data2.start.pos = as.numeric(sub(".*:(.*)-.*:.*", "\\1", table1.merge1.subset$Data2_Closest_Peak))
  data2.end.pos = as.numeric(sub(".*:.*-(.*):.*", "\\1", table1.merge1.subset$Data2_Closest_Peak))
  strand = sub(".*:.*-.*:(.*)", "\\1", table1.merge1.subset$Peak)

  peak.mapping.table = data.frame(Gene = table1.merge1.subset$GeneName, Strand = strand,
                                  Data1_Start = table1.merge1.subset$Start, Data2_Start = data2.start.pos,
                                  Data1_End = table1.merge1.subset$End, Data2_End = data2.end.pos)

  ## Calculate the minimum start point and maximum end point to create the merged peak
  min.start = apply(peak.mapping.table, 1, function(x) min(as.numeric(x[3]), as.numeric(x[4])))
  peak.mapping.table %>% mutate(Start_Min = min.start) -> peak.mapping.table
  max.end = apply(peak.mapping.table, 1, function(x) max(as.numeric(x[5]), as.numeric(x[6])))
  peak.mapping.table %>% mutate(End_Max = max.end) -> peak.mapping.table

  ## Add updated peak name to the table
  peak.mapping.table %>% mutate(MergedPeak = sprintf("%s:%d-%d:%s", Gene, Start_Min, End_Max, Strand),
                                Data1_Peak = table1.merge1.subset$Peak, Data2_Peak = table1.merge1.subset$Data2_Closest_Peak) ->
    peak.mapping.table

  ## Remove duplicated merged peaks
  peak.mapping.table %>% distinct(MergedPeak, .keep_all = TRUE) -> peak.mapping.table.unique

  ## Add all the non-overlapping peaks
  table1.distinct.subset = subset(apa.similarity.table, Matched_Similarity == FALSE)
  distinct.peaks = table1.distinct.subset$Peak

  ## Combined the merged and unique peaks
  combined.peaks = c(peak.mapping.table.unique$MergedPeak, table1.distinct.subset$Peak)

  return(combined.peaks)
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
#' @param apa.dataset.list a list of peak data-sets
#' @param self.merged.peaks.list the set of self-merged peaks from the reference data-set
#' @param sim.thresh The required similarity threshold for merging (default: 0.75)
#' @param allow.match.var The allowance for deviation from the sim.thresh for comparison peaks (default: 0.25)
#' @return a data-frame containing peaks, their class (merged or unique) and the original peak from the reference
#' @examples
#' generate_merged_peak_table(dataset.1, apa.dataset.list, self.merged.peaks.list)
generate_merged_peak_table <- function(dataset.1, apa.dataset.list, self.merged.peaks.list,
                                       sim.thresh = 0.75, allow.match.var = 0.25) {
  peak.set1 = self.merged.peaks.list[[dataset.1]]

  gene.names1 = sub("(*):.*", "\\1", peak.set1)

  ## remove any empty string gene names
  combined.peak.similarity.table = data.frame(Peak = peak.set1, stringsAsFactors = FALSE)
  rownames(combined.peak.similarity.table) = combined.peak.similarity.table$Peak

  similarity.labels = c()
  rev.similarity.labels = c()
  closest.peak.labels = c()
  match.labels = c()
  for (i in 1:length(apa.dataset.list)) {
    comparison.dataset = names(apa.dataset.list)[i]
    if (dataset.1 != comparison.dataset) {
      comparison.peak.set = self.merged.peaks.list[[comparison.dataset]]
      this.peak.similarity.table = generate_similarity_table(peak.set1, comparison.peak.set)
      table.subset = this.peak.similarity.table[combined.peak.similarity.table$Peak,
                                                c("Data2_Closest_Peak", "Data2_Similarity", "Data2_Data1_Similarity")]
      table.subset %>%
        mutate(Matched_Similarity =
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
    as.numeric(sub(".*:(.*)-.*:.*", "\\1", x))))
  colnames(start.positions) = c("Peak", closest.peak.labels)

  ## Go through all closest start sites and where matched similarities are true, take the minimum
  min.start.positions = vector("numeric", nrow(start.positions))
  for (i in 1:nrow(start.positions)) {
    start.set = start.positions[i, ][c(TRUE, as.logical(combined.peak.similarity.table.subset[i, match.labels]))]
    min.start.positions[i] = min(start.set)
  }

  ## Now calculate the maximum end positions for matched peaks
  end.positions = t(apply(combined.peak.similarity.table.subset[, c("Peak", closest.peak.labels)], 1, function(x)
    as.numeric(sub(".*:.*-(.*):.*", "\\1", x))))
  colnames(end.positions) = c("Peak", closest.peak.labels)

  ## Go through all closest start sites and where matched similarities are true, take the minimum
  max.end.positions = vector("numeric", nrow(end.positions))
  for (i in 1:nrow(end.positions)) {
    end.set = end.positions[i, ][c(TRUE, as.logical(combined.peak.similarity.table.subset[i, match.labels]))]
    max.end.positions[i] = max(end.set)
  }

  ## Now build a peak mapping table
  strand = sub(".*:.*-.*:(.*)", "\\1", combined.peak.similarity.table.subset$Peak)
  gene = sub("(.*):.*-.*:.*", "\\1", combined.peak.similarity.table.subset$Peak)

  peak.mapping.table = data.frame(OriginalPeak = combined.peak.similarity.table.subset$Peak, Gene = gene, Strand = strand,
                                  Start_Min = min.start.positions, End_Max = max.end.positions, stringsAsFactors = FALSE)


  ## Add updated peak name to the table
  peak.mapping.table %>% mutate(MergedPeak = sprintf("%s:%d-%d:%s", Gene, Start_Min, End_Max, Strand)) -> peak.mapping.table

  ## Remove duplicated merged peaks
  peak.mapping.table %>% distinct(MergedPeak, .keep_all = TRUE) -> peak.mapping.table.unique

  ## Now pull out the unique peaks
  peaks.unique = names(match.counts[which(match.counts == 0)])

  combined.peak.output.table = data.frame(Peak = c(peak.mapping.table.unique$MergedPeak, peaks.unique),
                                          Class = c(rep("Merged", nrow(peak.mapping.table.unique)),
                                                    rep(paste0("Unique_", dataset.1), length(peaks.unique))),
                                          OriginalPeak = c(peak.mapping.table.unique$OriginalPeak, peaks.unique),
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
#' @param apa.dataset.list a named list containing peaks called from different data-sets
#' @return a set of merged peaks
#' @examples
#' DoPeakMerging(apa.dataset.list)
do_peak_merging <- function(apa.dataset.list, sim.thresh = 0.75, allow.match.var = 0.25) {
  ### Generate the in-dataset similarity table
  self.merged.peaks.list = vector("list", length(apa.dataset.list))

  ## Go through the list and first generate self similarity tables to merge
  ## Return a list of peaks where similar peaks within each data-set have been merged
  for (i in 1:length(apa.dataset.list)) {
    dataset.name = names(apa.dataset.list)[i]
    this.peak.set = apa.dataset.list[[dataset.name]]
    print(paste0("Performing internal peak merging for ", dataset.name))
    this.peak.similarity.table = generate_self_similarity_table(this.peak.set)
    peak.set.merged = generate_self_merged_peaks(this.peak.similarity.table,
                                                 sim.thresh = sim.thresh, allow.match.var = allow.match.var)
    self.merged.peaks.list[[i]] = peak.set.merged
  }
  names(self.merged.peaks.list) = names(apa.dataset.list)

  ### Build a table of similarity between the reference data-set and remainder data-sets
  combined.merged.peaks.list = vector("list", length(apa.dataset.list))
  for (i in 1:length(apa.dataset.list)) {
    ref.dataset.name = names(apa.dataset.list)[i]
    print(paste0("Comparing peaks from ", ref.dataset.name, " to remaining data-sets"))
    this.merged.peak.table = generate_merged_peak_table(dataset.1 = ref.dataset.name, apa.dataset.list = apa.dataset.list,
                                                        self.merged.peaks.list = self.merged.peaks.list,
                                                        sim.thresh = sim.thresh, allow.match.var = allow.match.var)
    combined.merged.peaks.list[[i]] = this.merged.peak.table
  }
  names(combined.merged.peaks.list) = names(apa.dataset.list)

  ### Identif the peaks that are among all the merged and the ones that are unique to each data-set
  ## start with the merged peaks

  all.peaks = do.call(rbind, combined.merged.peaks.list)

  merged.peaks = subset(all.peaks, Class == "Merged")
  unique.peaks = subset(all.peaks, Class != "Merged")
  unique.peaks %>% distinct(Peak, .keep_all = TRUE) -> unique.peaks

  ## Remove duplicated merged peaks
  merged.peaks %>% distinct(Peak, .keep_all = TRUE) -> merged.peaks.unique

  ## Add the set of unique peaks that aren't already in the merged list
  unique.peaks = unique.peaks[which(!(unique.peaks$Peak %in% merged.peaks.unique$Peak)), ]

  final.merged.peaks = rbind(merged.peaks.unique, unique.peaks)

  return(final.merged.peaks)
}