#' Class transcript_archetypes
#'
#' Class \code{transcript_archetypes} holds lists of empirical transcript archetypes
#' and the matched data. Empirical profiles for flanking regions (width is user defined)
#' are included in the archetypes. Synthetic genome annotations, count data, and
#' transcript abundances can be simulated from this object and exported to standard
#' genomic formats.
#'
#' @slot transcripts a \code{\link[GenomicRanges]{GRanges-class}} that holds
#' all the transcript coordinates
#' @slot data a list containing run-length encodings of counts for both the sense and
#' anti-sense strand for the transcripts and their flanking regions
#' @slot abundance_estimate abundance estimates for each transcript
#' @slot configs a list of the configuration options used to generate this object
#'
#' @name transcript_archetypes-class
#' @rdname transcript_archetypes-class
#' @importClassesFrom GenomicRanges GRanges
#' @exportClass transcript_archetypes
methods::setClass("transcript_archetypes",
                  slots = c(transcripts = "GRanges",
                            data = "list",
                            abundance_estimate = "numeric",
                            configs = "list")
)

#' transcript_archetypes
#'
#' Contructs an object that holds the transcript archetypes models
#' @param transcripts a \link[GenomicRanges]{GRanges-class} object containing a
#' pre-filtered list of transcripts that show a clear signal that arises from only that
#' transcript with negligible signal from other sources
#' @param bigwig_plus polymerase density signal from the plus strand
#' @param bigwig_minus polymerase density signal from the minus strand
#' @param flank the length of the flanking region for which data is collected for
#' use in simulation
#' @param abundance_filter minimal abundance value for a viable archetype gene, all
#' transcripts with an estimated abundance less than this are removed (default = 100)
#' @param abundance_err_func Error function to use when estimating abundance (defaults
#' to identity)
#' @param mask a two element vector containing the number of bases to ignore at the head
#' and tail of the transcript when estimating abundance
#'
#' @return a \code{\link{transcript_archetypes-class}} object
#'
#' @export
transcript_archetypes <- function(transcripts, bigwig_plus, bigwig_minus, flank = 5e3,
                                  abundance_filter = 100,
                                  mask = c(1e3, 1e3),
                                  abundance_err_func = c("identity")) {
  # **Some checks prior to beginning construction**
  # Check for correct GRanges object metadata
  if (!is.GRanges(transcripts)) {
    stop("transcripts is not a GRanges object")
  }

  # Check that bigwig files exist
  if (!file.exists(bigwig_plus)) {
    stop("bigwig_plus file does not exist")
  }

  if (!file.exists(bigwig_minus)) {
    stop("bigwig_minus file does not exist")
  }

  # Check abundance filters
  if (!is.numeric(abundance_filter) | abundance_filter < 0 ) {
    stop("abundance_filter must be a numeric value > 0")
  }

  abundance_err_func <- abundance_err_func[1]
  if (!abundance_err_func %in%  c("identity")) {
    stop("invalid abundance_err_func specified")
  }

  # **End checks**

  # Construct query regions
  tmp <- construct_query_regions(transcripts, bigwig_plus, bigwig_minus, flank)
  transcripts <- tmp$tx
  query_regions <- tmp$q

  # Retrieve data for query regions on both strands
  pv_plus <- profile_views(bigwig_file = bigwig_plus, query_ranges = query_regions)
  pv_minus <- profile_views(bigwig_file = bigwig_minus, query_ranges = query_regions)

  # Extract the rle vectors and estimate abudance based on the sense strand, then keep
  # the entries that pass the abundance filter
  sense <- list()
  antisense <- list()
  abundance <- numeric(length(transcripts))
  sense_plus <- S4Vectors::decode(GenomicRanges::strand(transcripts) == "+")

  for(i in seq_along(transcripts)){
    if (sense_plus[i]) {
      sense[[i]] <- abs(pv_plus[i][[1]])
      antisense[[i]] <- abs(pv_minus[i][[1]])
    } else {
      sense[[i]] <- abs(S4Vectors::rev(pv_minus[i][[1]]))
      antisense[[i]] <- abs(S4Vectors::rev(pv_plus[i][[1]]))
    }
    abundance[i] <- estimate_abundance(sense[[i]], mask = mask, flank = flank,
                                       err_func = abundance_err_func)
  }

  final_transcripts <- transcripts[abundance >= abundance_filter]

  if( length(final_transcripts) == 0) {
    stop("No transcripts remaining after abundance filtering")
  }

  final_abundance <- abundance[abundance >= abundance_filter]
  final_qr <- query_regions[abundance >= abundance_filter]
  final_sense <- sense[abundance >= abundance_filter]
  final_antisense <- antisense[abundance >= abundance_filter]

  # Return transcript model object
  return(methods::new(Class = "transcript_archetypes",
                      transcripts = final_transcripts,
                      data = list(sense = final_sense,
                                  antisense = final_antisense),
                      abundance_estimate = final_abundance,
                      configs = list(bigwig_plus = bigwig_plus,
                                     bigwig_minus = bigwig_minus,
                                     flank = flank,
                                     abundance_filter = abundance_filter,
                                     mask = mask,
                                     abundance_err_func = abundance_err_func)
  ))
}

#' @title Estimate abundance
#'
#' @description Estimates abundance for transcript archetypes
#'
#' @param x A vector of counts
#' @param err_func error function (currently only identity supported)
#' @param mask a two element vector indicating the numbers of bases to skip at the head
#' and tail of the transcripts (respectively)
#' @param flank length of flanking regions present in the count vector
#'
#' @return A single abundance estimation
#'
#' @name estimate_abundance
estimate_abundance <- function(x, err_func, mask = c(0, 0), flank) {
  head_skip = mask[1] + flank
  tail_skip = mask[2] + flank
  if (err_func == "identity") {
    a <- mean(abs(x[(head_skip + 1):(length(x) - tail_skip)]))
  } else if (err_func == "log") {
    stop("log error not implemented yet")
  } else {
    stop("invalid err_func")
  }
  return(a)
}

#' @title Simulate annotations
#'
#' @description Simulates per loci annotation
#'
#' @param ta a \link{transcript_archetypes-class} object
#' @param n the number of loci to simulate annotations for
#' @param template can be a \code{\link[GenomicRanges]{GRanges-class}} that provides
#' template annotations as a basis for simulating loci with realistic properties.
#' @param method method to use for simulating intra-loci transcript properties. The two
#' methods are \code{neighbors} and \code{empirical}. More in information in details.
#' @param min_loci_tx minimum number of transcripts per loci
#' @param max_loci_tx maximum number of transcripts per loci
#' @param seed random seed to use for simulation
#' @param genome a genome name that will be used to retrieve chromosome lengths if
#' they are not already specified
#' @param tx_min_length do not simulate transcripts shorter than this length
#'
#' @details Stuff about methods
#'
#' @return A \link[GenomicRanges]{GRanges-class} of transcript annotations
#'
#' @name simulate_annotations
#' @rdname simulate_annotations
#'
#' @export
simulate_annotations <-  function(ta, n, template, genome = NULL, tx_min_length = 3e3,
                                  method = c("neighbor", "parametric"),
                                  min_loci_tx = 1, max_loci_tx = 40, seed = NULL) {
  # If seed is null choose a seed randomly
  if (is.null(seed)) {
    seed <- round(stats::runif(1) * 1e7)
  }
  set.seed(seed)

  # Must simulate at least 2 loci
  if (n < 2) {
    stop("n must be at least 2")
  }

  # Shorten method
  method <- method[1]

  # Chose style for seqinfo, give priority to ensembl, ncbi, and ucsc
  possible_styles <- GenomeInfoDb::seqlevelsStyle(template)
  possible_styles <-
    factor(possible_styles,
           levels = unique(c("Ensembl", "NCBI", "UCSC", possible_styles)))
  final_style <- as.character(sort(possible_styles)[1])

  # Check that all active seqlevels have a chromosomal length, if not a genome must
  # be specified that can be used to get those lengths
  if (any(is.na(GenomeInfoDb::seqlengths(template)))) {
    if (!is.null(genome)) {
      message("seqlengths not specified in annotation object, looking up seqlengths")
      genome_id <- GenomeInfoDb::mapGenomeBuilds(genome, "UCSC")$ucscID[1]
      if (!is.null(genome_id)) {
        if (utils::packageVersion("GenomeInfoDb") > 1.23) {
          chr_info <- GenomeInfoDb::getChromInfoFromUCSC(genome_id)
          # Create seqinfo for annotation
          si <- with(chr_info, GenomeInfoDb::Seqinfo(
            chrom, size, circular, genome_id))
        } else {
          chr_info <- GenomeInfoDb::fetchExtendedChromInfoFromUCSC(genome_id)
          # Create seqinfo for annotation
          si <- with(chr_info, GenomeInfoDb::Seqinfo(
            UCSC_seqlevel, UCSC_seqlength, circular, genome_id))
        }
      } else {
        stop("invalid genome")
      }
    } else {
      stop("invalid seqlengths specified without a genome being provided")
    }
  } else {
    si <- GenomicRanges::seqinfo(template)
  }

  # Handle various template forms
  if (is.GRanges(template)) {
    loci <- group_transcripts(template, distance = 0)
  } else {
    stop("template must be a GRanges object")
  }

  # Filter out short transcripts
  loci <- loci[GenomicRanges::width(loci) >= tx_min_length]

  # Filter for loci complexity
  loci <- loci[S4Vectors::elementNROWS(loci) >= min_loci_tx &
                 S4Vectors::elementNROWS(loci) <= max_loci_tx]

  # Create representation of loci that will be simulated with three vectors:
  # 1) List of which archetypes to use grouped by loci
  # 2) List of transcript offsets from 5' end of loci grouped by loci
  # 3) Loci strand
  # 4) The distance between loci

  # Now get #1 and #2 using either the neighbor method or the empirical method
  telomere_buffer <- 2e4 # buffer from tips of chromosome
  if(method == "neighbor") {
    target_loci <- loci[sample(seq_along(loci), size = n, replace = TRUE)]
    target_widths <- GenomicRanges::width(BiocGenerics::unlist(target_loci))
    archetype_widths <- GenomicRanges::width(ta@transcripts)
    # Calculate which archetype is closest to each transcript in length
    closest_archetype <- apply(as.matrix(target_widths), 1, function(x) {
      which.min(abs(x - archetype_widths))
    })

    # Split into per loci archetypes
    tx_archetypes <-
      S4Vectors::split(
        closest_archetype, rep(seq_len(n), S4Vectors::elementNROWS(target_loci)))
    # Get distance between the tss of the tx and start of loci
    loci_start <- GenomicRanges::flank(BiocGenerics::unlist(
      GenomicRanges::reduce(target_loci)), 1)
    tx_fiveprime_dist <- GenomicRanges::distance(BiocGenerics::unlist(target_loci),
                            loci_start[names(BiocGenerics::unlist(target_loci))])
    tx_fiveprime_dist <- split(tx_fiveprime_dist,
                               rep(seq_len(n), S4Vectors::elementNROWS(target_loci)))

  } else if(method == "empirical") {
    # Simulate inter-loci distances and loci strandedness
    stop("empirical method not yet implemented")
  } else {
    stop("Invalid simulation method specified")
  }

  # This is #3 & #4
  sim_loci_properties <- loci_property_sim(loci, n)
  loci_gaps <- sim_loci_properties$offset
  loci_strands <- sim_loci_properties$loci_strand

  # Generate annotations as GRanges object
  tx_gr <- annotate_loci(ta, tx_archetypes, tx_fiveprime_dist, loci_gaps, loci_strands,
                telomere_buffer, si)
  tx_gr <- GenomeInfoDb::sortSeqlevels(tx_gr)
  tx_gr <- GenomicRanges::sort(tx_gr, ignore.strand = T)
  GenomeInfoDb::seqlevelsStyle(tx_gr) <- final_style

  return(tx_gr)
}

#' @title Simulate loci level properties
#'
#' @description Simulates the distance between loci and the strandedness of those loci
#'
#' @param loci A \link[GenomicRanges]{GRangesList-class} of transcripts grouped by loci
#' @inheritParams simulate_annotations
#'
#' @return A list with two elements, the offsets between loci and the strands of the loci
#'
#' @name loci_property_sim
#' @rdname loci_property_sim
loci_property_sim <- function(loci, n) {
  # Calculate empirical distribution of distances between subsequent loci and whether
  # there is a change of strand or not
  loci_ranges <- unlist(GenomicRanges::reduce(loci))
  follow <- GenomicRanges::follow(loci_ranges, ignore.strand = TRUE, select = "all")
  # Add +1 to make it an offset as gap([1,2], [3,4]) = 0
  loci_offsets <- GenomicRanges::width(GenomicRanges::pgap(
    loci_ranges[S4Vectors::queryHits(follow)],
    loci_ranges[S4Vectors::subjectHits(follow)], ignore.strand = TRUE
  )) + 1
  loci_strand_swap <- GenomicRanges::strand(loci_ranges[S4Vectors::queryHits(follow)]) ==
    GenomicRanges::strand(loci_ranges[S4Vectors::subjectHits(follow)])

  # Simulate the distances between loci and their various strands
  offset_sample <- sample(seq_along(loci_offsets), size = n - 1, replace = TRUE)
  sim_offset <- loci_offsets[offset_sample]
  sim_switch <- BiocGenerics::as.vector(loci_strand_swap[offset_sample])

  # Choose random starting state
  init_state <- sample(c(0, 1), 1)
  loci_strand <- c(init_state, init_state + (cumsum(sim_switch) %% 2)) %% 2 + 1
  loci_strand <- c("+", "-")[loci_strand]
  return(list(offset = sim_offset, loci_strand = loci_strand))
}

#' @title Generate transcript annotations
#'
#' @description Generate a \link[GenomicRanges]{GRanges-class} of transcript annotations
#' from a variety of inputs describing the loci, the archetypes, and their relative
#' offsets
#'
#' @inheritParams simulate_annotations
#' @param tx_archetypes a list of which transcript archetypes go to each locus grouped
#' by locus
#' @param tx_fiveprime_dist a list of transcript offsets from 5' end of the locus grouped
#' by locus
#' @param loci_gaps a vector of offsets between loci
#' @param loci_strands a vector of the strandedness of each loci
#' @param telomere_buffer regions at the head and tail of each chromosome that do not
#' contain transcripts
#' @param seq_info a \code{seqinfo} object used to determine chromosome lengths
#'
#' @return A \link[GenomicRanges]{GRanges-class}
#'
#' @name annotate_loci
#' @rdname annotate_loci
annotate_loci <- function(ta, tx_archetypes, tx_fiveprime_dist, loci_gaps,
                          loci_strands, telomere_buffer, seq_info) {
  # Elements per loci use lengths()
  tx_per_loci <- unlist(lapply(tx_archetypes, length))
  # Get per trancript width
  tx_width <- GenomicRanges::width(ta@transcripts[unlist(tx_archetypes)])
  # Get max distance from tss to 3' end per transcript which is the total loci width
  loci_width <- unlist(lapply(S4Vectors::split(
    unlist(tx_fiveprime_dist, use.names = FALSE) + tx_width,
                   rep(seq_along(tx_per_loci), tx_per_loci)
    ), max
  ))

  # Get chromosomal lengths minus telomeres
  seq_lengths <- GenomeInfoDb::seqlengths(seq_info) - (2 * telomere_buffer)
  if (any(seq_lengths < 0)) {
    message("Removing ", sum(seq_lengths < 0), " chromosomes/scaffolds shorter than ",
            (2 * telomere_buffer), " from simulation")
    seq_lengths <- seq_lengths[seq_lengths > 0]
    seq_info <- seq_info[names(seq_lengths)]
  }

  cum_seq_lengths <- cumsum(as.numeric(seq_lengths))
  names(cum_seq_lengths) <- names(seq_lengths)

  # Get the end of each loci
  loci_end <- cumsum(c(0, loci_gaps) + loci_width)

  # Check if any loci go beyond total length of genome and remove them with a warning
  if (any(loci_end > cum_seq_lengths[length(cum_seq_lengths)])) {
    warning(sum(loci_end > cum_seq_lengths[length(cum_seq_lengths)]), " of ",
            length(loci_end), " loci go beyond the end of the genome and have been ",
            "removed. Consider simulating fewer loci.")
    keep_loci <- utils::tail(which(loci_end < cum_seq_lengths[length(cum_seq_lengths)]),
                             1)
    loci_end <- loci_end[seq_len(keep_loci)]
    tx_archetypes <- tx_archetypes[seq_len(keep_loci)]
    tx_fiveprime_dist <- tx_fiveprime_dist[seq_len(keep_loci)]
    loci_gaps <- loci_gaps[seq_len(keep_loci - 1)] # remove one off the end
    loci_strands <- loci_strands[seq_len(keep_loci)]
    tx_per_loci <- tx_per_loci[seq_len(keep_loci)]
    loci_width <- loci_width[seq_len(keep_loci)]
    # Get per trancript width
    tx_width <- GenomicRanges::width(ta@transcripts[unlist(tx_archetypes)])
  }

  # Convert strand to a per tx Rle for lookup
  tx_strands <- S4Vectors::Rle(loci_strands, tx_per_loci)

  # Adjust loci boundries to account for chromosomal lengths
  loci_seq_name <- cut(loci_end, breaks = c(0, cum_seq_lengths),
                       labels = names(cum_seq_lengths))
  loci_end <- loci_end - cum_seq_lengths[loci_seq_name] +
    seq_lengths[loci_seq_name] + telomere_buffer
  loci_start <- loci_end - loci_width

  # Rep loci end and start
  loci_end <- rep(loci_end, tx_per_loci)
  loci_start <- rep(loci_start, tx_per_loci)

  # Remove loci simulated on

  # Calculate the ranges of each transcript with the loci based offsets
  final_tx_start <- integer(length(tx_width))
  final_tx_end <- integer(length(tx_width))
  tx_plus <- BiocGenerics::as.vector(tx_strands  == "+")
  final_tx_start[tx_plus] <- loci_start[tx_plus] + unlist(tx_fiveprime_dist)[tx_plus]
  final_tx_start[!tx_plus] <- loci_end[!tx_plus] -
    (unlist(tx_fiveprime_dist)[!tx_plus] + tx_width[!tx_plus] - 1)
  final_tx_end[tx_plus] <- final_tx_start[tx_plus] + tx_width[tx_plus] - 1
  final_tx_end[!tx_plus] <- loci_end[!tx_plus] - unlist(tx_fiveprime_dist)[!tx_plus]

  # Create GRanges from transcript annotations
  # For now gene ids and loci ids are 1:1 but that may change in future
  tx_gr <- GenomicRanges::GRanges(
    rep(loci_seq_name, times = tx_per_loci),
    IRanges::IRanges(start = final_tx_start, end = final_tx_end),
    tx_strands,
    loci = rep(seq_along(tx_per_loci), tx_per_loci),
    gene_id = sprintf("G%08d", rep(seq_along(tx_per_loci), tx_per_loci)),
    transcript_id = sprintf("T%08d", seq_along(final_tx_start))
  )

  tx_gr <- tx_gr[GenomicRanges::start(tx_gr) > 0]
  GenomeInfoDb::seqinfo(tx_gr, pruning.mode=c("coarse")) <- seq_info

  return(tx_gr)
}

#' @title Simulate abundance values
#'
#' @description Take a \link[GenomicRanges]{GRanges-class} of transcript annotations
#' and return it annotated with abundance values of reach transcript
#'
#' @param annotations a \link[GenomicRanges]{GRanges-class} of transcript annotations
#' @param sampling_dist a function that generates a vector of random numbers greater than
#' or equal to zero. It also must have as an argument \code{n} which specifies the number
#' of observations to draw. (A number of the random number generators from the stats
#' package fits all these requirements; eg. rlnorm, rpois, rgamma, etc.). Additionally,
#' you may design your own custom function to pass here as long as it meets these
#' criteria.
#' @param seed random seed for reproducible sampling
#' @param ... any arguments to be passed to the \code{sampling_dist} function
#'
#' @return A \link[GenomicRanges]{GRanges-class} with a \code{score} column that holds
#' abundances
#'
#' @name simulate_abundances
#' @rdname simulate_abundances
#' @export
simulate_abundances <- function(annotations, sampling_dist = rgtex,
                                seed = NULL, ...) {
  # If seed is null choose a seed randomly
  if (is.null(seed)) {
    seed <- round(stats::runif(1) * 1e7)
  }
  set.seed(seed)

  # Check that sampling_dist is of a supported form
  if (class(sampling_dist) != "function") {
    stop("sampling_dist must be a function")
  }
  if (all(methods::formalArgs(sampling_dist) != "n")) {
    stop("sampling_dist function must have an argument \'n\' corresponding to the",
          "number of samples to be drawn")
  }
  # Check annotation is a GRanges
  if (!is.GRanges(annotations)) {
    stop("annotation must be a GRanges")
  }

  # Parse other args and extract any that should be
  args <- list(...)
  sampling_args <- formals(sampling_dist)
  sampling_args[intersect(methods::formalArgs(sampling_dist), names(args))] <-
    args[intersect(methods::formalArgs(sampling_dist), names(args))]
  sampling_args$n <- length(annotations)
  # Sample abundances
  annotations$score <- do.call(sampling_dist, sampling_args)

  if (any(annotations$abundance) < 0) {
    stop("abundances less than zero were generated from the simulating distribution")
  }
  return(annotations)
}

#' @title Simulate read data
#'
#' @description Take a \link[GenomicRanges]{GRanges-class} of transcript annotations
#' and a \code{transcript_archetypes} object and simulate read counts.
#' @importFrom utils setTxtProgressBar
#'
#' @param annotations a \code{GRanges} object with a \code{score} column indicating
#' abundances
#' @param ta a \link{transcript_archetypes-class} object
#'
#' @inheritParams simulate_annotations
#' @inheritParams simulate_abundances
#' @inheritParams resample_reads
#' @param show_progress show progress bar
#' @param library_depth library depth to simulate (approximate)
#'
#' @return A \link[GenomicRanges]{GRanges-class}
#'
#' @name simulate_data
#' @rdname simulate_data
#' @export
simulate_data <- function(ta, annotations, library_depth = 3e7, jitter = 0,
                          show_progress = F) {
  # Drop unused seqlevels
  GenomeInfoDb::seqlevels(annotations) <- GenomeInfoDb::seqlevelsInUse(annotations)

  # check that annotations have an score column
  if (!"score" %in% names(GenomicRanges::mcols(annotations))) {
    stop("annotations must have an \'score\' column")
  }

  # Check that all used seqinfos have  a length
  if (any(is.na(GenomeInfoDb::seqlengths(annotations)))) {
    stop("All seqlengths in use must have a finite, non-NA length")
  }

  # Create key from length to which archetype to use
  arch_key <- as.list(seq_along(ta@transcripts))
  names(arch_key) <- as.character(GenomicRanges::width(ta@transcripts))

  # Check that all widths in annotations have a matching width in the archetypes
  if (!all(as.character(GenomicRanges::width(ta@transcripts)) %in% names(arch_key))) {
    stop("There are annotations for which an archtype of a matched length does not ",
         "exist")
  }

  # Get mapping between transcripts and archetypes
  which_arch <- unlist(arch_key[as.character(GenomicRanges::width(annotations))])

  # Calculate reads per archetype
  arch_sense_read_depth <- unlist(lapply(ta@data$sense, sum))
  arch_antisense_read_depth <- unlist(lapply(ta@data$antisense, sum))

  # Get sampling depth for each annotation
  sense_sampling_depth <- (annotations$score / ta@abundance_estimate[which_arch]) *
    arch_sense_read_depth[which_arch]
  antisense_sampling_depth <- (annotations$score / ta@abundance_estimate[which_arch]) *
    arch_antisense_read_depth[which_arch]

  # Compute global scaling factor to get out correct total number of reads
  lib_size_scale_factor <- library_depth / sum(sense_sampling_depth, antisense_sampling_depth)
  sense_sampling_depth <- round(sense_sampling_depth * lib_size_scale_factor)
  antisense_sampling_depth <- round(antisense_sampling_depth * lib_size_scale_factor)

  # Iterate over seqlevels and simulate per chromosome
  ref_seqlengths <- GenomeInfoDb::seqlengths(annotations)
  flank_width <- ta@configs$flank
  plus <- IRanges::RleList()
  minus <- IRanges::RleList()
  if (show_progress) {
    pb <- utils::txtProgressBar(min = 0, max = sum(annotations$score > 0), style = 3)
  }
  pb_tracker <- 0
  for (seqlvl in names(ref_seqlengths)) {
    plus_sim <- integer(ref_seqlengths[seqlvl])
    minus_sim <- integer(ref_seqlengths[seqlvl])
    # Get list of transcripts with non-zero abundances
    anno_sub_idx <- which(as.vector(GenomicRanges::seqnames(annotations) == seqlvl) &
                            annotations$score > 0)
    # Iterate over transcripts with non-zero abundance and simulate data from archetypes
    for(i in anno_sub_idx) {
      adj_start <- GenomicRanges::start(annotations[i]) - flank_width
      adj_end <- GenomicRanges::end(annotations[i]) + flank_width
      strnd <- as.character(GenomicRanges::strand(annotations[i]))
      if (strnd == "+") {
        plus_sim[seq.int(adj_start, adj_end)] <-
          plus_sim[seq.int(adj_start, adj_end)] +
          resample_reads(ta@data$sense[[which_arch[i]]], sense_sampling_depth[i],
                         replace = T, jitter = jitter)
        minus_sim[seq.int(adj_start, adj_end)] <-
          minus_sim[seq.int(adj_start, adj_end)] -
          resample_reads(ta@data$antisense[[which_arch[i]]], antisense_sampling_depth[i],
                         replace = T, jitter = jitter)
      } else {
        minus_sim[seq.int(adj_end, adj_start)] <-
          minus_sim[seq.int(adj_end, adj_start)] -
          resample_reads(ta@data$sense[[which_arch[i]]], sense_sampling_depth[i],
                         replace = T, jitter = jitter)
        plus_sim[seq.int(adj_end, adj_start)] <-
          plus_sim[seq.int(adj_end, adj_start)] +
          resample_reads(ta@data$antisense[[which_arch[i]]],
                         antisense_sampling_depth[i], replace = T, jitter = jitter)
      }
      pb_tracker <- pb_tracker + 1
      if (show_progress) {
        setTxtProgressBar(pb, value = pb_tracker)
      }
    }
    plus[[seqlvl]] <- S4Vectors::Rle(plus_sim)
    minus[[seqlvl]] <- S4Vectors::Rle(minus_sim)
  }
  if (show_progress) {
    close(pb)
  }
  return(list(plus = plus, minus = minus))
}

#' @title Resample reads from 5' coverage vector
#'
#' @description Resample a vector of read coverage where each read only covers a single
#' point.
#'
#' @param x a numeric or \code{Rle} vector
#' @param size number of reads to sample
#' @param replace to sample with replacement (default: TRUE)
#' @param jitter amount by which reads can be shifted to the left or right (default: 0)
#'
#' @return A resampled numeric vector of 5' read counts
#'
#' @name resample_reads
resample_reads <- function(x, size, replace = T, jitter = 0) {
  out <- integer(length(x))
  pool <- which(S4Vectors::decode(x > 0))
  if (jitter == 0) {
    pos <- table(
      pool[sample.int(length(pool), size, replace = replace,
                      prob = S4Vectors::decode(x[x > 0]))]
    )
  } else {
    pos <- pool[sample.int(length(pool), size, replace = replace,
                    prob = S4Vectors::decode(x[x > 0]))] +
      sample.int(2 * jitter + 1, size, replace = T) - (jitter + 1)
    pos[pos < 1] <- 1
    pos[pos > length(out)] <- length(out)
    pos <- table(pos)
  }
  out[as.integer(names(pos))] <- as.integer(pos)
  return(out)
}

#' @title Exports simulation annotations and data
#'
#' @description Exports simulation annotations and data to bed and bigWig files
#' respectively
#'
#' @inheritParams simulate_data
#' @param data a two element list containing a \code{plus} and \code{minus}
#' \code{\link[IRanges]{RleList}} simulated from the \link{simulate_data}
#' function
#' @param directory directory to write output to
#' @param simulation_id character string to name simulation outputs with
#'
#' @return A resampled numeric vector of 5' read counts
#' @importFrom data.table :=
#' @name resample_reads
#' @export
export_simulation <- function(annotations, data, ta, directory = ".",
                              simulation_id = "sim") {
  ## Define data.table variables locally so that there is no R CMD CHECK complaint
  start <- NULL
  . <- NULL
  seqnames <- NULL
  end <- NULL
  strand <- NULL
  gene_id <- NULL
  transcript_id <- NULL
  sim_dat <- NULL

  ## Back to everything else
  dir.create(directory, showWarnings = F, recursive = T)
  bed_path <- file.path(directory, paste0(simulation_id, ".bed.gz"))
  plus_path <- file.path(directory, paste0(simulation_id, "_plus.bw"))
  minus_path <- file.path(directory, paste0(simulation_id, "_minus.bw"))
  ta_path <- file.path(directory, paste0(simulation_id, "_archetypes.RData"))
  bed <- data.table::as.data.table(annotations)
  bed[, start := start - 1]
  bed <- bed[, .(seqnames, start, end, '.', bed$score , strand, gene_id, transcript_id)]
  data.table::fwrite(file = bed_path, x = bed, compress = "gzip", sep = "\t",
                     col.names = F)
  rtracklayer::export.bw(sim_dat$plus, plus_path)
  rtracklayer::export.bw(sim_dat$minus, minus_path)
  saveRDS(ta, file = ta_path)
}

#' @title View archetype
#'
#' @description Plots archetype in sense/antisense orientation
#'
#' @inheritParams simulate_data
#' @param idx number id of archetype
#'
#' @name view_archetype
#' @export
view_archetype <- function(ta, idx) {
  # Create GRanges for data tracks
  sense <- iranges_to_granges(ta@data$sense[[idx]])
  antisense <- iranges_to_granges(ta@data$antisense[[idx]])
  GenomicRanges::score(antisense) <- -GenomicRanges::score(antisense)
  options(ucscChromosomeNames=FALSE)
  # Create scale track that labels flanking regions
  region_widths <- c(ta@configs$flank,
                     length(ta@data$sense[[idx]]) - 2 * ta@configs$flank,
                     ta@configs$flank)
  region_iranges <- IRanges::IRanges(start = cumsum(region_widths) - region_widths,
                   width = region_widths,
                   names = c("Flank", "Body", "Flank"))
  # Create plotting tracks
  axisTrack <- Gviz::GenomeAxisTrack(range = region_iranges)
  sense_track <- Gviz::DataTrack(range = sense, type = "h", name = "sense", col = "blue")
  antisense_track <- Gviz::DataTrack(range = antisense, type = "h", name= "antisense",
                                     col = "red")
  Gviz::plotTracks(list(axisTrack, sense_track, antisense_track), showId = TRUE)
}

iranges_to_granges <- function(x, chrom = "chr1", start = 0, strand = "+") {
  start_v <- (cumsum(S4Vectors::runLength(x)) - S4Vectors::runLength(x)) + start
  gr <- GenomicRanges::GRanges(chrom,
                               IRanges::IRanges(start = start_v,
                                                width = S4Vectors::runLength(x)),
                               strand = strand, score = S4Vectors::runValue(x)
  )
  return(gr)
}

