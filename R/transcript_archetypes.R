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
  sense_plus <- BiocGenerics::as.vector(GenomicRanges::strand(transcripts) == "+")

  for(i in seq_along(transcripts)){
    if (sense_plus[i]) {
      sense[[i]] <- as.vector(pv_plus[i])[[1]]
      antisense[[i]] <- as.vector(pv_minus[i])[[1]]
    } else {
      sense[[i]] <- S4Vectors::rev(as.vector(pv_minus[i]))[[1]]
      antisense[[i]] <- S4Vectors::rev(as.vector(pv_plus[i]))[[1]]
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
#' @param tx A \code{\link[GenomicRanges]{GRanges-class}} or
#' \code{\link[GenomicFeatures]{TxDb-class}}
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
