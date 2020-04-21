#' @title Group transcripts
#'
#' @description
#' Creates groups of transcripts on the same strand based on
#' proximity. The groups are constructed of connected transcripts
#' such that a pair of transcripts are considered connected if
#' they are within a given distance \emph{d} of each other. The
#' group then consists of all transcripts that are connected to
#' at least one other member.
#'
#' @param transcript_granges \code{\link[GenomicRanges]{GRanges-class}}
#' object
#' @param distance the distance within which two transcripts are
#' considered connected
#'
#' @return A \code{\link[GenomicRanges]{GenomicRangesList-class}}
#' object
#'
#' @rdname group_transcripts
#' @export
group_transcripts <- function(transcript_granges, distance = 0) {
  # Check that txdb is an actual GRanges object
  if (!methods::is(transcript_granges, "GRanges")) {
    stop("transcript is not a GRanges object")
  }
  # Ensure everything is sorted, otherwise later assignment assumptions fail
  transcript_granges <- GenomeInfoDb::sortSeqlevels(transcript_granges)
  transcript_granges <- GenomicRanges::sort(transcript_granges)
  # Assign a unique id to each row
  transcript_granges$unique_id <- seq_len(length(transcript_granges))
  # Create expanded granges
  tx_granges_expand <- GenomicRanges::resize(
    transcript_granges,
    width = GenomicRanges::width(transcript_granges) + floor(distance / 2),
    fix = "end")
  tx_granges_expand <- GenomicRanges::resize(
    tx_granges_expand,
    width = GenomicRanges::width(tx_granges_expand) + ceiling(distance / 2),
    fix = "start")
  tx_granges_expand <- GenomicRanges::trim(tx_granges_expand)
  # Reduce granges to get unions of overlapping ranges which will become the
  # transcript groups. The -/+ 1 bit is to assure that perfectly adjoining
  # groups are not merged
  GenomicRanges::end(tx_granges_expand) <- GenomicRanges::end(tx_granges_expand) - 1
  tx_reduce <- GenomicRanges::reduce(tx_granges_expand, ignore.strand = FALSE)
  GenomicRanges::end(tx_reduce) <- GenomicRanges::end(tx_reduce) + 1
  # Get the group assignments
  group_assignment <- GenomicRanges::findOverlaps(query = transcript_granges,
                                                  subject = tx_reduce)
  # Create group vector string
  uid <- paste0(S4Vectors::subjectHits(group_assignment), "_",
                GenomicRanges::strand(transcript_granges))

  # Split the transcripting into their groups
  gr_groups <- GenomicRanges::split(transcript_granges, f = uid)
  return(gr_groups)
}

#' @title Construct query regions
#'
#' @description Given a set of transcripts constructs query regions and removes
#' transcripts whose query regions go out of bounds
#'
#' @param transcripts A \code{\link[GenomicRanges]{GRanges-class}} object
#' @inheritParams transcript_archetypes
#'
#' @return A single abundance estimation
#' @include transcript_archetypes.R
#'
#' @name estimate_abundance
construct_query_regions <- function(transcripts, bigwig_plus, bigwig_minus,
                                    flank = 5e3) {
  # Get chrom information from chromosomes and set seqlevel style appropriately
  bw_info_p <- rtracklayer::seqinfo(rtracklayer::BigWigFile(bigwig_plus))
  bw_info_m <- rtracklayer::seqinfo(rtracklayer::BigWigFile(bigwig_plus))
  GenomeInfoDb::seqlevelsStyle(transcripts) <-
    GenomeInfoDb::seqlevelsStyle(bw_info_p)

  # Get active seqlevels
  active_seqlevels <- GenomeInfoDb::seqlevelsInUse(transcripts)

  # Check that bigwig plus and minus have same chromosome lengths
  if (any(GenomeInfoDb::seqlengths(bw_info_p[active_seqlevels]) !=
          GenomeInfoDb::seqlengths(bw_info_p[active_seqlevels]))) {
    stop("bigwig_plus and bigwig_minus file were constructed with differing chromosome ",
         "length specification for one or more of the chromosomes included in the ",
         "query regions")
  }

  # Create expanded query regions
  query_regions <- transcripts
  GenomicRanges::start(query_regions) <- GenomicRanges::start(query_regions) - flank
  GenomicRanges::end(query_regions) <- GenomicRanges::end(query_regions) + flank

  # Remove query regions which go out-of-bounds
  len_vec <- GenomeInfoDb::seqlengths(bw_info_p[active_seqlevels])
  tmp <- as.vector(GenomicRanges::seqnames(query_regions))
  remove_tx <- c(which(GenomicRanges::end(query_regions) > len_vec[tmp]),
                 which(GenomicRanges::start(query_regions) < 1))

  if (length(remove_tx) > 0) {
    query_regions <- query_regions[-remove_tx]
    transcripts <- transcripts[-remove_tx]
    message("Removing ", length(remove_tx), " out of bounds regions")
  }
  return(list(tx = transcripts, q = query_regions))
}

#' @title Filter annotations
#'
#' @description Gets a filtered candidate list of transcripts to use for simulation
#'
#' @param tx A \code{\link[GenomicRanges]{GRanges-class}} or
#' \code{\link[GenomicFeatures]{TxDb-class}}
#' @param min_distance minimum distance between candidate transcript and closest
#' transcript (strand agnostic).
#' @param min_length minimum length of candidate transcript
#' @param keep_chromosomes a character vector of chromosomes to select candidate
#' transcripts from. Defaults to all chromsomes in \code{tx} object
#' @param keep_tx_biotype a character vector of transcript types to select candidate
#' transcripts from. Defaults to all types.
#'
#' @return A list of vectors with each one corresponding to one set of bins and
#' each element of a vector corresponding to a bin
#'
#' @name filter_annotations
#' @rdname filter_annotations
#'
#' @export
filter_annotations <- function(tx, min_distance = 2e4, min_length = 3e3,
                               keep_chromosomes = NULL, keep_tx_biotype = NULL) {
  # Check that tx is a GRanges or txdb object
  if (!class(tx) %in% c("GRanges", "TxDb")) {
    stop("tx must be a GRanges or Txdb object")
  }

  # Extract transcripts from a txdb
  if (class(tx) == "TxDb") {
    tx_gr <- GenomicFeatures::transcripts(tx, columns=c("tx_name", "gene_id",
                                                        "tx_type"))
  }

  # Check that tx_gr has tx_type column if tx_type filter included
  if(!is.null(keep_tx_biotype) &&
     ! "tx_type" %in% colnames(GenomicRanges::mcols(tx_gr))) {
    stop("If tx_biotype filtering option is used then tx must contain a tx_type column")
  }

  # Keep specified chromosomes using filter
  if (!is.null(keep_chromosomes)) {
    keep_chromosomes <- as.character(keep_chromosomes)
    user_style <- GenomeInfoDb::seqlevelsStyle(keep_chromosomes)
    GenomeInfoDb::seqlevelsStyle(tx_gr) <- user_style[1]
    GenomeInfoDb::seqlevels(tx_gr, pruning.mode="coarse") <- keep_chromosomes
  }

  # Get distance to nearest transcript for all transcripts, note that strand is ignored
  dist_gr <- GenomicRanges::distanceToNearest(tx_gr, ignore.strand = TRUE)
  keep <- S4Vectors::queryHits(
    dist_gr[S4Vectors::elementMetadata(dist_gr)$distance >= min_length])
  candidates <- tx_gr[keep]

  # Now filter the candidates by length and txtype if specified
  candidates <- candidates[width(candidates) >= min_length]

  # If a tx_type filter is specified, apply it now
  if (!is.null(keep_tx_biotype)) {
    candidates <- candidates[candidates$tx_type %in% keep_tx_biotype]
  }
  return(candidates)
}
