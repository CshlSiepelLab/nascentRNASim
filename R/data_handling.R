#' @title Assert chromosome exists
#'
#' @description  Checks whether a chromosome exists in the target bigwig. Throws
#' error if it does not.
#' @param bigwig_file a string pointing to a bigwig file
#' @param chromosome chromosome name
#' @rdname assert_chromosome_exists
#' @return silent if true, else error
#' @export
assert_chromosome_exists <- function(chromosome, bigwig_file) {
  bw_seqnames <- GenomicRanges::seqnames(
    rtracklayer::seqinfo(rtracklayer::BigWigFile(bigwig_file))
  )
  pass <- all(chromosome %in% bw_seqnames)
  if (!pass) {
    not_incl <- setdiff(chromosome, bw_seqnames)
    stop(paste("chromosome",
               paste(not_incl, collapse = ", "),
               "does not exist in", bigwig_file))
  }
}

#' Class profile_views
#'
#' Class \code{profile_views} holds views of ranges on a dataset.
#'
#' @slot query_ranges a \code{\link[GenomicRanges]{GRanges-class}} that holds
#' all the query ranges
#' @slot profiles a \code{\link[IRanges]{ViewsList-class}} object that holds the data
#' views
#' @slot idx a lookup index used by acessors for the class
#'
#' @name profile_views-class
#' @rdname profile_views-class
#' @importClassesFrom GenomicRanges GRanges
#' @importClassesFrom IRanges ViewsList
#' @exportClass profile_views
methods::setClass("profile_views",
                  slots = c(query_ranges = "GRanges",
                            profiles = "ViewsList",
                            idx = "matrix")
)

#' @title extract regions from file
#'
#' @description  Extracts regions from a bigwig over a set of ranges from a
#' \code{\link[GenomicRanges]{GRanges-class}} object and build a
#' \code{\link{profile_views}} object
#'
#' @inheritParams assert_chromosome_exists
#' @param query_ranges A \code{\link[GenomicRanges]{GRanges-class}} object
#'
#' @importClassesFrom GenomicRanges GRanges
#'
#' @return A list of run length encodings
#'
#' @name profile_views
#' @rdname profile_views
#'
#' @export
profile_views <- function(bigwig_file, query_ranges) {
  # Restrict to active seqnames
  GenomeInfoDb::seqlevels(query_ranges) <- GenomeInfoDb::seqlevelsInUse(query_ranges)

  # Get unique chromosome names in bins
  chrom <-
    unique(S4Vectors::runValue(GenomicRanges::seqnames(query_ranges)))

  # Check that all query chromosomes are present in bigwig
  assert_chromosome_exists(chromosome = chrom, bigwig_file)

  # Split queries into list
  names(query_ranges) <- seq_along(query_ranges)
  query_list <- split(IRanges::IRanges(GenomicRanges::start(query_ranges),
                                       GenomicRanges::end(query_ranges),
                                       names = names(query_ranges)),
                      f = GenomicRanges::seqnames(query_ranges))

  # Import bigwig, only chromosomes to be queried
  imported_bw <- rtracklayer::import.bw(
    con = bigwig_file, which = query_ranges,
    as = "RleList")[names(query_list)]

  # Iterate over chromosomes
  tx_views <- IRanges::Views(
    subject = imported_bw,
    query_list)

  # Extract the originating row and final location in views object
  origin <- as.integer(unlist(lapply(query_list, names)))
  list_idx <- rep(seq_along(query_list), unlist(lapply(query_list, length)))
  row_idx <- unlist(lapply(query_list, seq_along))

  lookup_idx <- cbind(origin = origin, list_idx = list_idx, row_idx = row_idx)
  rownames(lookup_idx) <- NULL
  lookup_idx <- lookup_idx[order(lookup_idx[, 1]), ]

  return(methods::new(Class = "profile_views",
                      query_ranges = query_ranges,
                      profiles = tx_views,
                      idx = lookup_idx
  ))
}

#' @title subsetting
#'
#' @description Convenience subsetting for views object. Retrieves view based on index of
#' query range.
#'
#' @param x \link{profile_views-class} object
#' @param i subset entries of views object
#' @param j unsupported
#' @param drop unsupported
#' @param ... just for compatibility with S4 method
#'
#' @aliases [
#' @return A views object
#' @export
methods::setMethod("[", c(x = "profile_views", i = "numeric", j = "missing", drop = "ANY"),
          function(x, i, j, ..., drop=TRUE){
            i <- as.integer(i)
            return(x@profiles[[x@idx[i, 2]]][x@idx[i, 3]])
          }

)

#' @title show
#'
#' @description Display object
#'
#' @inheritParams methods::show
#' @importMethodsFrom methods show
#'
#' @export
methods::setMethod("show", c(object = "profile_views"),
          function(object) {
            print(paste("profile_views object with", length(object@query_ranges),
                        "views"))
          }

)
