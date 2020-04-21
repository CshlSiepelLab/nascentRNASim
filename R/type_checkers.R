#' is.GRanges
#'
#' Checks if object is a \code{\link[GenomicRanges]{GRanges-class}} object
#'
#' @return logical
#'
#' @export
is.GRanges <- function(x) {
  return(class(x) == "GRanges")
}
