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

#' is.GRangesList
#'
#' Checks if object is a \code{\link[GenomicRanges]{GRangesList-class}} object
#'
#' @return logical
#'
#' @export
is.GRangesList <- function(x) {
  return(class(x) == "GRangesList")
}
