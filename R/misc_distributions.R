#' @title Simulate from a zero-inflated exponential distribution
#'
#' @description Zero-inflated exponential distribution
#'
#' @param n number of observations. If length(n) > 1, the length is taken to be the
#' number required.
#' @param rate vector of rates
#' @param zero_weight probability of drawing a 0
#'
#' @return Numeric vector of length \code{n}
#'
#' @name rzexp
#' @rdname rzexp
#' @export
rzexp <- function(n, rate = 1, zero_weight = 0) {
  if (length(n) > 1) {
    n <- length(n)
  }
  x <- numeric(n)
  replace <- (stats::runif(n = n) > zero_weight)
  x[replace] <- stats::rexp(sum(replace), rate = rate)
  return(x)
}

#' @title List available density estimates from GTEx
#'
#' @description Lists available density estimates from GTEx for use with \code{rgtex}
#'
#' @param pattern tissue name pattern to match
#'
#' @return character vector of tissue names
#'
#' @name list_gtex_density
#' @rdname list_gtex_density
#' @export
list_gtex_density <- function(pattern = NULL) {
  kernel_dir <- system.file("extdata", "gtex_expression_kernels",
                            package = "nascentRNASim")
  files <- dir(path = kernel_dir, pattern = pattern)
  tissues <- gsub(".kernel.txt.gz","",files)
  return(tissues)
}

#' @title Sample from density estimates from GTEx
#'
#' @description Sample from probability density functions estimated from sex averaged
#' tissue specific transcript TPM from GTEx (with outliers removed).To view available
#' tissues use \code{list_gtex_density()}.
#'
#' @param n number of observations. If length(n) > 1, the length is taken to be the
#' number required.
#' @param tissue the tissue to simulate from (default: skeletal_muscle)
#'
#' @return Numeric vector of length \code{n}
#'
#' @name rgtex
#' @rdname rgtex
#' @export
rgtex <- function(n, tissue = "skeletal_muscle") {
  if (length(n) > 1) {
    n <- length(n)
  }
  if (length(tissue) > 1) {
    stop("Only one tissue may be specified at a time")
  }
  kernel_dir <- system.file("extdata", "gtex_expression_kernels",
                            package = "nascentRNASim")
  file_name <- paste0(tissue, ".kernel.txt.gz")
  density_kernel_file <- file.path(kernel_dir, file_name)
  if (!file.exists(density_kernel_file)) {
    stop(density_kernel_file, "does not exist. Please check that valid tissue name",
         "is specified")
  }
  dens <- data.table::fread(density_kernel_file)
  x <- sample(dens$x, size = n, replace = T, prob = dens$y)
  return(x)
}
