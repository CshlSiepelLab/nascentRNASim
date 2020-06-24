#' @title List available archetype sets
#'
#' @description Lists archetypes sets that are included in this package
#' @return A character vector
#'
#' @export
list_archetype_set <- function() {
  adir <- system.file("extdata", "archetype_sets", package = "nascentRNASim")
  return(dir(adir))
}

#' @title Load archetype set included in package
#'
#' @description Load archetype set included in package. To see available sets run
#' \code{list_archetype_set()}.
#' @param id the id of the archetype set to load
#'
#' @return A \code{transcript_archetypes} object
#'
#' @export
load_archetype_set <- function(id = "celastrol_dukler2017") {
  archive_dir <- system.file("extdata", "archetype_sets", package = "nascentRNASim")
  id_path <- file.path(archive_dir, id)
  if (!dir.exists(id_path)) {
    stop("Invalid id, run list_archetype_set() to view valid ids")
  }

  # Load archetype annotations
  bed_path <- file.path(id_path, paste0(id, "_archetypes.bed"))
  bed <- rtracklayer::import.bed(bed_path)

  # Load bigwigs
  bwp_path <- file.path(id_path, paste0(id, "_plus.bw"))
  bwm_path <- file.path(id_path, paste0(id, "_minus.bw"))

  # Create transcript_archetype object
  ta <- transcript_archetypes(transcripts = bed, bigwig_plus = bwp_path,
                              bigwig_minus = bwm_path,
                              flank = 6e3,
                              abundance_filter = 0.05,
                              mask = c(1e3, 1e3))
  return(ta)
}


#' @title Save an archetype set
#'
#' @description Save subset of data required to make an archetype set.
#'
#' @param bed_path path to bed file that holds archetype regions
#' @param bigwig_plus path to bigwig for plus strand
#' @param bigwig_minus path to bigwig for minus strand
#' @param out_dir directory to save files in
#' @param flank_length length of region flanking those in bed file to save in the
#' subsetted bigwig files
#'
#' @return NULL
#'
#' @export
store_reduced_archetype_bw <- function(bed_path, bigwig_plus, bigwig_minus,
                                       out_dir = ".", flank_length = 2e4) {
  # Read in bed file
  bed <- rtracklayer::import.bed(bed_path) + flank_length

  # Get id from bed file
  id <- gsub("_archetypes.bed", "",basename(bed_path))

  # Construct bigwig output filenames
  dir.create(out_dir, showWarnings = F, recursive = T)
  out_bwp <- file.path(out_dir, paste0(id, "_plus.bw"))
  out_bwm <- file.path(out_dir, paste0(id, "_minus.bw"))

  style <- GenomeInfoDb::seqlevelsStyle(rtracklayer::BigWigFile(bigwig_plus))[1]
  GenomeInfoDb::seqlevelsStyle(bed) <- style

  # Import subset of bigwig files
  imported_bwp <- rtracklayer::import.bw(
    con = bigwig_plus, which = bed,
    as = "RleList")
  imported_bwm <- rtracklayer::import.bw(
    con = bigwig_minus, which = bed,
    as = "RleList")

  # Export bigwigs
  rtracklayer::export.bw(imported_bwp, con = out_bwp)
  rtracklayer::export.bw(imported_bwm, con = out_bwm)
}
