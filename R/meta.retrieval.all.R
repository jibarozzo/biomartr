#' @title Perform Meta-Genome Retrieval of all organisms in all kingdoms of life
#' @description Download genomes, proteomes, cds, gff, rna, or assembly stats
#' files of individual species of all kingdoms of life.
#' @inheritParams meta.retrieval
#' @author Hajk-Georg Drost
#' @details This function aims to perform bulk retrieval of all genomes
#' of species for all kingdoms of life.
#' @examples
#' \dontrun{
#' # download all genomes from refseq
#' meta.retrieval.all(db = "refseq", type = "genome")
#' # download all vertebrate genomes from genbank
#' meta.retrieval.all(db = "genbank", type = "genome")
#' # download all vertebrate genomes from ensemblgenomes
#' meta.retrieval.all(db = "genbank", type = "ensemblgenomes")
#' }
#' @return a character vector storing the file paths of the retrieved files.
#' @family meta_retrival
#' @export
meta.retrieval.all <- function(db = "refseq", type = "genome", reference = FALSE,
                               release = NULL,
                               remove_annotation_outliers = FALSE,
                               analyse_genome = FALSE,
                               assembly_type = "toplevel",
                               skip_bacteria = FALSE,
                               mute_citation = FALSE) {
    message("Starting ", type, " meta retrieval process of all species individually from database: ", db," ...")
    # retrieve all files of filetype from all kingdoms of life
    paths <- unlist(lapply(getKingdoms(db = db),
                           function(x) meta.retrieval(x, type = type, db = db,
                                                      group = NULL,
                                                      reference = reference,
                                                      release = NULL,
                                                      remove_annotation_outliers = remove_annotation_outliers,
                                                      analyse_genome = analyse_genome,
                                                      assembly_type = assembly_type,
                                                      skip_bacteria = skip_bacteria,
                                                      mute_citation = TRUE)))
    message("Meta retrieval process... finished!")
    please_cite_biomartr(mute_citation)
    return(paths)
}
