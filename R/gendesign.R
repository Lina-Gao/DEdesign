#' Statistical design for RNA-seq experiments aimed at identifying differentially expressed genes
#'
#' Generates block design to assign samples to sequencing lanes and adaptors on an Illumina flow cell
#'
#'
#' @param treatments a character vector for treatment groups. If the experiment consists of two or more factors, provide a character vector representing all possible combinations of levels across these factors.
#' @param replicates an integer vector for number of replicates for each element of \code{treatments}. Length of \code{replicates} has to match length of \code{treatments}.
#' @param nperlane An integer for number of samples per lane for sequencing. When \code{nperlane} is not provided, a number between 4-6 will be selected based on number of treatments to generate efficient design
#' @param seed an integer initializing the random number generator. The default is \code{seed=1}
#'
#' @return \code{input} A list showing input parameters to the function
#' @return \code{Design} A list with two elements: \code{design} and \code{BlocksEfficiency}. \code{design} is a data frame giving lane and adaptor assignment for treatment groups. \code{BlocksEfficiency} is a data frame giving block efficiencies (D-Efficiencies) for lane and adaptor.
#' @return \code{suggestedDesign} When possible, \code{suggestedDesign} uses different number of samples per lane from the input to give more efficient block design; otherwise it is the same as \code{Design}.  \code{suggestedDesign} is a list with two elements: \code{design} and \code{BlocksEfficiency}. \code{design} is a data frame giving lane and adaptor assignment for treatment groups. \code{BlocksEfficiency} is a data frame giving block efficiencies (D-Efficiencies) for lane and adaptor.
#'
#'
#' @keywords RNA-seq, statistical experimental design, block design, Illumina flow cell
#'
#' @export
#'
#' @examples
#' gendesign(treatments=letters[1:4], replicates=rep(4,4), nperlane=NULL)
#' gendesign(treatments=letters[1:4], replicates=c(2,4,4,2), nperlane=3)
#' gendesign(treatments=letters[1:6], replicates=rep(3,6), nperlane=4)
#'
#' @import dplyr
#' @importFrom blocksdesign blocks
#'
#'
#'
#'
#'

gendesign <- function(treatments, replicates, nperlane = NULL, seed = 1) {
    if (class(treatments) != "character") {
        stop("treatments needs to be a character vector for treatment groups")
    }

    if (!((class(replicates) == "integer" & all(replicates > 0)) | (class(replicates) == "numeric" & all(replicates >
        0) & all(ceiling(replicates) == floor(replicates))))) {
        stop("replicates needs to be positive integers")
    }

    if (!is.null(nperlane)) {
        if (!((class(nperlane) == "integer" & nperlane > 0) | (class(nperlane) == "numeric" &
            nperlane > 0 & ceiling(nperlane) == floor(nperlane)))) {
            stop("if nperlane is not NULL, it needs be a positive integer")
        }
    }

    if (length(treatments) != length(replicates)) {
        stop("treatments and replicates need to be the same length")
    }


    out <- list(input = list(treatments = treatments, nperlane = nperlane, replicates = replicates, seed = seed), Design = list(),
        suggestedDesign = list())
    ntreatments <- length(treatments)

    # suggested design, basically suggesting nperlane
    if (ntreatments <= 6) {
        suggested.nperlane <- switch(ntreatments, 4, 4, 3, 4, 5, 6)
    } else {
        # ntreatments>6
        suggested.nperlane <- 4
        if (ntreatments%%5 == 0)
            suggested.nperlane <- 5
        if (ntreatments%%3 == 0 & ntreatments%%4 != 0)
            suggested.nperlane <- 6
    }

    if ((sum(replicates)%%suggested.nperlane) != 0) {
        # not filing a whole lane, add samples belong to ' ' to fill
        replicates <- c(replicates, (sum(replicates)%%suggested.nperlane))
        treatments <- c(treatments, " ")
    }

    suggested.des <- blocks(treatments = rep(1, length(replicates)), replicates = replicates, seed = seed,
        rows = suggested.nperlane, columns = ceiling(sum(replicates)/suggested.nperlane))
    out$suggestedDesign$design <- suggested.des$Design %>% transmute(flowcell = 1, lane = gsub("Cols_",
        "", Level_1.Cols) %>% as.integer, adaptor = gsub("Rows_", "", Level_1.Rows) %>% as.integer,
        treatment = treatments[Treatments]) %>% arrange(flowcell,lane,treatment,adaptor)

    #this is for one flowcell design, need to modify for multiple flowcells
    out$suggestedDesign$BlocksEfficiency <- suggested.des$BlocksEfficiency[1:2,3:4]
    names( out$suggestedDesign$BlocksEfficiency) <- c("numBlocks","Efficiency")
    row.names( out$suggestedDesign$BlocksEfficiency) <- c("Adaptors","Lanes")


    if (!is.null(nperlane)) {

        if ((sum(replicates)%%nperlane) != 0) {
            # not filing a whole lane, add samples belong to ' ' to fill
            replicates <- c(replicates, (sum(replicates)%%nperlane))
            treatments <- c(treatments, " ")
        }
        des <- blocks(treatments = rep(1, length(replicates)), replicates = replicates, seed = seed, rows = nperlane,
            columns = ceiling(sum(replicates)/nperlane))
        out$Design$design <- des$Design %>% transmute(flowcell = 1, lane = gsub("Cols_", "",
            Level_1.Cols) %>% as.integer, adaptor = gsub("Rows_", "", Level_1.Rows) %>% as.integer,
            treatment = treatments[Treatments]) %>% arrange(flowcell,lane,treatment,adaptor)
        out$Design$BlocksEfficiency <- des$BlocksEfficiency[1:2,3:4]
        names(out$Design$BlocksEfficiency) <- c("numBlocks","Efficiency")
        row.names(out$Design$BlocksEfficiency) <- c("Adaptors","Lanes")
    } else {
        # nperlane is null

        out$Design <- out$suggestedDesign  #when nperlane is NULL, optimal design is the same as suggested design
    }
    out  #output
}




