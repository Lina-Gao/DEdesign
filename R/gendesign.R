#' Statistical design for RNA-seq experiments aimed at identifying differentially expressed genes
#'
#' Generates block design to assign samples to sequencing lanes and adaptors on an Illumina flow cell
#'
#'
#' @param trts a character vector for treatment groups
#' @param reps an integer vector for number of replicates for each element of \code{trts}. Length of \code{reps} has to match length of \code{trts}.
#' @param nperlane An integer for number of samples per lane for sequencing. When \code{nperlane} is not provided, a number between 4-6 will be selected based on number of treatments to generate efficient design
#' @param seed an integer initializing the random number generator. The default is \code{seed=1}
#'
#' @return \code{input} A list showing input parameters to the function
#' @return \code{Design} A list with two elements: \code{design} and \code{BlocksEfficiency}. \code{design} is a data frame giving lane and adaptor assignment for treatment groups. \code{BlocksEfficiency} is a data frame giving block efficiencies (D-Efficiencies) for lane and adaptor.
#' @return \code{suggestedDesign} When \code{nperlane} is \code{NULL} or same as the value would be selected to give better efficiencey, \code{suggestedDesign} is the same as \code{Design}. Otherwise, \code{suggestedDesign} uses different number of samples per lane from the input to give more efficient block design. \code{suggestedDesign} is a list with two elements: \code{design} and \code{BlocksEfficiency}. \code{design} is a data frame giving lane and adaptor assignment for treatment groups. \code{BlocksEfficiency} is a data frame giving block efficiencies (D-Efficiencies) for lane and adaptor.
#'
#'
#' @keywords RNA-seq, statistical experimental design, block design, Illumina flow cell
#'
#' @export
#'
#' @examples
#' gendesign(trts=letters[1:4], reps=rep(4,4), nperlane=NULL)
#' gendesign(trts=letters[1:4], reps=c(2,4,4,2), nperlane=3)
#' gendesign(trts=letters[1:6], reps=rep(3,6), nperlane=4)
#'
#' @import dplyr
#' @importFrom blocksdesign blocks
#'
#'
#'
#'
#'

gen.design <- function(trts, reps, nperlane = NULL, seed = 1) {
    if (class(trts) != "character") {
        stop("trts needs to be a character vector for treatment groups")
    }

    if (!((class(reps) == "integer" & all(reps > 0)) | (class(reps) == "numeric" & all(reps >
        0) & all(ceiling(reps) == floor(reps))))) {
        stop("reps needs to be positive integers")
    }

    if (!is.null(nperlane)) {
        if (!((class(nperlane) == "integer" & nperlane > 0) | (class(nperlane) == "numeric" &
            nperlane > 0 & ceiling(nperlane) == floor(nperlane)))) {
            stop("if nperlane is not NULL, it needs be a positive integer")
        }
    }

    if (length(trts) != length(reps)) {
        stop("trts and reps need to be the same length")
    }


    out <- list(input = list(trts = trts, nperlane = nperlane, reps = reps, seed = seed), Design = list(),
        suggestedDesign = list())
    ntrts <- length(trts)

    # suggested design, basically suggesting nperlane
    if (ntrts <= 6) {
        suggested.nperlane <- switch(ntrts, 4, 4, 3, 4, 5, 6)
    } else {
        # ntrts>6
        suggested.nperlane <- 4
        if (ntrts%%5 == 0)
            suggested.nperlane <- 5
        if (ntrts%%3 == 0 & ntrts%%4 != 0)
            suggested.nperlane <- 6
    }

    if ((sum(reps)%%suggested.nperlane) != 0) {
        # not filing a whole lane, add samples belong to ' ' to fill
        reps <- c(reps, (sum(reps)%%suggested.nperlane))
        trts <- c(trts, " ")
    }

    suggested.des <- blocks(treatments = rep(1, length(reps)), replicates = reps, seed = seed,
        rows = suggested.nperlane, columns = ceiling(sum(reps)/suggested.nperlane))
    out$suggestedDesign$design <- suggested.des$Design %>% transmute(chip = 1, lane = gsub("Cols_",
        "", Level_1.Cols) %>% as.integer, adp = gsub("Rows_", "", Level_1.Rows) %>% as.integer,
        trt = trts[Treatments])

    #this is for one chip design, need to modify for multiple chips
    out$suggestedDesign$BlocksEfficiency <- suggested.des$BlocksEfficiency[1:2,3:4]
    names( out$suggestedDesign$BlocksEfficiency) <- c("numBlocks","Efficiency")
    row.names( out$suggestedDesign$BlocksEfficiency) <- c("Adaptors","Lanes")


    if (!is.null(nperlane)) {

        if ((sum(reps)%%nperlane) != 0) {
            # not filing a whole lane, add samples belong to ' ' to fill
            reps <- c(reps, (sum(reps)%%nperlane))
            trts <- c(trts, " ")
        }
        des <- blocks(treatments = rep(1, length(reps)), replicates = reps, seed = seed, rows = nperlane,
            columns = ceiling(sum(reps)/nperlane))
        out$Design$design <- des$Design %>% transmute(chip = 1, lane = gsub("Cols_", "",
            Level_1.Cols) %>% as.integer, adp = gsub("Rows_", "", Level_1.Rows) %>% as.integer,
            trt = trts[Treatments])
        out$Design$BlocksEfficiency <- des$BlocksEfficiency[1:2,3:4]
        names(out$Design$BlocksEfficiency) <- c("numBlocks","Efficiency")
        row.names(out$Design$BlocksEfficiency) <- c("Adaptors","Lanes")
    } else {
        # nperlane is null

        out$Design <- out$suggestedDesign  #when nperlane is NULL, optimal design is the same as suggested design
    }
    out  #output
}




