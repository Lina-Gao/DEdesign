#' Statistical design for RNA-seq experiments aimed at identifying differentially expressed genes
#'
#' Generates block design to assign samples to sequencing lanes and adapters for Illumina flow cell
#'
#'
#' @param treatments a character vector for treatment groups. If the experiment consists of two or more factors, provide a character vector representing all possible combinations of levels across these factors.
#' @param replicates an integer vector for number of replicates for each element of \code{treatments}. Length of \code{replicates} has to match length of \code{treatments}.
#' @param nperlane An integer for number of samples per lane for sequencing. When \code{nperlane} is not provided, a number between 4-6 will be selected based on number of treatments to generate efficient design
#' @param seed an integer initializing the random number generator. The default is \code{seed=1}. Designs can be rebuilt repeatedly using different seed to check that a near-optimum design has been found.
#' @param searches the maximum number of local optima searched at each stage of a treatment and block design optimization. The default depends on the design size. For optimum results, try large number of searches.
#'
#' @return \code{gendesign} returns an object of class "\code{DEdesign}", which is a list containing the following components:
#'
#' @return \code{input} A list showing input parameters to the function
#' @return \code{Design} A list with two elements: \code{design} and \code{BlocksEfficiency}. \code{design} is a data frame giving flowcell, lane and adapter assignment for treatment groups. \code{BlocksEfficiency} is a data frame giving block efficiencies (D-Efficiencies) for lane and adapter.
#' @return \code{suggestedDesign} When possible, \code{suggestedDesign} uses different number of samples per lane from the input to give more efficient block design; otherwise it is the same as \code{Design}.  \code{suggestedDesign} is a list with two elements: \code{design} and \code{BlocksEfficiency}. \code{design} is a data frame giving flowcell, lane and adapter assignment for treatment groups. \code{BlocksEfficiency} is a data frame giving block efficiencies (D-Efficiencies) for lane and adapter.
#'
#'
#' @keywords RNA-seq, statistical experimental design, block design, Illumina flow cell
#'
#'
#' @exportClass DEdesign
#'
#' @examples
#' gendesign(treatments=letters[1:4], replicates=rep(4,4), nperlane=NULL)
#' gendesign(treatments=letters[1:4], replicates=c(2,4,4,2), nperlane=3)
#' gendesign(treatments=letters[1:6], replicates=rep(3,6), nperlane=4)
#'
#' @import dplyr
#' @importFrom blocksdesign design
#'
#' @export
#'
gendesign <- function(treatments, replicates, nperlane = NULL, seed = 1, searches = NULL) {
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
            stop("if nperlane is not NULL, it needs to be a positive integer")
        }
    }

    if (length(treatments) != length(replicates)) {
        stop("treatments and replicates need to be the same length")
    }


    out <- structure(list(input = list(treatments = treatments, nperlane = nperlane, replicates = replicates, seed = seed, searches = searches),
                          Design = list(),
                          suggestedDesign = list()),
                     class = "DEdesign")

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
    # #don't need
    # if ((sum(replicates)%%suggested.nperlane) != 0) {
    #     # not filing a whole lane, add samples belong to ' ' to fill
    #     replicates <- c(replicates, (sum(replicates)%%suggested.nperlane))
    #     treatments <- c(treatments, " ")
    # }


    suggested.numlane <- ceiling(sum(replicates)/suggested.nperlane)
    trtforeachsample = as.factor(mapply(rep,treatments,replicates) %>% unlist)
    lane = factor(1:suggested.numlane)
    adapter = factor(1:suggested.nperlane)
    blocks = expand.grid(lane,adapter)
    names(blocks)=c("lane","adapter")

    if (nrow(blocks)>sum(replicates)) {
      ntoremove = nrow(blocks) - sum(replicates)
      blocks=blocks[!((blocks$lane==suggested.numlane) & (blocks$adapter %in% ((suggested.nperlane-ntoremove+1):suggested.nperlane))), ]
      # blocks=blocks %>% dplyr::filter_(!("lane"==suggested.numlane & "adapter" %in% ((suggested.nperlane-ntoremove+1):suggested.nperlane))) #not working
    }


    if (suggested.numlane <=8) {
      suggested.des=design(trtforeachsample,blocks,searches=searches,weighting=0, seed=seed)

      out$suggestedDesign$design <- suggested.des$design
      names(out$suggestedDesign$design)[names(out$suggestedDesign$design)=="TF"]="treatment"
      out$suggestedDesign$design$flowcell = 1
      out$suggestedDesign$design$lane = as.numeric(out$suggestedDesign$design$lane)
      out$suggestedDesign$design$adapter = as.numeric(out$suggestedDesign$design$adapter)
      out$suggestedDesign$design = out$suggestedDesign$design %>%
        select_("flowcell","lane","adapter","treatment") %>%
        arrange_("flowcell","lane","adapter","treatment")

      out$suggestedDesign$BlocksEfficiency <- suggested.des$blocks_model[,c(1,3)]
      names( out$suggestedDesign$BlocksEfficiency) <- c("Blocks","D-efficiency")
      out$suggestedDesign$BlocksEfficiency$Blocks <- c("lane","adapter")
    } else {
      blocks$flowcell= as.factor(ifelse(as.numeric(blocks$lane) %% 8!=0,as.numeric(blocks$lane) %/%8 +1,as.numeric(blocks$lane) %/%8))
      suggested.des=design(trtforeachsample,blocks,searches=searches,weighting=0, seed=seed)

      out$suggestedDesign$design <- suggested.des$design
      names(out$suggestedDesign$design)[names(out$suggestedDesign$design)=="TF"]="treatment"
      out$suggestedDesign$design$lane=ifelse(as.numeric(out$suggestedDesign$design$lane) %% 8!=0,as.numeric(out$suggestedDesign$design$lane) %% 8,8)
      out$suggestedDesign$design$adapter=as.numeric(out$suggestedDesign$design$adapter)
      out$suggestedDesign$design$flowcell=paste0("flowcell",as.numeric(out$suggestedDesign$design$flowcell))
      out$suggestedDesign$design = out$suggestedDesign$design %>%
        select_("flowcell","lane","adapter","treatment") %>%
        arrange_("flowcell","lane","adapter","treatment")

      out$suggestedDesign$BlocksEfficiency <- suggested.des$blocks_model[,c(1,3)]
      names( out$suggestedDesign$BlocksEfficiency) <- c("Blocks","D-efficiency")
      out$suggestedDesign$BlocksEfficiency$Blocks <- c("flowcell","lane","adapter")

    }

    if (!is.null(nperlane)) {

      # if ((sum(replicates)%%nperlane) != 0) {
      #   # not filing a whole lane, add samples belong to ' ' to fill
      #   replicates <- c(replicates, (sum(replicates)%%nperlane))
      #   treatments <- c(treatments, " ")
      # }

      numlane <- ceiling(sum(replicates)/nperlane)
      trtforeachsample = as.factor(mapply(rep,treatments,replicates) %>% unlist)
      lane = factor(1:numlane)
      adapter = factor(1:nperlane)
      blocks = expand.grid(lane,adapter)
      names(blocks)=c("lane","adapter")

      if (nrow(blocks)>sum(replicates)) {
        ntoremove = nrow(blocks) - sum(replicates)
        blocks=blocks[!((blocks$lane==numlane) & (blocks$adapter %in% ((nperlane-ntoremove+1):nperlane))), ]
        # blocks=blocks %>% dplyr::filter_(!("lane"==numlane & "adapter" %in% ((nperlane-ntoremove+1):nperlane)))
      }


      if (numlane <=8) {
        des=design(trtforeachsample,blocks,searches=searches,weighting=0, seed=seed)

        out$Design$design <- des$design
        names(out$Design$design)[names(out$Design$design)=="TF"]="treatment"
        out$Design$design$flowcell = 1
        out$Design$design$lane = as.numeric(out$Design$design$lane)
        out$Design$design$adapter = as.numeric(out$Design$design$adapter)
        out$Design$design = out$Design$design %>%
          select_("flowcell","lane","adapter","treatment") %>%
          arrange_("flowcell","lane","adapter","treatment")

        out$Design$BlocksEfficiency <- des$blocks_model[,c(1,3)]
        names( out$Design$BlocksEfficiency) <- c("Blocks","D-efficiency")
        out$Design$BlocksEfficiency$Blocks <- c("lane","adapter")
      } else {
        blocks$flowcell= as.factor(ifelse(as.numeric(blocks$lane) %% 8!=0,as.numeric(blocks$lane) %/%8 +1,as.numeric(blocks$lane) %/%8))
        des=design(trtforeachsample,blocks,searches=searches,weighting=0, seed=seed)

        out$Design$design <- des$design
        names(out$Design$design)[names(out$Design$design)=="TF"]="treatment"
        out$Design$design$lane=ifelse(as.numeric(out$Design$design$lane) %% 8!=0,as.numeric(out$Design$design$lane) %% 8,8)
        out$Design$design$adapter=as.numeric(out$Design$design$adapter)
        out$Design$design$flowcell=paste0("flowcell",as.numeric(out$Design$design$flowcell))
        out$Design$design = out$Design$design %>%
          select_("flowcell","lane","adapter","treatment") %>%
          arrange_("flowcell","lane","adapter","treatment")

        out$Design$BlocksEfficiency <- des$blocks_model[,c(1,3)]
        names( out$Design$BlocksEfficiency) <- c("Blocks","D-efficiency")
        out$Design$BlocksEfficiency$Blocks <- c("flowcell","lane","adapter")

      }

    } else {
      # nperlane is null
      out$Design <- out$suggestedDesign  #when nperlane is NULL,  design by input is the same as suggested design
    }
    out  #output
}


#' Extract design
#'
#' Extracts desgins from a \code{DEdesign} object
#'
#' @param x a \code{DEdesign} object (output from \code{gendesign}).
#' @param selection choose which design to output: \code{Design} based on input or \code{suggestedDesign}. Default is \code{Design}.
#'
#' @return a data frame giving flowcell, lane and adapter assignment for each experimental unit
#'
#'
#' @keywords RNA-seq, statistical experimental design, block design, Illumina flow cell
#'
#'
#' @examples
#' des <- gendesign(treatments=letters[1:6], replicates=rep(3,6), nperlane=4)
#' designDF(des,selection="Design")
#' designDF(des,selection="suggestedDesign")
#'
#' @export
designDF <- function(x, selection) UseMethod("designDF")

#' @export
designDF.DEdesign <- function(x, selection="Design") {
  DEdesignobj=x
  stopifnot(class(DEdesignobj) == "DEdesign")
  if (! (selection %in% c("Design","suggestedDesign"))) {
    stop("Please select which design to plot: Design based on input or suggestedDesign")
  }

  if (selection=="Design") {designtable=DEdesignobj$Design$design
  } else {designtable=DEdesignobj$suggestedDesign$design}
  designtable
}

#' Extract efficiency
#'
#' Extracts D-efficiency from a \code{DEdesign} object
#'
#' @param x a \code{DEdesign} object (output from \code{gendesign}).
#' @param selection choose which design to output: "\code{Design}" based on input or "\code{suggestedDesign}". Default is "\code{Design}".
#'
#' @return a data frame giving D-efficiencies of flowcell (If there are more than one flowcell), lane and adapter blocks
#'
#'
#' @keywords RNA-seq, statistical experimental design, block design, Illumina flow cell
#'
#'
#' @examples
#' des <- gendesign(treatments=letters[1:6], replicates=rep(3,6), nperlane=4)
#' efficiency(des,selection="Design")
#' efficiency(des,selection="suggestedDesign")
#'
#' @export
efficiency <- function(x, selection) UseMethod("efficiency")

#' @export
efficiency.DEdesign <- function(x, selection="Design") {
  DEdesignobj=x
  stopifnot(class(DEdesignobj) == "DEdesign")
  if (! (selection %in% c("Design","suggestedDesign"))) {
    stop("Please select which design to plot: Design based on input or suggestedDesign")
  }

  if (selection=="Design") {efficiencytable=DEdesignobj$Design$BlocksEfficiency
  } else {efficiencytable=DEdesignobj$suggestedDesign$BlocksEfficiency}
  efficiencytable
}

