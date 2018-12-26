#' Statistical design for RNA-seq experiments aimed at identifying differentially expressed genes
#'
#' Generates block design to assign samples to sequencing lanes and adapters for Illumina flow cell
#'
#'
#' @param treatments a \code{data.frame} showing number of samples per treatment group. For single factor experiment, \code{treatments} should contain two columns with one column showing levels of the factor and the other column named \code{"replicates"} for number of samples in each level. For factorial experiments, there should be one column for each factor, and the \code{"replicates"} column gives number of samples in each combination of factor levels. Please see \code{treatments.example} for a template. For the purpose of plotting, please use short factor names and in factorial design using numbers for factor levels is recommended.
#' @param nperlane An integer for number of samples per lane for sequencing. Default is 4. This is used to generate \code{Design} in results.
#' @param search.surrounding An non-negative integer with default 0, which means to only find design for number of samples per lane defined by \code{nperlane}. For other values, designs using \code{min(3,(nperlane-search.surrounding))} to \code{(nperlane+search.surrounding)} samples per lane will be compared to give \code{suggestedDesign}. \code{suggestedDesign} is the design with minimal number of flowcells among all candidate designs that gives the highest lane block efficiency (when more than one designs meet these criteria, design with highest nperlane is selected).
#' @param seed an integer initializing the random number generator. The default is \code{seed=1}. Designs can be rebuilt repeatedly using different seed to check that a near-optimum design has been found.
#' @param searches the maximum number of local optima searched at each stage of a treatment and block design optimization. The default depends on the design size. For optimum results, try large number of searches.
#'
#' @return \code{gendesign} returns an object of class "\code{DEdesign}", which is a list containing the following components:
#'
#' @return \code{input} A list showing input parameters to the function
#' @return \code{Design} A list with two elements: \code{design} and \code{BlocksEfficiency}. \code{design} is a data frame giving flowcell, lane and adapter assignment for treatment groups. \code{BlocksEfficiency} is a data frame giving block efficiencies (D-Efficiencies) for lane and adapter (and/or flowcell when applicable).
#' @return \code{suggestedDesign} When \code{search.surrounding >0}, a list with two elements are returned here: \code{design} and \code{BlocksEfficiency}. \code{design} is a data frame giving flowcell, lane and adapter assignment for treatment groups. \code{BlocksEfficiency} is a data frame giving block efficiencies (D-Efficiencies) for lane and adapter (and/or flowcell when applicable). Designs using \code{min(3,(nperlane-search.surrounding))} to \code{(nperlane+search.surrounding)} samples per lane will be compared to give \code{suggestedDesign}. \code{suggestedDesign} is the design with minimal number of flowcells among all candidate designs that gives the highest lane block efficiency (when more than one designs meet these criteria, design with highest nperlane is selected).
#'
#'
#' @keywords RNA-seq, statistical experimental design, block design, Illumina flow cell
#'
#'
#' @exportClass DEdesign
#'
#' @examples
#' gendesign(treatments= data.frame(trt = letters[1:4], replicates = rep(4,4)))
#' treatments <- data.frame(expand.grid(A=factor(1:2), B=factor(1:5)),replicates = 2)
#' gendesign(treatments=treatments, nperlane=4, search.surrounding = 2)
#'
#' @import dplyr
#' @importFrom purrr map
#' @importFrom blocksdesign design
#'
#' @export
#'
gendesign <- function(treatments, nperlane = 4, search.surrounding = 0, seed = 1, searches = NULL) {
  if (class(treatments) != "data.frame") {
    stop("treatments needs to be a data.frame showing number of samples per treatment group")
  }


  if (!((class(nperlane) == "integer" & nperlane > 0) | (class(nperlane) == "numeric" &
                                                         nperlane > 0 & ceiling(nperlane) == floor(nperlane)))) {
    stop("nperlane needs to be a positive integer")
  }

  if (!((class(search.surrounding) == "integer" & search.surrounding >= 0) | (class(search.surrounding) == "numeric" &
                                                                              search.surrounding >= 0 & ceiling(search.surrounding) == floor(search.surrounding)))) {
    stop("search.surrounding needs to be a non-negative integer")
  }

  if (!("replicates" %in% names(treatments))) {
    stop("input `treatments` data.frame has to contain a column named 'replicates' showing number of samples per treatment group")
  } else {
    if (!((class(treatments$replicates) == "integer" & all(treatments$replicates > 0)) | (class(treatments$replicates) == "numeric" &
                                                                                          all(treatments$replicates > 0) &
                                                                                          all(ceiling(treatments$replicates) == floor(treatments$replicates))))) {
      stop("`replicates` column in input `treatments` data.frame needs to be positive integers")
    }
  }


  out <- structure(list(input = list(treatments = treatments, nperlane = nperlane, search.surrounding = search.surrounding, seed = seed, searches = searches),
                        Design = list(),
                        suggestedDesign = list()),
                   class = "DEdesign")


  treatments = .reformat(treatments)

  nperlane_torun = max(3, (nperlane-search.surrounding)):(nperlane+search.surrounding)

  cachedRes = data.frame(nperlane_torun = nperlane_torun)
  cachedRes$res_list = map(cachedRes$nperlane_torun, .f = ~.onedesign (nperlane=.,treatments = treatments, seed=seed, searches=searches))

  out$Design = cachedRes$res_list[[which(cachedRes$nperlane_torun==nperlane)]]

  if(search.surrounding > 0) {
    # select nperlane that gives highest D-efficiency for lane (note lane is nested within flowcell, but it is not considered here, will improve in the future)
    out$suggestedDesign = cachedRes$res_list[[which(cachedRes$nperlane_torun==.whichHighestEffi(cachedRes, target = "lane"))]]
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
#' treatments <- data.frame(trt = letters[1:6], replicates = rep(3,6))
#' des <- gendesign(treatments=treatments, nperlane=4, search.surrounding = 2)
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
#' treatments <- data.frame(trt = letters[1:6], replicates = rep(3,6))
#' des <- gendesign(treatments=treatments, nperlane=4, search.surrounding = 2)
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

