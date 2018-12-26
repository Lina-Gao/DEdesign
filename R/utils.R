#' reformat input treatments data.frame
#'
#' @import dplyr
#' @importFrom purrr map map_dbl
#' @importFrom blocksdesign design
#'
#' @param treatments a \code{data.frame} showing number of samples per treatment group. For single factor experiment, \code{treatments} should contain two columns with one column showing levels of the factor and the other column named \code{"replicates"} for number of samples in each level. For factorial experiments, there should be one column for each factor, and the \code{"replicates"} column gives number of samples in each combination of factor levels. Please see \code{treatments.example} for a template. For the purpose of plotting, please use short factor names and in factorial design using numbers for factor levels is recommended.
#' @return a data frame of reformated input treatments dataframe to use as blocksdesign input, repeat each row "replicates" times, so nrow = total number of samples
.reformat = function (treatments){
  treatments[rep(seq_len(nrow(treatments)), times=treatments$replicates),names(treatments)!="replicates", drop = FALSE]
}

#' Calculate design for a specified nperlane
#'
#' @import dplyr
#' @importFrom blocksdesign design
#'
#' @param treatments reformated input treatments dataframe
#' @param nperlane An integer for number of samples per lane for sequencing. Default is 4. This is used to generate \code{Design} in results.
#' @param seed an integer initializing the random number generator. The default is \code{seed=1}. Designs can be rebuilt repeatedly using different seed to check that a near-optimum design has been found.
#' @param searches the maximum number of local optima searched at each stage of a treatment and block design optimization. The default depends on the design size. For optimum results, try large number of searches.
#' @return blocksdesign output for one nperlane value
.onedesign = function (treatments, nperlane, seed, searches) {
  res = list()
  numlane <- ceiling(nrow(treatments)/nperlane)
  lane = factor(1:numlane)
  adapter = factor(1:nperlane)
  blocks = expand.grid(lane,adapter)
  names(blocks)=c("lane","adapter")

  if (nrow(blocks)>nrow(treatments)) {
    ntoremove = nrow(blocks) - nrow(treatments)
    blocks=blocks[!((blocks$lane==numlane) & (blocks$adapter %in% ((nperlane-ntoremove+1):nperlane))), ]
    # blocks=blocks %>% dplyr::filter_(!("lane"==numlane & "adapter" %in% ((nperlane-ntoremove+1):nperlane)))
  }


  if (numlane <=8) {
    des=design(treatments,blocks,searches=searches,weighting=0, seed=seed)

    res$design <- des$design
    res$design$flowcell = 1
    res$design$lane = as.numeric(res$design$lane)
    res$design$adapter = as.numeric(res$design$adapter)
    # res$design = res$design %>%
    #   select_("flowcell","lane","adapter",
    #           names(res$design)[!names(res$design) %in% c("flowcell","lane","adapter")]) %>% #note only first element will be used
    #   arrange_("flowcell","lane","adapter")

    res$design = (res$design[,c("flowcell","lane","adapter",
                                names(res$design)[!names(res$design) %in% c("flowcell","lane","adapter")])]) %>%
      arrange_("flowcell","lane","adapter")

    res$BlocksEfficiency <- des$blocks_model[,c(1,3)]
    names( res$BlocksEfficiency) <- c("Blocks","D-efficiency")
    res$BlocksEfficiency$Blocks <- c("lane","adapter")
  } else {
    blocks$flowcell= as.factor(ifelse(as.numeric(blocks$lane) %% 8!=0,as.numeric(blocks$lane) %/%8 +1,as.numeric(blocks$lane) %/%8))
    des=design(treatments,blocks,searches=searches,weighting=0, seed=seed)

    res$design <- des$design
    res$design$lane=ifelse(as.numeric(res$design$lane) %% 8!=0,as.numeric(res$design$lane) %% 8,8)
    res$design$adapter=as.numeric(res$design$adapter)
    res$design$flowcell=paste0("flowcell",as.numeric(res$design$flowcell))
    res$design = (res$design[,c("flowcell","lane","adapter",
                                names(res$design)[!names(res$design) %in% c("flowcell","lane","adapter")])]) %>%
      arrange_("flowcell","lane","adapter")

    res$BlocksEfficiency <- des$blocks_model[,c(1,3)]
    names( res$BlocksEfficiency) <- c("Blocks","D-efficiency")
    res$BlocksEfficiency$Blocks <- c("flowcell","lane","adapter")

  }
  res
}

#' Find nperlane for suggestedDesign
#'
#' @importFrom purrr map map_dbl pluck
#' @param cachedRes a list data frame with columns "nperlane_torun" and "res_list" (design output based on nperlane given)
#' @param target which block efficiency should be used to select suggestedDesign, default is "lane"
#' @return \code{nperlane} value that gives the suggestedDesign

.whichHighestEffi = function(cachedRes, target = "lane"){

  cachedRes$effiDF = map(.x = cachedRes$res_list, .f = ~pluck( ., "BlocksEfficiency"))
  cachedRes$numflowcell = sapply(cachedRes$effiDF,nrow)

  cachedRes$min_numflowcell = (cachedRes$numflowcell==min(cachedRes$numflowcell))

  subuse = cachedRes[cachedRes$min_numflowcell,c("nperlane_torun","effiDF")]

  subuse$target_effi = map_dbl(subuse$effiDF, function(df){
    ((df$`D-efficiency`) %>% as.character() %>% as.numeric())[df$Blocks==target]}) #results from blocksdesign saved the efficiency as characters
  subuse$candidates = (subuse$target_effi==max(subuse$target_effi))
  max(subuse$nperlane_torun[subuse$candidates])
  # among min num_flowcell, returns nperlane_torun with the highest target
  # efficiency, when several nperlane gives the max, choose the largest nperlane
  # to use less lanes
}






