#' Field map of RNA-seq experimental design
#'
#' Plot the RNA-seq experimental design in the format of Illumina 8-lane flow cell
#'
#'
#' @param designtable a data frame for lane and adaptor assignment. Use output from \code{gen.design}.It should have 4 columns with column names \code{chip, lane, adp, trt}.
#'
#' @return None
#'
#'
#' @keywords RNA-seq, statistical experimental design, block design, Illumina flow cell
#'
#' @export
#'
#' @examples
#' des <- gendesign(trts=letters[1:6], reps=rep(3,6), nperlane=4)
#' plotdesign (des$Design$design)
#' plotdesign (des$suggestedDesign$design)
#'
#' @import dplyr
#' @import ggplot2
#' @importFrom RColorBrewer brewer.pal
#' @importFrom grDevices colorRampPalette
#'
#'
#'
#'


plotdesign <- function(designtable) {
  stopifnot(class(designtable) == "data.frame")
  inputs <- names(designtable)
  if (!("lane" %in% inputs)) {
    stop("need lane assignment in design")
  }
  if (sum(designtable$lane %in% 1:8) != length(designtable$lane)) {
    stop("lane needs to be 1 to 8")
  }
  if (!("adp" %in% inputs)) {
    stop("need adaptor assignment in design")
  }
  if (!("chip" %in% inputs)) {
    stop("need chip assignment in design")
  }

  nadp <- max(designtable$adp)
  lanelabel <- data.frame(lane = 1:8, adp = nadp + 1.2, label = paste0("lane", 1:8))
  adplabel <- data.frame(lane = rep(0, nadp), adp = (1:nadp) + 0.5,
                         label = paste(paste0("sample", 1:nadp),paste0("adp", 1:nadp),sep=":"))

  ntrt <- unique(designtable$trt) %>% length

  if (ntrt <= 12) {
    if (" " %in% designtable$trt) {
      colors <- c("grey100", brewer.pal(ntrt - 1, "Paired"))
    } else {
      colors <- brewer.pal(ntrt, "Paired")
    }
  } else {
    if (" " %in% designtable$trt) {
      colors <- c("grey100", colorRampPalette(rev(brewer.pal(n = 7, name = "RdYlBu")))(ntrt - 1))
    } else {
      colors <- colorRampPalette(rev(brewer.pal(n = 7, name = "RdYlBu")))(ntrt)
    }
  }

  p <- ggplot(designtable, aes(xmin = lane - 0.5, xmax = lane + 0.5, ymin = adp, ymax = adp + 1)) +
    xlim(-0.5, 8.5) + ylim(0, nadp + 2) +
    geom_rect(aes(fill = trt), colour = "grey50") +
    scale_fill_manual(values = colors) +
    geom_text(aes(x = lane, y = adp + 0.5, label = trt), show.legend = FALSE) +
    geom_text(data = lanelabel,
              aes(x = lane, y = adp, label = label)) +
    geom_text(data = adplabel, aes(x = lane, y = adp, label = label)) +
    theme_bw() +
    theme(panel.border = element_blank(), panel.grid.major.x = element_blank(), panel.grid.minor.x = element_blank(),
          panel.grid.major.y = element_blank(), panel.grid.minor.y = element_blank(), legend.position = "none",
          axis.title.x = element_blank(), axis.text.x = element_blank(), axis.ticks.x = element_blank(),
          axis.title.y = element_blank(), axis.text.y = element_blank(), axis.ticks.y = element_blank())

  if (length(unique(designtable$chip)) == 1) {
    p
  } else {
    p + facet_grid(chip ~ .)
  }
}
