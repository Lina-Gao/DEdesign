#' Field map of RNA-seq experimental design
#'
#' Plot the RNA-seq experimental design in the format of Illumina 8-lane flow cell
#'
#'
#' @param designtable a data frame for lane and adaptor assignment. Use output from \code{gen.design}.It should have 4 columns with column names \code{flowcell, lane, adaptor, treatment}.
#'
#' @return None
#'
#'
#' @keywords RNA-seq, statistical experimental design, block design, Illumina flow cell
#'
#' @export
#'
#' @examples
#' des <- gendesign(treatments=letters[1:6], replicates=rep(3,6), nperlane=4)
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
  if (!("adaptor" %in% inputs)) {
    stop("need adaptor assignment in design")
  }
  if (!("flowcell" %in% inputs)) {
    stop("need flowcell assignment in design")
  }

  nadaptor <- max(designtable$adaptor)
  lanelabel <- data.frame(lane = 1:8, adaptor = nadaptor + 1.2, label = paste0("lane", 1:8))
  adaptorlabel <- data.frame(lane = rep(0, nadaptor), adaptor = (1:nadaptor) + 0.5,
                         label = paste0("adaptor", 1:nadaptor))

  ntreatment <- unique(designtable$treatment) %>% length

  if (ntreatment <= 12) {
    if (" " %in% designtable$treatment) {
      colors <- c("grey100", brewer.pal(ntreatment - 1, "Paired"))
    } else {
      colors <- brewer.pal(ntreatment, "Paired")
    }
  } else {
    if (" " %in% designtable$treatment) {
      colors <- c("grey100", colorRampPalette(rev(brewer.pal(n = 7, name = "RdYlBu")))(ntreatment - 1))
    } else {
      colors <- colorRampPalette(rev(brewer.pal(n = 7, name = "RdYlBu")))(ntreatment)
    }
  }

  p <- ggplot(designtable, aes_(xmin = ~lane - 0.5, xmax = ~lane + 0.5, ymin = ~adaptor, ymax = ~adaptor + 1)) +
    xlim(-0.5, 8.5) + ylim(0, nadaptor + 2) +
    geom_rect(aes_(fill = ~treatment), colour = "grey50") +
    scale_fill_manual(values = colors) +
    geom_text(aes_(x = ~lane, y = ~adaptor + 0.5, label = ~treatment), show.legend = FALSE) +
    geom_text(data = lanelabel,
              aes_(x = ~lane, y = ~adaptor, label = ~label)) +
    geom_text(data = adaptorlabel, aes_(x = ~lane, y = ~adaptor, label = ~label)) +
    theme_bw() +
    theme(panel.border = element_blank(), panel.grid.major.x = element_blank(), panel.grid.minor.x = element_blank(),
          panel.grid.major.y = element_blank(), panel.grid.minor.y = element_blank(), legend.position = "none",
          axis.title.x = element_blank(), axis.text.x = element_blank(), axis.ticks.x = element_blank(),
          axis.title.y = element_blank(), axis.text.y = element_blank(), axis.ticks.y = element_blank())

  if (length(unique(designtable$flowcell)) == 1) {
    p
  } else {
    p + facet_grid(flowcell ~ .)
  }
}
