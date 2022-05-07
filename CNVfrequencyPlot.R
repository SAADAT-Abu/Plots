#' circos plot of CNVs.
#'
#' The function plots frequency of CNVs.
#'
#' @md
#' @param data CNV data should be a list of dataframe with columns for chr, start, end and freq. first dataframe in list should be for deletion second for duplication
#' @export
#' @return None.
#' @examples
#' circos.cnvplot(cnv)
#' Author: Abu Saadat

circos.cnvplot <- function(data)
{
  for(p in c("circlize","ComplexHeatmap")) {
    if (length(grep(paste("^package:", p, "$", sep=""), search())) == 0) {
      if (!requireNamespace(p, quietly = TRUE))
        warning(paste("circos.cnvplot needs package `", p, "' to be fully functional; please install", sep=""))
    }
  }
  requireNamespace("circlize")
  circlize::circos.par(start.degree = 50, track.height = 0.3, cell.padding = c(0, 0, 0, 0))
  circlize::circos.initializeWithIdeogram(species="hg19", track.height = 0.05, ideogram.height = 0.06)
  text(0, 0, "CNV Frequency", cex = 2.0)
  circlize::circos.genomicTrackPlotRegion(cnv.data, ylim = c(0,1),numeric.column = c(4,4), panel.fun = function(region,value,...) {
    i<-getI(...)
    with(cbind(region,value),circlize::circos.segments(start,freq,end,freq,col=i+1,lwd=3))
  })
  lgd_lines = Legend(at = c("Deletion", "Dupication"), type = "lines", grid_height = unit(15, "mm"),
                     grid_width = unit(15, "mm"),
                     legend_gp = gpar(col = c(2,3), lwd = 3), title_position = "topleft", 
                     title_gp = gpar(fontsize = 18, fontface = "bold"), labels_gp = gpar(fontsize = 15),
                     title = "CNV: ")
  draw(lgd_lines, x = unit(90, "mm"), y = unit(60, "mm"), just = c("left", "bottom"))
  circlize::circos.clear()
}
