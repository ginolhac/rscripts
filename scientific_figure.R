library("grid")
library("gridBase")


# Draws 1 or more figures in a format suitable for submission, in which each
# subfigure (if more than one) is labeled as a, b, c, d, and a figure title is
# added if not set to NULL. The height of the title is set to 1 by convention,
# and heights of other rows are set relative to this.
#
# Example:
#   pdf(sprintf("figure_1.pdf"), width=14, height=6)
#   sci.plot(function() plot(rnorm(100)),
#            ggplot(...) + ...,
#            title="Figure 1",
#            widths=c(5, 5),
#            ncol=2)
#   dev.off()
#
sci.plot <- function(..., nrow=1, ncol=1, widths=5, heights=5,
                     numbering=letters, title="Figure 1",
                     xadj=0, yadj=0)
{
    stopifnot(nrow >= 1)
    stopifnot(ncol >= 1)
    stopifnot(length(widths) == ncol)
    stopifnot(length(heights) == nrow)
    stopifnot(length(numbering) >= (nrow * ncol))

    plots <- list(...)
    stopifnot(length(plots) == (nrow * ncol))

    extra.row = ifelse(is.null(title), 0, 1)

    if (!is.null(title)) {
        heights <- c(1, heights)
    }

    plot.new()
    layout <- grid.layout(nrow=nrow + extra.row, ncol=ncol,
                          widths=widths,
                          heights=heights)

    if (!is.null(title)) {
        vp <- viewport(layout.pos.col=1, layout.pos.row=1)
        pushViewport(vp)
        par(new=TRUE, fig=gridFIG())
        grid.text(title,
                  hjust=0, vjust=2,
                  x=unit(0.00, "npc"),
                  y=unit(1.00, "npc"),
                  gp=gpar(fontsize=30))
        popViewport()
    }

    pushViewport(viewport(layout=layout))
    for (row in 1:nrow) {
        for (col in 1:ncol) {
            vp <- viewport(layout.pos.col=col, layout.pos.row=(row + extra.row))
            pushViewport(vp)

            idx <- nrow * (row - 1) + col
            func <- plots[[idx]]
            par(new=TRUE, fig=gridFIG())

            if(is.function(func)) {
                pp <- func()
            } else {
                pp <- func
            }

            if (!is.null(pp)) {
                print(pp, newpage=FALSE)
            }

            if (length(plots) > 1) {
                grid.text(numbering[[idx]],
                          hjust=0, vjust=2,
                          x=unit(0.00 + xadj, "npc"),
                          y=unit(1.02 + yadj, "npc"),
                          gp=gpar(fontsize=20))
            }

            popViewport()
        }
    }

    popViewport(1)
}
