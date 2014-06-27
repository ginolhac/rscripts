#!/usr/bin/Rscript
# Wrapper around the 'ape' package, to simplify drawing of phylogenetic trees

library(ape)

# Enable full backtraces on errors
on_error <- function(e)
{
  traceback()
  quit(status = 1)
}
options(error = on_error)


get.clade.node <- function(tr, taxa)
{
    # Get node with lowest depth
    nodes <- tr$edge[which.edge(tr, taxa)]
    depths <- node.depth(tr)[nodes]
    node <- nodes[order(depths, decreasing=TRUE)][1]
    return(node)
}


get.leaf.nodes <- function(tr, leafs)
{
    nodes <- NULL
    for (leaf in leafs) {
        nodes <- append(nodes, get.clade.node(tr, leaf))
    }
    return(nodes)
}


set.tip.label <- function(tr, old, new, color=NULL)
{
  if (is.null(tr$verbose.tip.label)) {
    tr$verbose.tip.label <- tr$tip.label
  }

  if (!is.null(color)) {
    if (is.null(tr$verbose.tip.color)) {
        tr$verbose.tip.color <- rep("black", length(tr$tip.label))
    }

    tr$verbose.tip.color[tr$verbose.tip.label == old] <- color
  }

  tr$verbose.tip.label[tr$verbose.tip.label == old] <- new
  return(tr)
}


set.node.label <- function(tr, taxa, new)
{
  if (is.null(tr$verbose.node.label)) {
    tr$verbose.node.label <- tr$node.label
  }

  node <- get.clade.node(tr, taxa) - length(tr$tip.label)
  tr$verbose.node.label[node] <- new
  return(tr)
}


set.clade.defaults <- function(tr,
                               edge.color="black",
                               edge.width=2)
{
  if (is.null(tr$clade.edge.color)) {
    tr$clade.edge.color <- rep(edge.color, length(tr$edge))
  }

  if (is.null(tr$clade.edge.width)) {
    tr$clade.edge.width <- rep(edge.width, length(tr$edge))
  }

  return(tr)
}


set.clade.settings <- function(tr, taxa,
                               edge.color="black",
                               edge.width=2,
                               mark.clade=TRUE)
{
  tr <- set.clade.defaults(tr)

  for (ntaxa in 1:length(taxa)) {
    for (selection in combn(taxa, ntaxa, simplify=FALSE)) {
      edges <- which.edge(tr, selection)
      tr$clade.edge.width[edges] <- edge.width
      tr$clade.edge.color[edges] <- edge.color
    }
  }

  if (mark.clade) {
    node <- get.clade.node(tr, taxa)
    tr$clade.mark <- append(tr$clade.mark, node)
    tr$clade.mark.color <- append(tr$clade.mark.color, edge.color)
  }

  return(tr)
}


do.plot.tree <- function(tr, ...)
{
  if (is.null(tr$clade.edge.color)) {
    tr$clade.edge.color <- "black"
  }

  if (is.null(tr$clade.edge.width)) {
    tr$clade.edge.width <- 2
  }

  if (!is.null(tr$verbose.tip.label)) {
      tr$tip.label <- tr$verbose.tip.label
  }

  plot(tr, no.margin=TRUE,
       label.offset=0.00005,
       edge.width=tr$clade.edge.width,
       edge.color=tr$clade.edge.color,
       tip.color=tr$verbose.tip.color,
       ...)
}

do.plot.nodelabels <- function(tr, nodes = NULL, exclude = NULL, underscore=FALSE,
                               cex=0.75, frame="none", col="darkgrey",
                               ...)
{
  if (is.null(nodes)) {
      nodes <- unique(tr$edge[,1])
  }

  if (is.null(tr$verbose.node.label)) {
      text <- tr$node.label
  } else {
      text <- tr$verbose.node.label
  }

  if (!is.null(text)) {
    nodes <- setdiff(nodes, exclude)
    text <- text[nodes - length(tr$tip.label)]
    if (!underscore) {
        text <- gsub("_", " ", text)
    }

    nodelabels(text, nodes, cex=cex, frame=frame, col=col, ...)
  }
}


do.plot.clademarks <- function(tr, ...)
{
  if (!is.null(tr$clade.mark)) {
    nodelabels(rep(" ", length(tr$clade.mark)), cex=0.5,
               frame="circle",
               node=tr$clade.mark,
               col=tr$clade.mark.color,
               bg=tr$clade.mark.color,
               ...)
  }
}


do.add.scale.bar <- function(...)
{
    add.scale.bar(...)
}


if (!exists("standalone")) {
  args <- commandArgs(trailingOnly=TRUE)
  if (length(args) == 2)
  {
    tr <- read.tree(args[1])
    template <- file(paste(args[2], "R", sep="."), "w")
    cat("standalone <- NULL\n", file=template)
    cat("source(\"", sub(".*=", "", commandArgs()[4]), "\")\n", sep="", file=template)
    cat("# Update tip labels\n", file=template)
    cat("tr <- read.tree(\"", args[1], "\")\n\n", sep="", file=template)

    cat("tr <- set.clade.defaults(tr, edge.color=\"black\", edge.width=2)\n", sep="", file=template)

    cat("\n# Update tip labels\n", file=template)
    for (name in tr$tip.label) {
      cat("# tr <- set.tip.label(tr, \"", name, "\", \"", name, "\")\n", sep="", file=template)
    }

    cat("\n# Mark clades (specified as list of taxa names)\n", file=template)
    cat("# tr <- set.clade.settings(tr, c(\"name_1\", \"name_2\", ...), edge.width=1, edge.color=\"black\", mark.clade=TRUE)\n",
        file=template)

    cat("\n# Draw phylogeny\n", file=template)
    cat("pdf(\"", paste(args[2], "pdf", sep="."), "\", height=3.5)\n",
        sep="", file=template)

    cat("do.plot.tree(tr)\n", file=template)
    cat("do.plot.nodelabels(tr, adj=c(1.25, 1.5))\n", file=template)
    cat("do.plot.clademarks(tr)\n", file=template)
    cat("add.scale.bar()\n", file=template)
    cat("dev.off()\n", file=template)
    close(template)

    cat("Template written!\n")
    cat("Update and execute", paste(args[2], "R", sep="."), "to plot phylogeny:\n")
    cat("$ Rscript ", paste(args[2], "R", sep="."), "\n", sep="")
  } else {
    cat("Usage: Rscript plot_tree.R <tree.newick> <output>\n")
    quit(status=1)
  }
}
