library(ape)

read.tree.tips <- function(file) {
	tree <- read.tree(file)
	tree$tip.dates <- as.numeric(gsub("(.+)_([0-9]+)$", "\\2", tree$tip.label, perl=T))
	return (tree)
}

plot.regression <- function(tree, tip.dates=NULL) {
	distances <- node.depth.edgelength(tree)[1:(tree$Nnode + 1)]

	dates <- tree$tip.dates[1:(tree$Nnode + 1)]
	if(!is.null(tip.dates)) {
		dates <- tip.dates[1:(tree$Nnode + 1)]
	}

	m <- lm(distances ~ dates)
	plot(distances, dates)
	abline(m)
}