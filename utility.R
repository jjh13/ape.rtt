library(ape)

read.tree.tips <- function(file) {
	tree <- read.tree(file)
	tree$tip.dates <- as.numeric(gsub("(.+)_([0-9]+)$", "\\2", tree$tip.label, perl=T))
	return (tree)
}

plot.regression <- function(tree) {
	distances <- node.depth.edgelength(tree)[1:(t$Nnode + 1)]
	dates <- tree$tip.dates

	m <- lm(distances ~ dates)
	plot(distances, dates)
	abline(m)
}