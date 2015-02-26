library(ape)


# Basically, this function builds a distribution of
# possible dates for a tip by back propogating 
# times from 
date.tips <- function(t, tip.dates, edge.clocks, criterion='MEAN') {

	# Find the tips of the tree that're missing dates
	missing.indices <- which(is.na(tip.dates))
	valid.indices <- which(!is.na(tip.dates))

	# Setup a new tree with edge lenths as units of
	# time scaled by the local clock over each edge.
	t.delta <- t
	t.delta$edge.length <- t.delta$edge.length / edge.clocks

	# This gives the time from 'root' to each tip
	root.to.tip.time <- node.depth.edgelength(t.delta)

	# This is the offset each tip is 'wrong' by
	delta.t <- tip.dates - root.to.tip.time[valid.indices] 

	t.dist <- t(matrix(rep(delta.t, length(missing.indices)), nrow=length(valid.indices)))
	t.dist <- t.dist + matrix(rep(root.to.tip.time[missing.indices], lenths(valid.indices)), nrow=length(missing.indices))

	tip.dates[missing.indices] <- apply(t.dist, 1, mean)
	tip.dates
}