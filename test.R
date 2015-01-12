test.rtt <- function(t, tip.dates, remove = 1, random = T) {
	t<-unroot(t)
	to.remove <- if(!random) 
	remove
	else
	sample(length(tip.dates), remove, replace=F)
	mdates <- tip.dates
	mdates[to.remove] <- (NA)
	t <- rtt.missing.dates(t, mdates)
	sum((t$tip.dates - tip.dates)^2)/length(tip.dates)
}