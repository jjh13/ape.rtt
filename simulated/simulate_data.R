# Generates all the molecular 
# data for the simulated results

library(ape)
library(TreeSim)
library(phangorn)
library(NELSI)

date.branches <- function(s.tree) {
	tree <- s.tree[[1]]
	subs <- s.tree[[2]][, 6]
	times <- s.tree[[2]][, 7]

	print(cbind(subs, times))

	dates <- unlist(Map(toString, times))[1:n.tips]
	tree$tip.label <- paste(tree$tip.label, dates, sep='_')
	tree
}

write.sequences <- function(seq.n){
	seq <- seq.n[[1]]
	i <- seq.n[[2]]
	write.dna(as.DNAbin(seq), sprintf("HIV_SIM_%d.phy", i))
}

# set the seed for this run
set.seed(19902302)

# Number of guide trees to create
n.trees <- 2
n.tips <- 50
n.seqlen <- 630

# Parameters for the nucleotide simulation
# HKY85 with rate = kappa
kappa <- 6.0
Q <- c(kappa, 1, 1, 1, 1, kappa)
pi.vec <- c(0.42, 0.15, 0.15, 0.28)
seq.rate <- 2.0

# Parameters for initial guide tree(s)
# These are somewhat abitrarily set
lambda <- c(2)
mu <- c(2)
sampprob <-c(0.5)
times<-c(0)

# Parameters for the clock of the tree
sim.params <- list(rate = 1.0, noise = 0)
sim.clockmodel <- simulate.clock

#trees <- apply(matrix(rep(n.tips,n.trees)), 1, sim.bdsky.stt, lambdasky=lambda, deathsky=mu, timesky=times, sampprobsky=sampprob,rho=0,timestop=0)
#trees <- lapply(trees, function(x) {x[[1]]}) # undo the strange indexing that sim.bdsky forces on you

trees <- apply(matrix(rep(n.tips,n.trees)), 1, rtree)

sim.trees <- lapply(trees, sim.clockmodel, params=sim.params)
trees <- lapply(sim.trees, date.branches)
seqs <- lapply(trees, simSeq, l=n.seqlen, Q=Q, bf=pi.vec, rate=seq.rate)

apply(cbind(seqs, 1:n.trees), 1, write.sequences)