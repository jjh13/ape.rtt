# Generates the guide trees for the simulated
# molecular data. Once this is done, an  
# indelible control file is output

library(ape)
library(TreeSim)

#source('../rtt.R')
#source('../test.R')

# Number of guide trees to create
n.trees <- 50
n.partitions <- 50
n.replicates <- 2
n.tips <- 50


# These are somewhat abitrarily set
lambda <- c(2)
mu <- c(2)
sampprob <-c(0.5)
times<-c(0)

date_branches <- function(tree) {
	tree <- unroot(tree[[1]])

	tip.dates <- node.depth.edgelength(tree)[1:n.tips]
#	tree$edge.length <- tree$edge.length * 100
	dates <- unlist(Map(toString, node.depth.edgelength(tree)[1:n.tips]))

	# Quick hack for now, I know I can do this in one or two lines
	for(i in 1:n.tips)
		tree$tip.label[i] <- paste(tree$tip.label[i], dates[i], sep='_')

	print(tree)
	tree #rtt(tree, tip.dates)
}

trees <- apply(matrix(rep(n.tips,n.trees)), 1, sim.bdsky.stt, lambdasky=lambda, deathsky=mu, timesky=times, sampprobsky=sampprob,rho=0,timestop=0)
#trees <- apply(matrix(rep(n.tips,n.trees)), 1, rtree)
for(i in 1:n.trees){
	trees[[i]] <- date_branches(trees[[i]])

}

indel_control <- 
"
[TYPE] NUCLEOTIDE 2

[SETTINGS]
  [output]                   FASTA 

[MODEL]    HKY_HIV
  [submodel]  HKY 9.5               
  [statefreq] 0.42 0.15 0.15 0.28

"

for(i in 1:n.trees) {
	tree_dat <- write.tree(trees[[i]])
	#tree_dat <- substr(tree_dat, 1, nchar(tree_dat))
	print(tree_dat)
	indel_control <- paste0(indel_control, sprintf("[TREE] tree_%d %s \n", i, tree_dat))
}

index <- sample(1:n.trees, n.partitions, replace = (n.partitions > n.trees))
for(i in 1:n.partitions) {
	indel_control <- paste0(indel_control, sprintf("[PARTITIONS] pHKY_%d [tree_%d HKY_HIV 900] \n", i, index[i]))
}

indel_control <- paste0(indel_control, "[EVOLVE] \n")
for(i in 1:n.partitions) {
	indel_control <- paste0(indel_control, sprintf("    pHKY_%d %d HIV_%d_out \n", i, n.replicates, i))
}

write(indel_control, 'control.txt')
