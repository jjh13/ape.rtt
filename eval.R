library(ape)
library(phyloch)

#source('ape.patches.R')
source('rtt.R')
source('test.R')

trim <- function (x) gsub("^\\s+|\\s+$", "", x)
extract_dates <- function(x) as.numeric(gsub("(.+)_([0-9\\.]+)$", "\\2", x, perl=T))

n.simulated <- 50
mean_rmse <- rep(0, n.simulated)
mean_rmse5 <- rep(0, n.simulated)
mean_rmse10 <- rep(0, n.simulated)
nfp <- c(1, 2, 2, 3, 2, 3, 3, 4, 4, 5, 5, 6, 5, 6, 6, 7, 5, 6, 6, 7, 6, 7, 7, 8, 9, 10, 10, 11)

for(i in 1:n.simulated) {
	
	unlink("__tmp.dna_phyml_tree.txt", force = T)
	unlink("__tmp.dna_phyml_stats.txt", force = T)

	simenv.dna <- read.FASTA(sprintf("simulated/HIV_%d_out_TRUE.fas", i))
	uid <- sample(1:length(names(simenv.dna)), length(names(simenv.dna)), replace=F)
	# read.FASTA is buggy, and creates the DNABin structure incorrectly
	# BUT, it's so broken that we can mess around with the names
	for(j in 1:length(names(simenv.dna))) { 
		# give each bin a unique name why not,
		names(simenv.dna)[j] <- sprintf("u%s_%s", uid[j], names(simenv.dna)[j])
	}

	# then output this file
	write.dna(simenv.dna, "__tmp.dna")
	# and re-read it as a proper DNA Bin
	simenv.dna <- read.dna("__tmp.dna")

	phyml.env <- phymltest("__tmp.dna", execname = "~/Binaries/phyml")
	aic <- 2 * (nfp - phyml.env)
	plot(phyml.env)
	trees <- read.tree("__tmp.dna_phyml_tree.txt")

	# Choose the tree with the min aic
	tree <- trees[[which(aic==min(aic))]]

	plot(tree)
	write.tree(tree, sprintf("simulated/tree/HIV_ml_%d_out.nwk", i))

	tip.dates <- trim(tree$tip.label)
	tip.dates <- extract_dates(tree$tip.label)

	mean_rmse[i] <- test.rtt(tree, tip.dates, remove = 1, random = T)
	mean_rmse5[i] <- test.rtt(tree, tip.dates, remove = 5, random = T)
	mean_rmse10[i] <- test.rtt(tree, tip.dates, remove = 10, random = T)
}

print(mean_rmse)
print(sprintf("Mean: %s", mean(mean_rmse)))
write(mean_rmse, 'error')
write(mean_rmse5, 'error5')
write(mean_rmse10, 'error10')
