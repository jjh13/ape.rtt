library(ape)
library(phyloch)
library(parallel)

#source('ape.patches.R')
source('rtt.R')
source('test.R')

trim <- function (x) gsub("^\\s+|\\s+$", "", x)
extract_dates <- function(x) as.numeric(gsub("(.+)_([0-9\\.]+)$", "\\2", x, perl=T))


n.simulated <- 2
nfp <- c(1, 2, 2, 3, 2, 3, 3, 4, 4, 5, 5, 6, 5, 6, 6, 7, 5, 6, 6, 7, 6, 7, 7, 8, 9, 10, 10, 11)
exclude <- c(
	"JC69", 
	"JC69+I", "JC69+G", "JC69+I+G",
      "K80", "K80+I", "K80+G", "K80+I+G",
      "F81", "F81+I", "F81+G", "F81+I+G",
      "F84", "F84+I", "F84+G", "F84+I+G",
      #"HKY85", 
      #"HKY85+I", 
      "HKY85+G", "HKY85+I+G",
      "TN93", "TN93+I", "TN93+G", "TN93+I+G",
      "GTR", "GTR+I", "GTR+G", "GTR+I+G")

ml.tree <- function(i){
	dna.file <- sprintf("__tmp%d.dna",i)
	unlink(sprintf("__tmp%d.dna_phyml_tree.txt", i), force = T)
	unlink(sprintf("__tmp%d.dna_phyml_stats.txt", i), force = T)
	unlink(dna.file)

#	simenv.dna <- read.FASTA(sprintf("simulated/HIV_%d_out_TRUE.fas", i))
	simenv.dna <- read.dna(sprintf("simulated/HIV_SIM_%d.phy", i))
#	uid <- sample(1:length(names(simenv.dna)), length(names(simenv.dna)), replace=F)
	# read.FASTA is buggy, and creates the DNABin structure incorrectly
	# BUT, it's so broken that we can mess around with the names, then...
#	for(j in 1:length(names(simenv.dna))) { 
		# give each bin a unique name 
#		names(simenv.dna)[j] <- sprintf("u%s_%s", uid[j], names(simenv.dna)[j])
#	}

	# ... output a file (which works for some reason)
	write.dna(simenv.dna, dna.file)
	# and re-read it as a properly formatted DNA Bin
	simenv.dna <- read.dna(dna.file)

	# only to never re-use that dna bin ever again...
	phyml.env <- phymltest(dna.file, execname = "~/Binaries/phyml")#, exclude=exclude)

	aic <- 2 * (nfp - phyml.env)
	trees <- read.tree(sprintf("__tmp%d.dna_phyml_tree.txt", i))

	# Choose the tree with the min aic
	tree <- trees[[which(aic==min(aic))]]

	write.tree(tree, sprintf("simulated/tree/HIV_ml_%d_out.nwk", i))
	unlink(sprintf("__tmp%d.dna_phyml_tree.txt", i), force = T)
	#unlink(sprintf("__tmp%d.dna_phyml_stats.txt", i), force = T)
	unlink(dna.file)
	tree
}

trees.read <- function(base.path) {
	ml.tree.read <- function(i, path){
		read.tree(sprintf("%s/HIV_ml_%d_out.nwk", path, i))
	}
	mclapply(1:n.simulated, ml.tree.read, base.path)
}

trees <- mclapply(1:n.simulated, ml.tree, mc.cores=12)
#trees <- trees.read("simulated/tree/")
run_name <- "coal"
species <- "Simulated"
n.runs <- 1

for(remove in c(1,5,10)){
	pdf(sprintf('%s_%d.pdf', run_name, remove), width=11.5, height=8.5)

	mse <- rep(0, n.simulated*n.runs)
	for(i in 1:n.simulated){
		tip.dates <- extract_dates(trees[[i]]$tip.label)
		test.n.plot.rtt(trees[[i]], tip.dates, remove, name=sprintf("%s: %s-no.%d (Remove %d)", run_name, species, i, remove))
		for(j in 1:n.runs)
		 	mse[(j-1)*n.simulated + i] <- test.rtt(trees[[i]], tip.dates, remove)
	}

	par(mfrow = c(1, 1))
	mean_mse <- mean(mse)
	hist(mse, breaks=100,  main="Histogram of MSE", sub=sprintf("(Average MSE [in red]: %s)", signif(mean_mse, digits=6)))
	abline(v = mean_mse, col = "red", lwd = 2)
	dev.off()
}


