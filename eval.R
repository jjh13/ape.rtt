library(ape)

#source('ape.patches.R')
source('rtt.R')
source('test.R')

trim <- function (x) gsub("^\\s+|\\s+$", "", x)
extract_dates <- function(x) as.numeric(gsub("(.+)_([0-9\\.]+)$", "\\2", x, perl=T))

n.simulated <- 1
mean_rmse <- rep(0, n.simulated)

# Clean up old stats files
unlink("__export_tmp.txt_phyml_tree.txt", force = T)
unlink("__export_tmp.txt_phyml_stats.txt", force = T)

for(i in 1:n.simulated) {
	simenv.dna <- read.FASTA(sprintf("simulated/HIV_%d_out_TRUE.fas", i))
	write.dna(simenv.dna , "__export_tmp.txt")
	phyml.env <- phymltest("__export_tmp.txt", execname = "~/Binaries/phyml")
	# , exclude=c("JC69+G", "JC69+I+G",
 #      "K80", "K80+I", "K80+G", "K80+I+G",
 #      "F81", "F81+I", "F81+G", "F81+I+G",
 #      "F84", "F84+I", "F84+G", "F84+I+G",
 #      "HKY85", "HKY85+I", "HKY85+G", "HKY85+I+G",
 #      "TN93", "TN93+I", "TN93+G", "TN93+I+G", "JC69", "JC69+I"
 #      , "GTR+I", "GTR+G", "GTR+I+G"))
	trees <- read.tree("__export_tmp.txt_phyml_tree.txt")
#	print trees
	tree <- trees[[28]]

	write.tree(tree, sprintf("simulated/tree/HIV_ml_%d_out.nwk", i))

	tip.dates <- trim(tree$tip.label)
	tip.dates <- extract_dates(tree$tip.label)

	mean_rmse[i] <- test.rtt(tree, tip.dates, remove = 20, random = T)

	unlink("__export_tmp.txt_phyml_tree.txt", force = T)
	unlink("__export_tmp.txt_phyml_stats.txt", force = T)
}

print(mean_rmse)
print(sprintf("Mean: %s", mean(mean_rmse)))
write(mean_rmse, 'error')
