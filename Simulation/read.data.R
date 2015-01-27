dat <- read.phyDat("outputname_TRUE.fas", format="fasta", type="AA")
dists <- dist.ml(dat, model="HIVw")


