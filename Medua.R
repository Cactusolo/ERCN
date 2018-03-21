#install.packages("devtools")
library("devtools")
install_github("josephwb/turboMEDUSA/MEDUSA", dependencies=TRUE)

library("MEDUSA")

richness <- read.csv("data/Genus.Richness.csv", header=TRUE)
tree <- read.tree("results/Date_CN_ToL_genus_pruned_mono.tre")
is.ultrametric(tree)

library("phytools")

tree.ult <- force.ultrametric(tree, method="extend")
tree.nnl <- force.ultrametric(tree, method="nnls")
is.ultrametric(tree.ult)
write.tree(tree.ult, "results/Dated_CN_ToL_genus_pruned_mono_ultra.tre")
par(mfrow=c(1,3))
plot.phylo(tree, show.tip.label = FALSE, main="Orignal")
plot.phylo(tree.ult, show.tip.label = FALSE, main="extend")
plot.phylo(tree.nnl, show.tip.label = FALSE, main="nnls")

rosids_medusa <- MEDUSA(phy=tree.ult, richness = richness)
rosids_medusa2 <- MEDUSA(phy=tree, richness = richness,criterion="aic",stop="threshold")

summ <- medusaSummary(rosids_medusa)
summ2 <- medusaSummary(rosids_medusa2)

saveRDS(rosids_medusa,"results/MEDUSA.rds")
saveRDS(sum,"results/CNTree_MEDUSA_sum.rds")
