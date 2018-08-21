#install.packages("devtools")
library("devtools")
install_github("josephwb/turboMEDUSA/MEDUSA", dependencies=TRUE)

library("MEDUSA")

richness <- read.csv("data/Genus.Richness.csv", header=TRUE)
tree <- read.tree("results/Date_CN_ToL_genus_pruned_mono.tre")
is.ultrametric(tree) #the phy tree need to be ultrametric

library("phytools")

tree.ult <- force.ultrametric(tree, method="extend")
tree.nnl <- force.ultrametric(tree, method="nnls")
is.ultrametric(tree.ult)
write.tree(tree.ult, "results/Dated_CN_ToL_genus_pruned_mono_ultra.tre")

par(mfrow=c(1,3))
plot.phylo(tree, show.tip.label = FALSE, main="Orignal")
plot.phylo(tree.ult, show.tip.label = FALSE, main="extend")
plot.phylo(tree.nnl, show.tip.label = FALSE, main="nnls")

#CNTree_medusa <- MEDUSA(phy=tree.ult, richness = richness)
CNTree_medusa <- readRDS("results/CNTree_MEDUSA.rds")
summ <- medusaSummary(CNTree_medusa)

save(tree.ult, richness, CNTree_medusa, summ, file="results/CNTree_MEDUSA.rds")


#this figure is still ugly
plotPrettyTree(summ, show.tip.label = FALSE, time=FALSE, node.labels = TRUE)
axisPhylo()

# and now we can summarize as though it were a distribution
results <- list(CNTree_medusa, CNTree_medusa)
summ.rates <- multiMedusaSummary(results, tree.ult)
write.csv(summ.rates$summary.tree$rates, "./results/CNTree_Medusa_rates.csv")
