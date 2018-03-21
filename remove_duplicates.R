rm(list=ls())
library("ape")

#define a function that take a tree and list with duplicate names 
#return a table with list of that me with loger branches.
Dup_name_branch <- function(tree, list){
  tip <- list
  index<- which(tree$tip.label %in% tip)
  Dup_List <- NULL  
  for(i in index){
    if(i %in% tree$edge[, 2]){
      tip <- paste0(tree$tip.label[i], spe="_", i)
      matrix <- cbind.data.frame(seq(1,length(tree$edge.length)),tree$edge,tree$edge.length)
      colnames(matrix) <- c("Index", "Pos_1", "Pos_2", "brlen")
      Loc <- matrix$Index[matrix$Pos_2==i]
      result <- cbind.data.frame(Tip=tip, brLegth=tree$edge.length[Loc])
    }
    Dup_List <-rbind.data.frame(Dup_List, result)
  }
    
  #rmtips <- Dup_List$Tip[Dup_List$brLegth==max(Dup_List$brLegth)]
 write.csv(Dup_List,"results/duplicate_names_with_branches2.csv")
 #return(rmtips)
#return(Dup_List)
}

CNtree <- ladderize(read.tree("data/Date_CN_ToL_genus_pruned.tre"))
du.names <- read.csv("data/360_duplicate.txt",stringsAsFactors = FALSE,  header=FALSE)
du.names <- du.names$V1
Dup_name_branch(CNtree, du.names)

data.file <- read.csv("results/duplicate_names_with_branches2.csv", header=TRUE)

#rm.list <- NULL
#keep.list <- NULL
for (name in du.names){
  Tt <- which.max(data.file$brLegth[grep(name, data.file$Tip)])
  #good.list <- as.character(data.file$Tip[grep(name, data.file$Tip)[Tt]])[[1]]
  bad.list <- as.character(data.file$Tip[grep(name, data.file$Tip)[-Tt]])
  write(bad.list, "results/duplcate_names_removed.txt", append=TRUE)
  #keep.list <- rbind(keep.list, good.list)
}

#write(rm.list, "results/duplcate_names_removed.txt")

#write(keep.list, "results/duplcate_names_keep.txt")

new.tips <- paste0(CNtree$tip.label, sep="_", seq(1,length(CNtree$tip.label)))

new.tree <- CNtree
new.tree$tip.label <- new.tips
write.tree(new.tree, "results/Date_CN_ToL_genus_pruned_seq.tre")

remve_tips <- read.csv("results/duplcate_names_removed.txt", header=FALSE)
remove_tips <- as.character(remve_tips$V1)

Tree_mono <- drop.tip(new.tree,remove_tips)
as.phylo(Tree_mono)

write.tree(Tree_mono, "results/Date_CN_ToL_genus_pruned_mono_seq.tre")
