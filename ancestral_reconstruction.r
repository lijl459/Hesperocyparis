library(ape)
library(phytools)
setwd("F:/Desktop/shift")
getwd()

###########
# Environmental PCA

trait = read.csv("AI.csv", header =T, sep = ",")
head(trait)
trait.vector = trait$AI
names(trait.vector) <- trait$species

tree = read.tree("17sp_time_tree.tre")
tree <- ladderize(tree)

figure = contMap(tree, trait.vector, plot = FALSE)
plot(figure, type = "phylogram",legend=0.7*max(nodeHeights(tree)),fsize=c(0.7,0.9))

# To show tip labels 
plot(figure, type = "phylogram",legend=0.7*max(nodeHeights(tree)),fsize=c(0.7,0.9))
# To show nice small node labels
nodelabels(text=tree$edge,node=tree$edge, cex = 0.05, frame = "none")

#############
# Find major climatic shifts by percentile
# tree$edge has columns: parent, child

anc_rec = fastAnc(tree, trait.vector)
# Don't use ape -- node numbering doesn't match this

edges = data.frame(tree$edge)
	
diffs <- data.frame(NA,NA)
names(diffs) <- c("node", "value")

# Excludes terminals because these return NAs which we redefine as zero jump
i = 1
# doesn't work due to node order
for (i in 1:length(edges[,1])) {
	#print(edges[i,1])
	#print(edges[i,2])
	#print(anc_rec[toString(edges[i,2])])
	#print(anc_rec[toString(edges[i,1])])
	shifts = c(toString(edges[i,2]), abs(as.numeric(anc_rec[toString(edges[i,2])]) - as.numeric(anc_rec[toString(edges[i,1])])))
	print(shifts)
	diffs <- rbind(diffs, shifts)
	i = i + 1
	}
	
#diffs <- na.omit(diffs)
diffs[is.na(diffs)] <- 0
diffs <- as.data.frame(matrix(unlist(diffs), nrow = length(unlist(diffs[1]))), stringsAsFactors = FALSE)
names(diffs) <- c("node", "value")
row.names(diffs) = diffs$node
library(gtools)
attach(diffs)
diffs <- diffs[mixedsort(node), ] # Sort strings numerically using gtools mixedsort
detach(diffs)
diffs$node <- NULL
# E.g., slice by row name: diffs["2909",]


diffs1 = diffs[,1]
diffs1 = as.numeric(diffs1)
names(diffs1) = row.names(diffs)

# To prove this difference is correct, these should produce the same result:
# FOR CLIMATE MATCHED TREE ONLY
anc_rec["2673"] - anc_rec["2670"] # Node 2673 minus its parent
diffs1["2673"] # Absolute value of difference of node 2673 from its parent
# FOR CLIMATE AND PHENOTYPE MATCHED TREE ONLY
anc_rec["2702"] - anc_rec["2701"] # Node 2673 minus its parent
diffs1["2702"] # Absolute value of difference of node 2673 from its parent
anc_rec["2679"] - anc_rec["2678"] # Node 2673 minus its parent
diffs1["2679"] # Absolute value of difference of node 2673 from its parent



# Plot top 95% of shifts
diffs_top95 <- which(diffs1 > quantile(diffs1, .95))
figure = contMap(tree, trait.vector, plot = FALSE)

col <- colorRampPalette(c("#FF0000","#EE166E","#DA229D","#8D38EB","#0042FF"))(5)

obj2 <- setMap(figure, col);

plot(obj2, type = "phylogram",legend=0.7*max(nodeHeights(tree)),fsize=c(0.7,0.9))
nodelabels(node = diffs_top95, pch = 21, cex=0.2)


plot(obj2, type = "phylogram",legend=0.7*max(nodeHeights(tree)),fsize=c(0.7,0.9))
nodelabels(node = diffs_top95, pch = 19, cex=1)


# To check that this is behaving correctly:
# print quantile
quantile(diffs1, .95)
# The values returned below should be above the quantile
for (x in diffs_top95) {
print(diffs1[x])
}

#######################################
# Get timing of top 95% shifts

times <- data.frame(NA,NA)
names(times) <- c("node", "nodeheight")
for (i in diffs_top95) {
	print(i)
	item = c(i, nodeheight(tree, i))
	print(item)
	times <- rbind(times, item)
	}
times <- na.omit(times)
times$nodeheight = 47 - times$nodeheight

write.table(times, "AI_niche_shift_timings.txt")
# Adjust format manually

phylosig(tree, AI, test=TRUE, nsim=10000)
