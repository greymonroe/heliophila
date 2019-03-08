#ancestral state

require(phytools)
require(ape)
require(geiger)
require(extRemes)

tree <- read.nexus("../data/TrimmedMCCwoPolytomies.nex")
ladderize(tree) -> tree
rescale(tree, "depth", 1) -> tree

data <- read.csv("../data/species_means_table.csv", header=T)
data$A_P_num<-as.numeric(data$LH)-1

data.frame(data$A_P_num) -> df
rownames(df) <- gsub(" ","_", data$Species)
colnames(df) <- "Habit"
attach(df)

names(Habit) <- rownames(df)
name.check(tree, df)

mtreesER <- make.simmap(tree, Habit, model = "ER")
mtreesARD <- make.simmap(tree, Habit, model = "ARD")

lr.test(mtreesER$logL, mtreesARD$logL, df=100)

Maps <- make.simmap(tree, Habit, model = "ER", nsim=10000)
XX <- describe.simmap(Maps, plot = F)
pies<-cbind(data.frame(node=42+(1:nrow(XX$ace))), XX$ace)
write.csv(pies, "../data/ancestral_state.csv")

XX
         