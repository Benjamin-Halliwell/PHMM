# implements the idea in Price for simulating niche conservatism

library(ape)
library(phangorn)
library(MASS)
library(apTreeshape)

SPECIES = 15 #number of species
Sigma <- matrix(c(1,0.8,0.8,1),2,2) #parameters for mvnorm
reps = 2
p.trees <- list()
niches <- list()

for (j in 1:reps){

# randomly generate 2D points from a multivariate normal using Sigma as VCV
# matrix and set mean to 0
niche.space <- mvrnorm(n = SPECIES, rep(0, 2), Sigma) 
niche.dist <- as.matrix(dist(niche.space)) # euclidean distance between species in niche space

# for each entry find the tip with a lower index that is closest
# start at 3 as initial tree must be a cherry
parent <- rep(0,SPECIES)
for (i in 3:SPECIES){
  temp <- niche.dist[i,1:(i-1)]
  parent[i] <- which.min(temp)
  }

#init a cherry tree
atree <- rtree(2)
atree$edge.length <- c(0,0)
num.tips = 2

while(num.tips < SPECIES){
  # wait some random amount of time
  # i.e. grow all tips
  growth <- rep(0, nrow(atree$edge))
  growth[atree$edge[,2]<=num.tips] <-rep(rexp(1), num.tips)
  atree$edge.length <- atree$edge.length + growth
  #plot(atree)
  
  #pick a random node
  #rnode <- sample(1:(atree$Nnode+1),1)
  rnode <- parent[num.tips + 1] # euc dists and closest species are generated above
  
  #add a tip
  newlabel <- paste("t",as.character(num.tips + 1), sep="")
  atree <- add.tips(atree, newlabel, rnode)
  num.tips = num.tips + 1

  }

niches[[j]] <- niche.space
p.trees[[j]] <- atree

}

t <- rtree(2); t$edge.length <- c(1,1)
a <- p.trees[[1]] # clade A
b <- p.trees[[2]] # clade B
b$tip.label <- paste0("t", 16:30)

plot(t)

# re-scale both clades to height 1
a$edge.length<-a$edge.length/max(nodeHeights(a)[,2]);nodeHeights(a)
b$edge.length<-b$edge.length/max(nodeHeights(b)[,2]);nodeHeights(b)

# duplicate backbone and re-scale to short and long
t.1 <- t;t.1$edge.length<-t.1$edge.length/max(nodeHeights(t.1)[,2])*0.1
t.2 <- t;t.2$edge.length<-t.2$edge.length/max(nodeHeights(t.2)[,2])

# bind clades onto backbones
t1 <- bind.tree(t.1, a, where = 1, position = 0)
t1 <- bind.tree(t1, b, where = 1, position = 0)
t2 <- bind.tree(t.2, a, where = 1, position = 0)
t2 <- bind.tree(t2, b, where = 1, position = 0)

plot(t1);plot(data.frame(min(min(niches[[1]][,1]), min(niches[[2]][,1])),), type="n"); text(rbind(niches[[1]],niches[[2]]), t1$tip.label, cex=0.8)
plot(t2);plot(niches[[1]], type="n"); text(rbind(niches[[1]],niches[[2]]), t1$tip.label, cex=0.8)

par(mfrow = c(1,2))
plot(atree)
plot (niche.space, type="n")
text(niche.space, atree$tip.label, cex=0.8)
atree$tip.label


