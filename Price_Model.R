# implements the idea in Price for simulating niche conservatism

library(ape)
library(phangorn)
library(MASS)
library(apTreeshape)
library(tidyverse)
library(mvtnorm)

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

SPECIES=2
niche.space <- mvrnorm(n = SPECIES, rep(0, 2), Sigma)

plot(atree)
t <- data.frame(x=niche.space[,1],y=niche.space[,2])
mut <- mvrnorm(n = n, rep(0,n), 1)





atree <- rtree(2);atree$edge.length <- c(0,0) # backbone tree
num.tips = 2 # start with 2 species
m <- 2 # number of traits
N <- 1000000 # total final species count
t <- matrix(NA,N,2) # species trait values
t[1,] <- mvrnorm(1, rep(0, m), Sigma)
t[2,] <- mvrnorm(1, rep(0, m), Sigma)

# grow tree to size N by adding one tip per time step. mutate tip states at each time step and weight candidate
# states (species) by
for (i in 2:N) {

e = rnorm(2*num.tips,0,1) %>% matrix(num.tips) # mutation (independent additive genetic variances)
t_star = t[1:num.tips,] + e # create candidate species by mutating current species
w <- dmvnorm(t_star, rep(0,m), Sigma, log=FALSE) # calculate weights for candidates from density of MVN(0,Sigma)
new.sp <- sample(1:nrow(t_star), 1, replace = T, prob=w) # take weighted sample from candidates
t[num.tips + 1,] <- t_star[new.sp,] # add new species trait values to t

  # wait some random amount of time
  # i.e. grow all tips
  growth <- rep(0, nrow(atree$edge))
  growth[atree$edge[,2]<=num.tips] <-rep(rexp(1), num.tips)
  atree$edge.length <- atree$edge.length + growth
  
  # add new species as a bifurcation
  newlabel <- paste("t",as.character(num.tips + 1), sep="")
  atree <- add.tips(atree, newlabel, new.sp)
  num.tips = num.tips + 1

}

par(mfrow = c(1,2))
plot(atree);plot(t, type="n"); text(t, atree$tip.label, cex=0.8)

cor(t)

# calculate bivariate normal density (Sigma) for t_star to give weights
# take a weighted random choice of the Nk candidates
# iterate


