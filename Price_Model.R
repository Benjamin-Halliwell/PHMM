
## PRICE MODEL
# implements the 'adaptive radiation' model in Price 1997 for simulating the evolution of
# bivariate trait data with niche conservatism. covariance Sigma defines the ellipse in 
# trait space representing competent combinations of trait1 and trait2

#----------------------------------------------------------------------

library(ape)
library(phangorn)
library(MASS)
library(apTreeshape)
library(tidyverse)
library(mvtnorm)
library(phytools)
library(mixtools)


force.ultrametric<-function(tree,method=c("nnls","extend")){
  method<-method[1]
  if(method=="nnls") tree<-nnls.tree(cophenetic(tree),tree,
                                     rooted=TRUE,trace=0)
  else if(method=="extend"){
    h<-diag(vcv(tree))
    d<-max(h)-h
    ii<-sapply(1:Ntip(tree),function(x,y) which(y==x),
               y=tree$edge[,2])
    tree$edge.length[ii]<-tree$edge.length[ii]+d
  } else 
    cat("method not recognized: returning input tree\n\n")
  tree
}

#----------------------------------------------------------------------

# TWO DIFFERENT METHODS:

## 1. NICHE OPENING METHOD ##

# grow tree to size N by adding one tip per time step. Draw state of new tip (open new niche) from mvnorm with 
# covariance defined by Sigma. calculate euclidean distance between tip states. attach new tip to existing tip 
# with closest state (simulates niche conservatism). N.B. all tips states drawn and euc dists calculated
# up front, but same result as drawing from mvnorm and calculating min dist at each time step.

N = 15 # number of species
Sigma <- matrix(c(1,0.8,0.8,1),2,2) # parameters for mvnorm
reps = 2 # generate 2 clades to stitch together
trees <- list() # list of trees generated
t.list <- list() # list of species trait data (niches)

for (j in 1:reps){

# randomly generate 2D points from a multivariate normal using Sigma as VCV
# matrix and set mean to 0
niche.space <- mvrnorm(n = N, rep(0, 2), Sigma) 
niche.dist <- as.matrix(dist(niche.space)) # euclidean distance between N in niche space

# for each entry find the tip with a lower index that is closest. 
# tips added sequentially, start at 3 as initial tree must be a cherry
parent <- rep(0,N)
for (i in 3:N){
  temp <- niche.dist[i,1:(i-1)]
  parent[i] <- which.min(temp)
  }

# init a cherry tree
atree <- rtree(2)
atree$edge.length <- c(0,0)
num.tips = 2

while(num.tips < N){
  
  # wait some random amount of time, i.e. grow all tips
  growth <- rep(0, nrow(atree$edge))
  growth[atree$edge[,2]<=num.tips] <-rep(rexp(1), num.tips)
  atree$edge.length <- atree$edge.length + growth

  # add a tip (niche filled by speciation from existing tip with most similar niche)
  node <- parent[num.tips + 1] # parent node defined by euc dists generated above
  newlabel <- paste("t",as.character(num.tips + 1), sep="")
  atree <- add.tips(atree, newlabel, node)
  num.tips = num.tips + 1

  }

t.list[[j]] <- niche.space
trees[[j]] <- atree

}

# create backbone
phy <- rtree(2); phy$edge.length <- c(1,1)
a <- trees[[1]] # clade A
b <- trees[[2]] # clade B
b$tip.label <- paste0("t", (N+1):(2*N))
# merge trait data
t <- rbind(t.list[[1]],t.list[[2]])

# re-scale both clades to height 1
a$edge.length<-a$edge.length/max(nodeHeights(a)[,2]);nodeHeights(a)
b$edge.length<-b$edge.length/max(nodeHeights(b)[,2]);nodeHeights(b)

# duplicate backbone and re-scale to short and long
phy.1 <- phy;phy.1$edge.length<-phy.1$edge.length/max(nodeHeights(phy.1)[,2])*0.1
phy.2 <- phy;phy.2$edge.length<-phy.2$edge.length/max(nodeHeights(phy.2)[,2])

# bind clades onto backbones
phy1 <- bind.tree(phy.1, a, where = 1, position = 0);phy1 <- bind.tree(phy1, b, where = 1, position = 0)
phy2 <- bind.tree(phy.2, a, where = 1, position = 0);phy2 <- bind.tree(phy2, b, where = 1, position = 0)

# plot trees and traits
par(mfrow = c(1,2))
plot(phy1);plot(t, type="n");text(t, phy1$tip.label, cex=0.8) # short basal branches
# plot(phy2);plot(t, type="n");text(t, phy1$tip.label, cex=0.8) # long
cor(t)



## 2. MUTATION METHOD ##

# grow tree to size N by adding one tip per time step. mutate tip states at each time step and 
# weight candidate states (species) by joint probability density. sig.scale scales the magnitude 
# of mutational genetic (co)variance, determines the speed and pattern by which trait space is 
# explored from starting species values

# Two ways of incorporating extinction into either method:
# A = with probability p at each time step, random tip pruned (goes extinct)
# B = tree pruned to size N*(1-prop) at the end
p <- 0.05 # prob of extinction event each time step
prop <- 0.2 # prop of tips randomly pruned

N <- 25 # final species count
m <- 2 # number of traits
reps <- 2 # 2 clades
Sigma <- matrix(c(1,0,0,1),2,2) # parameters for mvnorm
sig <- matrix(c(1,-0.75,-0.75,1),2,2) # non-zero off-diagonals to simulate genetic covariance between t1 and t2
sig.scale <- 0.75 # scale genetic (co)variance term
trees <- list()
t.list <- list()

for (j in 1:reps){

  t <- matrix(NA,N,2) # mat to store species trait values
  t[1,] <- mvrnorm(1, rep(0, m), Sigma) # draw trait values for first two species
  t[2,] <- mvrnorm(1, rep(0, m), Sigma)
  atree <- rtree(2);atree$edge.length <- c(0,0) # backbone tree
  num.tips = 2 # start with 2 species
  
while (num.tips < N) {

  e = mvrnorm(num.tips, rep(0, m), sig) # mutation (additive (co)variance)
  e <- e * sig.scale
  t_star = t[1:num.tips,] + e # create candidate species by mutating current species
  w <- dmvnorm(t_star, rep(0,m), Sigma, log=FALSE) # calculate weights for candidates from density of MVN(0,Sigma)
  new.sp <- sample(1:nrow(t_star), 1, replace = T, prob=w) # take weighted sample from candidates (i.e. penalise mutant phenotypes that fall too far out of the ellipse)
  t[num.tips + 1,] <- t_star[new.sp,] # add new species trait values to t

  # wait some random amount of time, i.e. grow all tips
  growth <- rep(0, nrow(atree$edge))
  growth[atree$edge[,2]<=num.tips] <-rep(rexp(1), num.tips)
  atree$edge.length <- atree$edge.length + growth

  # add new species as a bifurcation to parent node
  newlabel <- paste("t",as.character(num.tips + 1), sep="")
  atree <- add.tips(atree, newlabel, where = new.sp)
  num.tips = num.tips + 1

  ## Extinction A
  # if(rbinom(1,1,prob=p)==1){
  #   drop <- sample(atree$tip.label, 1, prob=rep(0.01, length(atree$tip.label)))
  #   atree <- drop.tip(atree, drop, trim.internal = T)
  # } else NULL
  
  atree <- force.ultrametric(atree)
  num.tips = length(atree$tip.label)
}

  t.list[[j]] <- t
  trees[[j]] <- atree
  
}

## Extinction B
# drop <- sample(atree$tip.label, N*prop, prob=rep(0.01, length(atree$tip.label)))
# drop2 <- as.numeric(gsub("t", "", drop))
# rtree <- drop.tip(atree, drop, trim.internal = T)
# par(mfrow = c(1,2))
# plot(rtree, cex=0.5);plot(t, type="n"); text(t[-c(drop2),], rtree$tip.label, cex=0.5)
# cor(t[-c(drop2),])

# stitch trees together to create long/short basal branches
phy <- rtree(2); phy$edge.length <- c(1,1)
a <- trees[[1]] # clade A
b <- trees[[2]]# clade B
b$tip.label <- paste0("t", (N+1):(2*N))
t <- rbind(t.list[[1]],t.list[[2]])

# re-scale both clades to height 1
a$edge.length<-a$edge.length/max(nodeHeights(a)[,2])
b$edge.length<-b$edge.length/max(nodeHeights(b)[,2])

# duplicate backbone and re-scale to short and long
phy.1 <- phy;phy.1$edge.length<-phy.1$edge.length/max(nodeHeights(phy.1)[,2])*0.1
phy.2 <- phy;phy.2$edge.length<-phy.2$edge.length/max(nodeHeights(phy.2)[,2])

# bind clades onto backbones
phy1 <- bind.tree(phy.1, a, where = 1, position = 0);phy1 <- bind.tree(phy1, b, where = 1, position = 0)
phy2 <- bind.tree(phy.2, a, where = 1, position = 0);phy2 <- bind.tree(phy2, b, where = 1, position = 0)

# plot trees and traits
par(mfrow = c(1,2))
plot(phy1, cex=0.5, label.offset = 0.01); plot(matrix(c(3,-3,3,-3),2,2), type="n");text(t, phy1$tip.label, cex=0.8) # short basal branches
# plot(phy2);plot(t, type="n");text(t, phy1$tip.label, cex=0.8) # long

# plot ellipses representing constraints from ecological (black) and genetic (red) covariances
ellipse(0, Sigma, alpha = 0.05, npoints = 250, col = "black")
ellipse(0, sig.scale*sig, alpha = 0.05, npoints = 250, col = "red")


cor(t)




