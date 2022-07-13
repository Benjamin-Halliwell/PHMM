library(geiger)
library(MASS)
library(phytools)
library(phangorn)
library(phylolm)
library(brms)

# return tree and traits as list
sim_Price <- function (N, Sigma) {
  
  # randomly generate 2D points from a multivariate normal using Sigma as VCV
  # matrix and set mean to 0
  niche.space <- mvrnorm(n = N, rep(0, 2), Sigma) 
  niche.dist <- as.matrix(dist(niche.space)) # euclidean distance between N in niche space
  
  # for each entry find the tip with a lower index that is closest. 
  # tips added sequentially, start at 3 as initial tree must be a cherry
  parent <- rep(0,N)
  for (k in 3:N){
    temp <- niche.dist[k,1:(k-1)]
    parent[k] <- which.min(temp)
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
  
  # wait some random amount of time, i.e. grow all tips
  growth <- rep(0, nrow(atree$edge))
  growth[atree$edge[,2]<=num.tips] <-rep(rexp(1), num.tips)
  atree$edge.length <- atree$edge.length + growth
  
  list(niche.space, atree)

  
}

N <-  10
Sigma <- matrix(c(1,0.5,0.5,1),2,2)

tt <- sim_Price(N, Sigma)

tree <- tt[[2]]
trait <- as.data.frame(tt[[1]]);names(trait) <- c("x","y");rownames(trait) <- tree$tip.label


#### FITS ####

## PGLS
fit.pgls <- phylolm(y ~ x, phy=tree, model = "lambda", data=trait)


## BRMS
# global pars
m=2
beta = c(0,0)
k <- 1

vars <- data.frame(mod.evo = c("BM1","BM2","BM3","BM4","Price_0.75"),
                   b11 = c(0.75,0.75,0.75,0.75, NA),
                   b22 = c(0,0.75,0.75,0.75, NA),
                   b12_rho = c(0,0,0.5,0.5, NA),
                   c11 = c(0.25,0.25,0.25,0.25, NA),
                   c22 = c(0.25,0.25,0.25,0.25, NA),
                   c12_rho = c(0,0,0,0.5,NA))


sig.B <- c(b11 = vars[vars$mod.evo==vars$mod.evo[k],]$b11, b22 = vars[vars$mod.evo==vars$mod.evo[k],]$b22)
Bcor <- matrix(c(c(1,vars[vars$mod.evo==vars$mod.evo[k],]$b12_rho),c(vars[vars$mod.evo==vars$mod.evo[k],]$b12_rho,1)),m,m, byrow = T) 
B <- matrix(kronecker(sig.B, sig.B),m,m)*Bcor # VCV (point-wise product). 

# residual covariance matrix
sig.C <- c(c11 = vars[vars$mod.evo==vars$mod.evo[k],]$c11, c22 = vars[vars$mod.evo==vars$mod.evo[k],]$c22)
Ccor <- matrix(c(c(1,vars[vars$mod.evo==vars$mod.evo[k],]$c12_rho),c(vars[vars$mod.evo==vars$mod.evo[k],]$c12_rho,1)),m,m, byrow = T) 
C <- matrix(kronecker(sig.C, sig.C),m,m)*Ccor 

# sim dat
A.mat <- vcv.phylo(tree, corr = T) 

# # ensure A.mat is positive definite (need to extend to ensure matrix is non-singular?)
# while (dim(table(eigen(A.mat)$values > 0))==2) { # should have all eigenvalues > 0
#   A.mat <-  A.mat + diag(1e-6,nrow(A.mat)) } # if there are zero values, adjust by adding a small constant to the diagonal

# A = phylogenetic taxon-level VCV matrix
A = A.mat
# number of species
n = nrow(A)
# identity matrix
I = diag(n)

# simulate phylogenetic random effects for all traits as one draw from a MVN
a <- mvrnorm(1,rep(0,n*m),kronecker(B,A))
a1 <- a[1:n]
a2 <- a[1:n + n]

# simulate residuals
e <- mvrnorm(1,rep(0,n*m),kronecker(C,I))
e1 <- e[1:n]
e2 <- e[1:n + n]

# construct response traits from each linear predictor
y1 <- beta[1] + a1 + e1
y2 <- beta[2] + a2 + e2

# generate df
d <- data.frame(animal=tree$tip.label,y1,y2)

fit.brms <-  brm(
                bf(mvbind(y1, y2) ~ (1|p|gr(animal, cov = A))) + set_rescor(TRUE),
                data = d,
                data2 = list(A = A),
                family = gaussian(),
                fit = fit.brms,
                chains = 1)

### END ####
               
fit.pgls
fit.brms
