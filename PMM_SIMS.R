library("MCMCglmm");library("brms");library("tidyverse");library("ape");
library("phytools");library("MASS");library("plyr");library("bindata");library('phangorn')


##------------------------------ SIMULATE TREES ---------------------------------##

n=50 # n taxa

trees <- list()
for (i in 1:2){
  
  t <- geiger::sim.bdtree(b=1, d=0, stop="taxa", n=2, extinct=FALSE) # backbone
  a <- geiger::sim.bdtree(b=1, d=0, stop="taxa", n=n, extinct=FALSE) # clade A
  b <- geiger::sim.bdtree(b=1, d=0, stop="taxa", n=n, extinct=FALSE) # clade B
  
  b$tip.label <- paste0("s",51:100)
  
  # re-scale both clades to height 1
  a$edge.length<-a$edge.length/max(nodeHeights(a)[,2])*0.90
  b$edge.length<-b$edge.length/max(nodeHeights(b)[,2])*0.90
  
  # duplicate backbone and re-scale to short and long
  t.1 <- t;t.1$edge.length<-t.1$edge.length/max(nodeHeights(t.1)[,2])*0.1
  t.2 <- t;t.2$edge.length<-t.2$edge.length/max(nodeHeights(t.2)[,2])*1.1
  
  # bind clades onto backbones
  t1 <- bind.tree(t.1, a, where = 1, position = 0)
  t1 <- bind.tree(t1, b, where = 1, position = 0)
  t2 <- bind.tree(t.2, a, where = 1, position = 0)
  t2 <- bind.tree(t2, b, where = 1, position = 0)
  
  trees[[i]] <- list(t1,t2)
}


##------------------------------ SIMULATE DATA ----------------------------------##

## define parameters

# constants
beta = c(0,0) # intercepts
m = 2 # m traits

# variables
vars <- data.frame(mod.evo = c("BM1","BM2"),
                   b11 = c(0.75,0.75),
                   b22 = c(0.75,0.75),
                   b12_rho = c(0,0.5),
                   c11 = c(0.25,0.25),
                   c22 = c(0.25,0.25),
                   c12_rho = c(0,0.5))


## create list to store model fits
# list nesting = [[tree.rep]][[tree.type]][[model.evo]] - HOW TO NAME LIST ELEMENTS PROPERLY?
fits <- list(rep(NA,length(trees)))
for (i in 1:length(trees)){fits[[i]] <- list(NA,NA)}
for (i in 1:length(trees)){fits[[i]][[1]] <- list(NA,NA);fits[[i]][[2]] <- list(NA,NA)}

## LOOP OVER TREE REPS, TREE TYPES and MODEL OF EVOLUTION
for (i in 1:length(trees)){
  for (j in 1:length(trees[[1]])){
    for (k in 1:length(vars$mod.evo)){
      
      tree <- trees[[i]][[j]]
      
      # Create VCV matrix from tree, scale to correlation matrix for BRMS
      A.mat <- ape::vcv.phylo(tree, corr = T) # scale with corr = T
      
      sig.B <- c(b11 = vars[vars$mod.evo==mod.evo[k],]$b11, b22 = vars[vars$mod.evo==mod.evo[k],]$b22)
      Bcor <- matrix(c(c(1,vars[vars$mod.evo==mod.evo[k],]$b12_rho),c(vars[vars$mod.evo==mod.evo[k],]$b12_rho,1)),m,m, byrow = T) 
      B <- matrix(kronecker(sig.B, sig.B),m,m)*Bcor # VCV (point-wise product). 
      
      sig.C <- c(c11 = vars[vars$mod.evo==mod.evo[k],]$c11, c22 = vars[vars$mod.evo==mod.evo[k],]$c22)
      Ccor <- matrix(c(c(1,vars[vars$mod.evo==mod.evo[k],]$c12_rho),c(vars[vars$mod.evo==mod.evo[k],]$c12_rho,1)),m,m, byrow = T) 
      C <- matrix(kronecker(sig.C, sig.C),m,m)*Ccor 
      
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
      species <- trees[[1]][[1]]$tip.label
      clade <- c(rep("A",n/2),rep("B",n/2))
      d <- data.frame(species,clade,y1,y2)
      d$animal <- d$species # "animal" reserved term in MCMCglmm 
      d$obs <- 1:nrow(d)
      
      fit <- brm(
        mvbind(y1, y2) ~ (1|p|gr(animal, cov = A)),
        data = d,
        data2 = list(A = A.mat),
        family = gaussian(),
        cores = 4,
        chains = 4, iter = 1000, thin = 1
      )
      
      fits[[i]][[j]][[k]] <- fit
      
    }
  }
}


# saveRDS(fits, "fits.rds")
fits <- readRDS("fits.rds")


##------------------------------ RESULTS ----------------------------------##

## LEGEND
# list nesting = [[tree.rep]][[tree.type]][[model.evo]]
# tree.rep = 1:ntrees
# tree.type = c(early split, balanced)
# model.evo = c(BM1, BM2) (PRICE to be added)

# generating parameters for reference
vars
## compare fits to BM1 and BM2 from the same tree
summary(fits[[1]][[1]][[1]])
summary(fits[[1]][[1]][[2]])
## compare fits to BM1 from balanced and early split trees
summary(fits[[1]][[1]][[1]]);plot(trees[[1]][[1]]) # balanced
summary(fits[[1]][[2]][[1]]);plot(trees[[1]][[2]]) # early split
## compare fits to BM2 from balanced and early split trees
summary(fits[[1]][[1]][[2]])
summary(fits[[1]][[2]][[2]])


## TO DO / RESOLVE ##

# 1. WHY ARE MEANS FOR CLADES A AND B NOT MORE DIFFERENT IN OUR SIM DATA?

# means of clade A and clade B
mean(d[d$clade=="A",]$y1);mean(d[d$clade=="B",]$y1)
# sig difference between means of clade A and clade B?
t.test(d[d$clade=="A",]$y1,d[d$clade=="B",]$y1)
# compare to difference when BM simulated directly on the trees - sig2 != b_11?
x<-fastBM(trees[[1]][[1]],sig2=1)
mean(x[1:50])-mean(x[51:100])
y<-fastBM(trees[[1]][[2]],sig2=1)
mean(y[1:50])-mean(y[51:100])



##--------------------------- NOTES --------------------------------##

# final tree heights and plot first pair of trees in list
nodeheight(trees[[1]][[1]],1);plot(trees[[1]][[1]])
nodeheight(trees[[1]][[2]],1);plot(trees[[1]][[2]])

# # match heights of clades by adding difference to shorter of the pair
# x <- c(nodeheight(t1,1)-nodeheight(t1,n+1),nodeheight(t1,n+1)-nodeheight(t1,1))
# if (which(x<=0)==1){t1$edge.length[1] <- t1$edge.length[1]+x[which(x>=0)]} else {t1$edge.length[100] <- t1$edge.length[100]+x[which(x>=0)]}

