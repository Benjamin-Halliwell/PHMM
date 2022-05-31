library("MCMCglmm");library("brms");library("tidyverse");library("ape");
library("phytools");library("MASS");library("plyr");library("bindata");library('phangorn')


##---------------------------- ~ --------------------------------##

#### DEFINE SIMULATION PARAMETERS ####


## SIMULATE TREE ##

# simulate tree
seed <- 58198 # sample(1e6,1); seed
# seed 58198
t1 <- geiger::sim.bdtree(b=1, d=0, stop="taxa", n=50, seed = seed, extinct=FALSE) # sim until we get one with two major clades
plot(t1)
#t1 <- multi2di(t1)

# create duplicate with long basal edges
t2 <- t1
adj <- t2$edge.length[1]-t2$edge.length[50]
t2$edge.length[50] <- 5-adj
t2$edge.length[1] <- 5

# compare trees
plot(t1)
plot(t2)


## SIMULATE DATA ##

# Create VCV matrix from tree, scale to correlation matrix for BRMS
A.mat <- ape::vcv.phylo(t2, corr = T) # scale with corr = T
A.mat

# ## alternative method 
# inv.phylo <- MCMCglmm::inverseA(t.toy, nodes = "TIPS", scale = TRUE)
# A <- solve(inv.phylo$Ainv)
# rownames(A) <- rownames(inv.phylo$Ainv)
# A.mat <- A

# number of traits
k = 2

# intercepts for each trait - plogis(beta[1]) to get on probability scale
beta = c(1,2)

## B is the phylogenetic trait-level VCV matrix (Sigma_{phy} in the HTML). B specifies the phylogenetic
## variance in each trait (diagonals) as well as the phylogenetic covariance between traits (off-diagnonals).
## Each unique element in B is to be estimated by the model.
sig.B <- c(b11 = 0.4, b22 = 0.6) # standard deviation of diagonal components (i.e. sqrt of the phylogenetic variance for each trait)
b12_rho = 0.5

Bcor <- matrix(c(c(1,b12_rho), # correlation matrix
                 c(b12_rho,1)),k,k, byrow = T) 

B <- matrix(kronecker(sig.B, sig.B),k,k)*Bcor # VCV (point-wise product). 
B
# N.B. Kronecker used here just for matrix formatting. Do not confuse with kronecker discussed in the text / tutorial.

## C is the residual trait-level VCV matrix (Sigma_{res} in the HTML).
sig.C <- c(c11 = 0.1, c22 = 0.1) # standard deviation of diagonal components (i.e. sqrt of the residual variance for each trait)
c12_rho = 0 # off-diagonal correlation coefficients

Ccor <- matrix(c(c(1,c12_rho), # correlation matrix
                 c(c12_rho,1)),k,k, byrow = T) 

C <- matrix(kronecker(sig.C, sig.C),k,k)*Ccor # VCV
C

# The phylogenetic taxon-level VCV matrix, A, is supplied. It is therefore treated as fixed and known without error, although may represent a transformation (see above).
A = A.mat

# number of species
n = nrow(A); n

# identity matrix
I = diag(n)

# simulate phylogenetic random effects for all traits as one draw from a MVN:
a <- mvrnorm(1,rep(0,n*k),kronecker(B,A)) # In the Kronecker, trait-level covariance captured by B, taxon-level covariance captured by A. 

# extract random effects for each trait from the vector, a.
a1 <- a[1:n]
a2 <- a[1:n + n]

# simulate residuals (on the link scale)
e <- mvrnorm(1,rep(0,n*k),kronecker(C,I))
e1 <- e[1:n]
e2 <- e[1:n + n]

# construct response traits from each linear predictor
g1 <- beta[1] + a1 + e1 # gaussian
g2 <- beta[2] + a2 + e2
b1 <- rbinom(n,1,plogis(beta[1] + a1 + e1)) # binomial
b2 <- rbinom(n,1,plogis(beta[2] + a2 + e2))

# generate df
species <- t.toy$tip.label
d.toy <- data.frame(species,g1,g2,b1,b2)
d.toy$animal <- d.toy$species # "animal" is a reserved term in MCMCglmm used to identify taxa in quantitative and phylogenetic models
d.toy$obs <- 1:nrow(d.toy)
head(d.toy)

## plot gaussians
# plot(d.toy$g1, d.toy$g2)
# table binomials
table(d.toy$b1, d.toy$b2)
table(d.toy$b1, d.toy$b2)/250 # prop

## take a look at probabilities on the inverse-link scale
# plogis(beta[1])
# plogis(beta[1] + a1 + e1)[1:10]

### END ####


##------------------- MULTIVARIATE GAUSSIAN ---------------------##

m.1 <- readRDS("m.1.rds")
b.1 <- readRDS("b.1.rds")

### FIT MODELS ####

### MCMCglmm
p <- list(R = list(R1=list(V = diag(2), nu=2)),
          G = list(G1=list(V = diag(2), nu=2)))
m.1 <- MCMCglmm(cbind(g1, g2) ~ trait-1, 
                random = ~us(trait):animal, 
                rcov = ~us(trait):units, 
                pedigree=t.toy,
                family = c("gaussian","gaussian"), 
                nodes="ALL", data = d.toy, prior=p,
                nitt=600000, burnin=100000, thin=500, 
                pr=TRUE,verbose = FALSE)

## DIAGNOSTICS
# view chain
# plot(m.1$VCV)
# calculate autocorrelation among samples
# autocorr(m.1$VCV)


## For good practice, let's explicitly state the default priors

### BRMS
b.1 <- brm(
  mvbind(g1, g2) ~ (1|p|gr(animal, cov = A)), 
  data = d.toy,
  data2 = list(A = A.mat),
  family = gaussian(),
  cores = 4,
  chains = 4, iter = 2000, thin = 2
)
b.1
## DIAGNOSTICS
pairs(b.1)
b.1 %>% plot
b.1 %>% pp_check(resp = "g1",ndraws = 100)
b.1 %>% pp_check(resp = "g2",ndraws = 100)
conditional_effects(b.1)


## COMPARISON BETWEEN METHODS
# N.B. Location estimates are directly comparable between MCMCglmm and BRMS (assuming no variances have fixed OR fixed at same level). However, for variance 
# components, MCMCglmm reports (co)variances while brms reports standard deviations and trait-level correlations. Therefore, some re-scaling is necessary to 
# to compare results between model fits. We have chosen to re-scale the estimates from MCMCglmm, as stdevs and corrs are more natural to interpret.
#
# Presumably, any differences are due to the differing priors.. do we want to confirm this or discuss merits of brms defaults?

# summary outputs
summary(m.1) # MCMCglmm
summary(b.1) # brms

# make nice plots to compare estimates and simulated values (i.e., reflected density plots)

# intercepts
summary(m.1)$solutions
summary(b.1)[["fixed"]]

# transform posterior columns first, then compute summary statistics and density plots etc.
# phylogenetic variances (elements of sig.B) and phylogenetic correlation (b12_rho)  # Q. SQRT OF CI CORRECT CI FOR VARIANCES?
sqrt(summary(m.1)$Gcovariances)[c(1,nrow(summary(m.1)$Gcovariances)),];data.frame(row.names = "traitg1:traitg2.animal", post.mean = mean(m.1$VCV[,"traitg1:traitg2.animal"]/sqrt(m.1$VCV[,"traitg1:traitg1.animal"]*m.1$VCV[,"traitg2:traitg2.animal"])))
summary(b.1)[["random"]]

# residual covariances (elements of sig.C) and residual correlation (c12_rho) 
sqrt(summary(m.1)$Rcovariances[c(1,nrow(summary(m.1)$Gcovariances)),]);data.frame(row.names = "traitg1:traitg2.units", post.mean = mean(m.1$VCV[,"traitg1:traitg2.units"]/sqrt(m.1$VCV[,"traitg1:traitg1.units"]*m.1$VCV[,"traitg2:traitg2.units"])))
summary(b.1)[["spec_pars"]];summary(b.1)[["rescor_pars"]]


### END ####
