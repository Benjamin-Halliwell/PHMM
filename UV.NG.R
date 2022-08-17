library("MCMCglmm");library("brms");library("tidyverse");library("ape");
library("phytools");library("MASS");library("plyr");library("bindata");

##------------- UNIVARIATE NON-GAUSSIAN (BERNOULLI) -------------##


#### DEFINE SIMULATION PARAMETERS ####


t.toy <- geiger::sim.bdtree(b=1, d=0, stop="taxa", n=250, extinct=FALSE) # may be best not to use pure birth trees (Adams and Collyer 2018)
t.toy <- multi2di(t.toy)
plot(t.toy)
# save(t.toy, file="t.toy")

# Create VCV matrix from tree, scale to correlation matrix for BRMS
A.mat <- ape::vcv.phylo(t.toy, corr = T) # scale with corr = T

# ## alternative method 
# inv.phylo <- MCMCglmm::inverseA(t.toy, nodes = "TIPS", scale = TRUE)
# A <- solve(inv.phylo$Ainv)
# rownames(A) <- rownames(inv.phylo$Ainv)
# A.mat <- A

# number of traits
k = 2

# intercepts for each trait - HOW TO INTERPRET FOR BINOMIAL DATA?
beta = c(1,2)

## B is the phylogenetic trait-level VCV matrix (Sigma_{phy} in the HTML). B specifies the phylogenetic
## variance in each trait (diagonals) as well as the phylogenetic covariance between traits (off-diagnonals).
## Each unique element in B is to be estimated by the model.
sig.B <- c(b11 = 2, b22 = 4) # standard deviation of diagonal components (i.e. sqrt of the phylogenetic variance for each trait)
b12_rho = 0

Bcor <- matrix(c(c(1,b12_rho), # correlation matrix
                 c(b12_rho,1)),k,k, byrow = T) 

B <- matrix(kronecker(sig.B, sig.B),k,k)*Bcor # VCV (point-wise product). 
# N.B. Kronecker used here just for matrix formatting. Do not confuse with kronecker discussed in the text / tutorial.


## C is the residual trait-level VCV matrix (Sigma_{res} in the HTML).
sig.C <- c(c11 = 1, c22 = 2) # standard deviation of diagonal components (i.e. sqrt of the residual variance for each trait)
c12_rho = 0.5 # off-diagonal correlation coefficients


Ccor <- matrix(c(c(1,c12_rho), # correlation matrix
                 c(c12_rho,1)),k,k, byrow = T) 

C <- matrix(kronecker(sig.C, sig.C),k,k)*Ccor # VCV


# The phylogenetic taxon-level VCV matrix, A, is supplied. It is therefore treated as fixed and known without error, although may represent a transformation (see above).
A = A.mat

# number of species
n = nrow(A)

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
y1 <- rbinom(n,1,plogis(beta[1] + a1 + e1))
y2 <- rbinom(n,1,plogis(beta[2] + a2 + e2))

# generate df
species <- t.toy$tip.label
d.toy <- data.frame(species,y1,y2)
d.toy$animal <- d.toy$species # "animal" is a reserved term in MCMCglmm used to identify taxa in quantitative and phylogenetic models
d.toy$obs <- 1:nrow(d.toy)

table(d.toy$y1, d.toy$y2)


### END ####


### FIT MODLES ####

### MCMCglmm 

## PRIORS
# p <- list(G = list(G1=list(V=1, nu=0.002)), # Inv Wish for random intercepts (MCMCglmm default)
#           R = list(V = 1, fix = 1))
p <- list(G = list(G1=list(V=1, nu=1000, alpha.mu=0, alpha.V=1)), # parameter expanded G priors increase ESS of intercept and phylogenetic variance
          R = list(V = 1, fix = 1))
# p <- list(G = list(G1=list(V=1, nu=1, alpha.mu=0, alpha.V=1)), # reducing nu reduces ESS and gives inflated estimates of phylogenetic variance
#           R = list(V = 1, fix = 1))
# p <-  list(B = list(mu = c(0, 0), V = diag(2) * (1 + pi^2/3)), # prior to use when incomplete separation is a problem for fixed effects
#                    R = list(V = 1, fix = 1))
# p <- list(G = list(G1=list(V=1, nu=1, alpha.mu=0, alpha.V=1000)), # parameter expanded prior recommended for binary data in MCMCMglmm course notes reduces ESS and gives inflated estimates of phylogenetic variance
#           R = list(V = 1,  fix = 1))


# intercept only model
m.1<-MCMCglmm(y1 ~ 1,
              random = ~animal,
              pedigree=t.toy, 
              family = c("categorical"), 
              data = d.toy, prior = p,
              nitt=350000, burnin=100000, thin=300,
              pr=T, verbose = F)
summary(m.1)

# y2 included as fixed effect
m.2<-MCMCglmm(y1 ~ y2,
              random = ~animal,
              pedigree=t.toy, 
              family = c("categorical"), 
              data = d.toy, prior = p,
              nitt=350000, burnin=100000, thin=300,
              pr=T, verbose = F)
summary(m.2)

## DIAGNOSTICS
# view chain
plot(m.1$VCV)
plot(m.2$VCV)

# calculate autocorrelation among samples. All values above lag 0 should be close to 0
autocorr(m.1$VCV)
autocorr(m.2$VCV)


## RESCALE BASED ON RESIDUAL VARIANCE ? ##

# For binomial, family variance only makes sense in context of the residual variance.
# Use intra-class correlation coefficient to rescale the phylogenetic variance 
# according to the assumed residual variance (fixed at 1)
Int.1 <- m.1$VCV[, 1]/(rowSums(m.1$VCV) + pi^2/3) # rowSums used to add res var
mean(Int.1);HPDinterval(Int.1)


# rescale estimates by the estimated residual variance in order to obtain the
# posterior distributions of the parameters under the assumption that the 
# actual residual variance is equal to some other value.
c2 <- ((16 * sqrt(3))/(15 * pi))^2 # constant for logit link, c = 1 for probit
Int.2 <- m.1$Sol/sqrt(1 + c2 * m.1$VCV[, 2])
mean(Int.2[,"(Intercept)"]);HPDinterval(Int.2[,"(Intercept)"])


#---------------------------------------------------------------------#

## BRMS

# see which priors brms is using by default
priors_b.1<-get_prior(y1 ~ 1 + (1 | gr(animal, cov = A)),
                        data = d.toy,
                        data2 = list(A = A.mat),
                        family = bernoulli())  

b.1 <- brm(
  y1 ~ 1 + (1|gr(animal, cov = A)), 
  data = d.toy, 
  family = bernoulli(),
  control = list(adapt_delta = 0.9),
  data2 = list(A = A.mat),
  cores=4,
  chains=4, iter = 6000, thin = 3
)
b.1

b.2 <- brm(
  y1 ~ y2 + (1|gr(animal, cov = A)), 
  data = d.toy, 
  family = bernoulli(),
  control = list(adapt_delta = 0.9),
  data2 = list(A = A.mat),
  cores=4,
  chains=4, iter = 6000, thin = 3
)
b.2


## DIAGNOSTICS
b.1 %>% plot
b.1 %>% pp_check(nsamples = 100)
pairs(b.1)

b.2 %>% plot
b.2 %>% pp_check(nsamples = 100)
pairs(b.2) # sd_animal and b_y2 positively correlated

# If divergent transitions error thrown. Explain adapt_delta and re-fit.


## COMPARISON BETWEEN METHODS

# summary outputs
summary(m.1) # MCMCglmm
summary(b.1) # brms

## Parameter estimates:
# N.B. Intercept estimates are directly comparable between MCMCglmm and BRMS but the variance components reported
# by MCMCglmm need to be sqrt() to compare to the SDs reported by BRMS (as well as the generating parameters for our sims).

# intercept
summary(m.1)$solutions
summary(b.1)[["fixed"]]

summary(m.2)$solutions
summary(b.2)[["fixed"]]

# phylogenetic variance (sig.B)
sqrt(summary(m.1)$Gcovariances) # N.B. effective sample size is now incorrect
summary(b.1)[["random"]]

sqrt(summary(m.2)$Gcovariances) 
summary(b.2)[["random"]]

# residual variance (sig.C)
sqrt(summary(m.1)$Rcovariances)
summary(b.1)[["spec_pars"]]

sqrt(summary(m.2)$Rcovariances)
summary(b.2)[["spec_pars"]]
sqrt(2.6)


## NOTES
#
# Both methods provide similar point estimates and CIs
#
# Both methods tend to overestimate phylo effects but this could be due to:
#   1) Covariances with other traits involved in the sim but not included in the model (i.e. un-observed variables)
#   2) Trade off between sig.B and sig.C
#   3) Relatively high residual variance (values of sig.C) in a trait?


### END ####

