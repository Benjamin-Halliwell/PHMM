library("MCMCglmm");library("brms");library("tidyverse");library("ape");
library("phytools");library("MASS");library("plyr");library("bindata");


#### SIMULATE DATA AND TREE ####

## simulate tree
t.toy <- geiger::sim.bdtree(b=1, d=0, stop="taxa", n=200, extinct=FALSE) # may be best not to use pure birth trees (Adams and Collyer 2018)
t.toy <- multi2di(t.toy)
plot(t.toy, cex=0.4, label.offset = 0.05)
# save(t.toy, file="t.toy")

# Create VCV matrix from tree
A.mat <- ape::vcv.phylo(t.toy) # BM


### DEFINE SIMULATION PARAMETERS 

# number of traits
k = 2

# intercepts for each trait
beta = c(1,2)

## B is the phylogenetic trait-level VCV matrix (Sigma_{phy} in the HTML). B specifies the phylogenetic
## variance in each trait (diagonals) as well as the phylogenetic covariance between traits (off-diagnonals).
## Each unique element in B is to be estimated by the model.
sig.B <- c(b11 = 2, b22 = 4) # standard deviation of diagonal components (i.e. sqrt of the phylogenetic variance for each trait)
b12_rho = 0.75

Bcor <- matrix(c(c(1,b12_rho), # correlation matrix
                 c(b12_rho,1)),k,k, byrow = T) 

B <- matrix(kronecker(sig.B, sig.B),k,k)*Bcor # VCV (point-wise product). 
# N.B. Kronecker used here just for matrix formatting. Do not confuse with kronecker discussed in the text / tutorial.


## C is the residual trait-level VCV matrix (Sigma_{res} in the HTML).
sig.C <- c(c11 = 1, c22 = 2) # standard deviation of diagonal components (i.e. sqrt of the residual variance for each trait)
c12_rho = 0.25# off-diagonal correlation coefficients


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
y1 <- beta[1] + a1 + e1 # Gaussian variables constructed directly
y2 <- beta[2] + a2 + e2

# generate df
species <- t.toy$tip.label
d.toy <- data.frame(species,y1,y2)
d.toy$animal <- d.toy$species # "animal" is a reserved term in MCMCglmm used to identify taxa in quantitative and phylogenetic models
d.toy$obs <- 1:nrow(d.toy)

str(d.toy)

### END ####



p <- list(R = list(R1=list(V = diag(2), nu=2)),
       G = list(G1=list(V = diag(2), nu=2)))
m.1 <- MCMCglmm(cbind(y1, y2) ~ trait:parity-1, 
                    random = ~us(trait):animal, 
                    rcov = ~us(trait):units, 
                    pedigree=tree.u,
                    family = c("gaussian","gaussian"), 
                    nodes="ALL", data = d, prior=p,
                    nitt=350000, burnin=100000, thin=300, pr=TRUE,verbose = FALSE) 
summary(m.1)


# BRMS
b.1 <- brm(
  mvbind(y1, y2) ~ (1|p|gr(animal, cov = A)), 
  data = d.toy, 
  control = list(adapt_delta = 0.9),
  family = gaussian(), 
  data2 = list(A = A.mat),
  cores = 4,
  chains = 4, iter = 6000, thin = 3
)
summary(b.21)
pairs(b.1)
b.1 %>% plot
b.1 %>% pp_check(resp = "y1",ndraws = 100)
b.1 %>% pp_check(resp = "y2",ndraws = 100)

# saveRDS(b.1, "PMM_brms_2_VAR.rds")







### IS PGLS EQUIVALENT TO PARTIALLING FOR PHYLOGENY? ###

library(nlme); library(phytools)

# Consider the model,
#
#   y = XB + Zb + e
#
# where XB represents a vector of fixed effects including an intercept = 1, and  Z represents the design matrix for the random effects b where,
#
#   b ~ N(0, sigma^2_b*I_g))
#   e ~ N(0, sigma^2*I)
#

### TREE ###

# generate tree from overall covariance of errors, Zb + e (random + residual)
n = 10 # number of species in tree
k = 2 # group levels (number of clades)

sigma = 1 # residual variance
sigma_b = 3 # variance of group (phylo) random effects

I = diag(n) # obs level identity
I_g = diag(k) # group level identity

Z = matrix(c(rep(c(1,c(rep(0,k-1))),n/k),rep(c(0,c(rep(1,k-1))),n/k)),n,k, byrow = T) # design matrix for group level (phylo) random effects
Sigma = sigma_b * I_g # group level covariance structure

A = Z%*%Sigma%*%t(Z) + sigma*I # overall covariance matrix (Zb + e)
dimnames(A) <- list(LETTERS[1:n],LETTERS[1:n])
A # covariance matrix capturing group and residual error

# create tree from VCV and plot
tree <- vcv2phylo(A)

# simple tree with two polytomous clades separated by deep time. This tree structure allows us to
# represent trait change owing to shared ancestry as a simple vector to be partialled (see below)
plot(tree)


### TRAITS ###

# create trait x1 by simulating BM on tree
x1 <- fastBM(tree, a = 0 , sig2 = 1)

# extract the extent of change in x1 along shared branches for each clade. This is equivalent to the value of x1 at the 
# crown node of each clade, which under BM is simply the mean of x1 for each clade. Because y is directly related to x1
# (see formula for y below), x2 represents the phylogenetically conservative portion of the relationship 
# between y and x1. Partialling out x2 = partialling out phylogeny
x2 <-  ifelse(names(x1) %in% LETTERS[1:(n/2)], mean(x1[LETTERS[1:(n/2)]]), mean(x1[LETTERS[(n/2):n]]))

# x3 represents observation level (residual) error in x1, i.e. the change in x1 along terminal branches. This is the variation in x1 that 
# is related to y but not to phylogeny and should therefore be equal to the remainder after partialling out x2, e.g. x3 = res(lm(x1 ~ x2))
x3 <- x1 - x2 

## construct response variable y from intercept and x1. Add some residual error.
y = 1 + 1*x1 + rnorm(n,0,0.1)

# strong relationship between y and x1
plot(y,x1)

# other quantities to note
y_cor = y - x2 # y corrected for changes in x1 due to shared ancestry (i.e. corrected for x2)
y - 1 - x2 - x3 # residual error in y (i.e., y - intercept - phylo correlation with x1 - residual correlation with x1)
y - y_cor # the part of y that is phylogenetically determined via association with x1
round((y - y_cor),3) == round(x2, 3) # equivalent to x2


### MODELS ###

dat <- data.frame(y=y,x1=x1,x2=x2)

# fit lm with both x1 and x2 as predictors
fit.lm <- lm(y ~ x1 + x2, dat); summary(fit.lm); cor.test(x1, x2)

# partial for x1
x1.lm <- lm(x1 ~ x2, dat) # 1. regress x1 on x2
x1.par <- x1.lm$residuals # 2. extract residuals

# x1.par = x3 for clade A-E but only approximate for clade F-J, WHY? - !! RESOLVE !!
round(x3,2) == round(x1.par,2)
round(x3,2);round(x1.par,2)

fit.lm.partial <- lm(y ~ x1.par, dat); summary(fit.lm.partial) # 3. re-fit using residuals of x1.par as predictor

# fit PGLS with x1 as only predictor
f1 <- gls(y ~ x1, dat, corBrownian(1, tree), method='ML'); summary(f1)
# f1 <- gls(y ~ x1, dat, corPagel(1,tree,fixed=T), method='ML'); summary(f1) # fix / estimate lambda instead

# compare coefficient estimates between PGLS and LM partialling out x2
round(fit.lm.partial$coefficients,4);round(f1$coefficients,4)

# Beta1 coefficients very close between the two methods, indicating PGLS is indeed equivalent to 'partialling for phylogeny'. 


