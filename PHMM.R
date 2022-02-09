library("MCMCglmm");library("brms");library("tidyverse");library("ape");
library("phytools");library("MASS");library("plyr");library("bindata");


load(file="t.toy")
load(file="d.toy")


## QUESTIONS ##

# Are we able to define interactions between response variables (y1 and y2) on
# a third response variable (y3) when using conditioning to predict?


#### SIMULATE DATA AND TREE ####

## simulate tree
t.toy <- geiger::sim.bdtree(b=1, d=0, stop="taxa", n=100, extinct=FALSE) # may be best not to use pure birth trees (Adams and Collyer 2018)
t.toy <- multi2di(t.toy)
plot(t.toy)
# save(t.toy, file="t.toy")

# Create VCV matrix from tree
A.mat <- ape::vcv.phylo(t.toy) # BM

# transform A.mat according to different models of evolution:
A.LA <- "lambda"
A.WN <- "white noise"
A.OU <- "OU"


## DEFINE SIMULATION PARAMETERS ##

# number of traits
k = 5
# intercepts for each trait
beta = c(1,2,1,-2,-1)


## B is the phylogenetic trait-level VCV matrix (Sigma_{phy} in the HTML). B specifies the phylogenetic
## variance in each trait (diagonals) as well as the phylogenetic covariance between traits (off-diagnonals).
## The value of each element in B is to be estimated by the model.
sig.B <- c(b11 = 1, b22 = 2, b33 = 1, b44 = 1, b55 = 1) # standard deviation of diagonal components (i.e. sqrt of the phylogenetic variance for each trait)
b12_rho = 0.75; b13_rho = 0.25; b14_rho = 0.25; b15_rho = 0.25 # off-diagonal correlation coefficients
b23_rho = 0.25; b24_rho = 0.25; b25_rho = 0.25; 
b34_rho = 0.25; b35_rho = 0.25; b45_rho = 0.25

Bcor <- matrix(c(c(1,b12_rho,b13_rho,b14_rho,b15_rho), # correlation matrix
                 c(b12_rho,1,b23_rho,b24_rho,b25_rho),
                 c(b13_rho,b23_rho,1,b34_rho,b35_rho),
                 c(b14_rho,b24_rho,b34_rho,1,b45_rho),
                 c(b15_rho,b25_rho,b35_rho,b45_rho,1)),k,k, byrow = T) 

B <- matrix(kronecker(sig.B, sig.B),k,k)*Bcor # VCV (point-wise product)


## C is the residual trait-level VCV matrix (Sigma_{res} in the HTML).
sig.C <- c(c11 = 2, c22 = 1, c33 = 1, c44 = 1, c55 = 1) # standard deviation of diagonal components (i.e. sqrt of the residual variance for each trait)
c12_rho = 0.25; c13_rho = 0.25; c14_rho = 0.25; c15_rho = 0.25 # off-diagonal correlation coefficients
c23_rho = 0.25; c24_rho = 0.25; c25_rho = 0.25; 
c34_rho = 0.25; c35_rho = 0.25; c45_rho = 0.75

Ccor <- matrix(c(c(1,c12_rho,c13_rho,c14_rho,c15_rho), # correlation matrix
                 c(c12_rho,1,c23_rho,c24_rho,c25_rho),
                 c(c13_rho,c23_rho,1,c34_rho,c35_rho),
                 c(c14_rho,c24_rho,c34_rho,1,c45_rho),
                 c(c15_rho,c25_rho,c35_rho,c45_rho,1)),k,k, byrow = T) 

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
a3 <- a[1:n + 2*n]
a4 <- a[1:n + 3*n]
a5 <- a[1:n + 4*n]

# simulate residuals (on the link scale)
e <- mvrnorm(1,rep(0,n*k),kronecker(C,I))
e1 <- e[1:n]
e2 <- e[1:n + n]
e3 <- e[1:n + 2*n]
e4 <- e[1:n + 3*n]
e5 <- e[1:n + 4*n]

# construct response traits from each linear predictor
y1 <- beta[1] + a1 + e1 # Gaussian variables constructed directly
y2 <- beta[2] + a2 + e2
y3 <- rpois(n,exp(beta[3] + a3 + e3)) # non-Gaussian variables constructed using inverse link functions
y4 <- rbinom(n,1,plogis(beta[4] + a4 + e4))
y5 <- rbinom(n,1,plogis(beta[5] + a5 + e5))

# generate df
species <- t.toy$tip.label
d.toy <- data.frame(species,y1,y2,y3,y4,y5)
d.toy$animal <- d.toy$species # "animal" is a reserved term in MCMCglmm used to identify taxa in quantitative and phylo- genetic models
d.toy$y5 <- as.factor(y5) # second Bernoulli variable as factor in order to estimate group level VCV matrices ("at.level"/"by" coding).
d.toy$y5 <- revalue(d.toy$y5, c("0"="L1", "1"="L2"))
d.toy$obs <- 1:nrow(d.toy)
# save(d.toy, file="d.toy")
head(d.toy)
summary(d.toy)

table(d.toy$y4, d.toy$y5)

#### END ####


#### PLOTTING ####

# plot data
plot(y1~y2, d.toy)
plot(y1~y3, d.toy)

# phylo dist of y1
names(y1) <- t.toy$tip.label
plotTree.barplot(t.toy, y1)
# phylo dist of y2
names(y2) <- t.toy$tip.label
plotTree.barplot(t.toy, y2)
# phylo dist of y3
names(y3) <- t.toy$tip.label
par(mfrow=c(1,1), mar=c(2,2,2,2))
dotTree(t.toy, y3)

# strength of trait-level correlation is reduced by taxon-level and residual covariance
par(mfrow=c(1,1), mar=c(5,5,2,2))
plot(y1,y2); cor(y1, y2)
# If instead we estimate trait level correlation across multiple draws for the same taxon (i.e. taxon is fixed):
traitT <- mvrnorm(100,rep(0,n*nrow(B)),kronecker(B,A))
# 100 samples of both y1 and y2 for a given taxon reveals a phylogenetic trait-level 
# correlation closer to the generating value of 0.75:
plot(traitT[,n+1],traitT[,1]) # column 1 contains samples of t1 for s1, column n+1 contains samples of t2 for s1, etc.
cor(traitT[,n+1],traitT[,1]) # estimated cor is close to supplied value: b12_rho = 0.75

#### END ####


#### PHYLOGENETIC SIGNAL ####

# HYP: y1 should have less signal than y2 (b11 = 1, b22 = 2)

# lambda more sensitive to phylogenetic variance in traits, prone to type 1 error?
# k only significant at relatively high levels of phylogenetic variance, prone to type 2 error?

# Ives 2018 shows that p-values from this test are unreliable (need to bootstrap)

phylosig(t.toy, x = d.toy$y1, method= "lambda", test=T, nsim=10000)
phylosig(t.toy, x = d.toy$y1, method= "K", test=T, nsim=10000)

phylosig(t.toy, x = d.toy$y2, method= "lambda", test=T, nsim=10000)
phylosig(t.toy, x = d.toy$y2, method= "K", test=T, nsim=10000)


#### END ####


#### MODELS ####


##------------- UNIVARIATE GAUSSIAN -------------##

## MCMCglmm 
p <- list(G = list(G1 = list(V = 10, nu = 0.002)), 
          R = list(V = 10, nu = 0.002))

m.1<-MCMCglmm(y1 ~ 1,
                random = ~animal,
                pedigree=t.toy, 
                family = c("gaussian"), 
                data = d.toy, prior = p,
                nitt=250000, burnin=50000, thin=200,
                pr=TRUE,verbose = FALSE)
summary(m.1)

## DIAGNOSTICS
par(mar=(c(2,2,2,2)))
plot(m.1$Sol)
plot(m.1$VCV)
# calculate autocorrelation among samples. All values above lag 0 should be close to 0
autocorr(m.1$VCV)

# EXTRACT OUTPUTS
# calculate different statistics from the posterior dist.
mean(m.1$VCV[,"animal"])
posterior.mode(m.1$VCV)
median(m.1$VCV[,"animal"])
# credible intervals
HPDinterval(m.1$VCV)

#--- ~ ---#


## BRMS
b.1 <- brm(
  y1 ~ 1 + (1|gr(animal, cov = A)), 
  data = d.toy, 
  family = gaussian(), 
  data2 = list(A = A.mat),
  cores=4,
  chains=4, iter = 6000, thin = 3
)
b.1

## DIAGNOSTICS
b.1 %>% plot
b.1 %>% pp_check(nsamples = 100)
pairs(b.1) # sd_animal and sigma are trading off

# divergent transitions error thrown. Explain adapt_delta and re-fit.


## COMPARISON BETWEEN METHODS

# summary outputs
summary(m.1) # MCMCglmm
summary(b.1) # brms

# N.B. Intercept estimates are directly comparable between MCMCglmm and BRMS but the variance components reported
# by MCMCglmm need to be sqrt() to compare to the SDs reported by BRMS (as well as the generating parameters for our sims).

## Parameter estimates:
# intercept
summary(m.1)$solutions
summary(b.1)[["fixed"]]
# phylogenetic variance (sig.B)
sqrt(summary(m.1)$Gcovariances) # N.B. effective sample size is now incorrect
summary(b.1)[["random"]]
# residual variance (sig.C)
sqrt(summary(m.1)$Rcovariances)
summary(b.1)[["spec_pars"]]


## NOTES
#
# Both methods provide similar point estimates and CIs
#
# Both methods tend to overestimate phylo effects but this could be due to:
#   1) Covariances with other traits involved in the sim but not included in the model (i.e. un-observed variables)
#   2) Trade off between sig.B and sig.C
#   3) Relatively high residual variance (values of sig.C) in a trait?



##---------- MULTIVARIATE GAUSSIAN -----------##

## MCMCglmm
p2 <- list(G = list(G1 = list(V = diag(10,2), nu = 1.002)), 
          R = list(V = diag(10,2), nu = 1.002))

m.2<-MCMCglmm(cbind(y1, y2) ~ trait-1,
                random = ~us(trait):animal,
                rcov = ~us(trait):units,
                pedigree=t.toy,
                family = c("gaussian","gaussian"), 
                data = d.toy, prior=p2, 
                nitt=250000, burnin=50000, thin=200,
                pr=TRUE,verbose = FALSE)
summary(m.2)
plot(m.2$VCV)
autocorr(m.2$VCV)[,,1]


# Various statistics can be calculated as ratios of the variance explained in each trait, 
# or that of the entire model, by phylogenetic and residual components: 

# PHYLOGENETIC HERITABILITY (sig.B^2/(sig.B^2 + sig.C^2)
posterior.mode(m.2$VCV[,'traity1:traity1.animal']/(m.2$VCV[,'traity1:traity1.animal']+m.2$VCV[,'traity1:traity1.units'])) # y1
posterior.mode(m.2$VCV[,'traity2:traity2.animal']/(m.2$VCV[,'traity2:traity2.animal']+m.2$VCV[,'traity2:traity2.units'])) # y2

# A KIND OF PARTIAL (PHYLO) R^2?
posterior.mode((m.2$VCV[,'traity1:traity1.animal']+m.2$VCV[,'traity2:traity2.animal'])/
                 (m.2$VCV[,'traity1:traity1.animal']+m.2$VCV[,'traity1:traity1.units']+
                    m.2$VCV[,'traity2:traity2.animal']+m.2$VCV[,'traity2:traity2.units']))


## Calculate (co)variances directly from the posterior of each parameter

# PHYLOGENETIC VARIANCE (sig.B)
paste0(round(sqrt(mean(m.2$VCV[,"traity1:traity1.animal"])),2)," (",round(sqrt(HPDinterval(m.2$VCV[,"traity1:traity1.animal"]))[1],2),", ",round(sqrt(HPDinterval(m.2$VCV[,"traity1:traity1.animal"]))[2],2),")")
paste0(round(sqrt(mean(m.2$VCV[,"traity2:traity2.animal"])),2)," (",round(sqrt(HPDinterval(m.2$VCV[,"traity2:traity2.animal"]))[1],2),", ",round(sqrt(HPDinterval(m.2$VCV[,"traity2:traity2.animal"]))[2],2),")")

# RESIDUAL VARIANCE (sig.C)
paste0(round(sqrt(mean(m.2$VCV[,"traity1:traity1.units"])),2)," (",round(sqrt(HPDinterval(m.2$VCV[,"traity1:traity1.units"]))[1],2),", ",round(sqrt(HPDinterval(m.2$VCV[,"traity1:traity1.units"]))[2],2),")")
paste0(round(sqrt(mean(m.2$VCV[,"traity2:traity2.units"])),2)," (",round(sqrt(HPDinterval(m.2$VCV[,"traity2:traity2.units"]))[1],2),", ",round(sqrt(HPDinterval(m.2$VCV[,"traity2:traity2.units"]))[2],2),")")

# PHYLOGENETIC CORRELATION (b_{ij}_rho)
posterior.mode(m.2$VCV[,'traity1:traity2.animal']/sqrt((m.2$VCV[,'traity1:traity1.animal']*m.2$VCV[,'traity2:traity2.animal'])))
HPDinterval(m.2$VCV[,'traity1:traity2.animal']/sqrt((m.2$VCV[,'traity1:traity1.animal']*m.2$VCV[,'traity2:traity2.animal'])))

# RESIDUAL CORRELATION (c_{ij}_rho)
posterior.mode(m.2$VCV[,'traity1:traity2.units']/sqrt((m.2$VCV[,'traity1:traity1.units']*m.2$VCV[,'traity2:traity2.units'])))
HPDinterval(m.2$VCV[,'traity1:traity2.units']/sqrt((m.2$VCV[,'traity1:traity1.units']*m.2$VCV[,'traity2:traity2.units'])))

#--- ~ ---#


## BRMS
b.2 <- brm(
  mvbind(y1, y2) ~ (1|p|gr(animal, cov = A)), 
  data = d.toy, 
  control = list(adapt_delta = 0.9),
  family = gaussian(), 
  data2 = list(A = A.mat),
  cores = 4,
  chains = 4, iter = 6000, thin = 3
)
summary(b.2)
pairs(b.2)
b.2 %>% plot
b.2 %>% pp_check(resp = "y1",ndraws = 100)
b.2 %>% pp_check(resp = "y2",ndraws = 100)

# saveRDS(b.2, "PMM_brms.rds")


#--- ~ ---#


## COMPARE MODEL OUTPUTS

# summaries
summary(m.2)
summary(b.2)

## Extract parameter estimates
# intercepts (beta[i])
summary(m.2)$solutions
summary(b.2)$fixed
# phylogenetic (co)variances (sig.B)
summary(m.2)$Gcovariances[c(1,4),]
m.2.Gcov <- c(posterior.mode(m.2$VCV[,'traity1:traity2.animal']/sqrt((m.2$VCV[,'traity1:traity1.animal']*m.2$VCV[,'traity2:traity2.animal']))),HPDinterval(m.2$VCV[,'traity1:traity2.animal']/sqrt((m.2$VCV[,'traity1:traity1.animal']*m.2$VCV[,'traity2:traity2.animal'])))[,]); names(m.2.Gcov)[1] <- "y1:y2.animal"; m.2.Gcov
summary(b.2)$random
# residual (co)variances (sig.C)
sqrt(summary(m.2)$Rcovariances)[c(1,4),]
m.2.Gcov <- c(posterior.mode(m.2$VCV[,'traity1:traity2.units']/sqrt((m.2$VCV[,'traity1:traity1.units']*m.2$VCV[,'traity2:traity2.units']))),HPDinterval(m.2$VCV[,'traity1:traity2.units']/sqrt((m.2$VCV[,'traity1:traity1.units']*m.2$VCV[,'traity2:traity2.units'])))[,]); names(m.2.Gcov)[1] <- "y1:y2.units"; m.2.Gcov
summary(b.2)$spec_pars;summary(b.2)$rescor_pars


## NOTES
#
# MCMCglmm does better job at estimating residual covariance (c_{ij}_rho)
#



##---------- MULTIVARIATE NON-GAUSSIAN -----------##

## MCMCglmm
p3 <- list(G = list(G1 = list(V = diag(1, 4), nu = 3.002)), # NOTE - NEED TO USE DIFFERENT PRIORS FOR NON-GAUSSIAN RESPONSES
           R = list(V = diag(1, 4), nu = 3.002))

m.3 <- MCMCglmm(cbind(y1, y2, y3, y4) ~ trait-1,
                random = ~us(trait):animal,
                rcov = ~us(trait):units,
                pedigree = t.toy,
                family = c("gaussian","gaussian","poisson", "categorical"), 
                data = d.toy, prior=p3, 
                nitt=250000, burnin=50000, thin=200,
                pr=TRUE,verbose=FALSE)
summary(m.3)

#--- ~ ---#

## BRMS
# for MV non_gaussian models, we can use the bf() function in brms to
# define separate formulas and error terms for each response variable.
# In this example, we set sigma to a nominal value of 0.01 for our 
# Gaussian variables as a trick to circumvent XYZ. It is possible to 
# avoid this by editing the STAN code under the hood, but for the sake
# of simplicity we avoid that here.

bf_y1 <- bf(y1 ~ 1 + (1|p|gr(animal, cov = A)) + (1|q|obs), sigma = 0.01) + gaussian()
bf_y2 <- bf(y2 ~ 1 + (1|p|gr(animal, cov = A)) + (1|q|obs), sigma = 0.01) + gaussian()
bf_y3 <- bf(y3 ~ 1 + (1|p|gr(animal, cov = A)) + (1|q|obs)) + poisson()
bf_y4 <- bf(y4 ~ 1 + (1|p|gr(animal, cov = A)) + (1|q|obs)) + bernoulli()

b.3 <- brm(bf_y1 + bf_y2 + bf_y3 + bf_y4 + set_rescor(FALSE),
           data = d.toy, 
           family = gaussian(), 
           data2 = list(A = A.mat),
           cores = 4,
           chains = 4, iter = 6000, thin = 3
)
b.3

## COMPARE MODEL OUTPUTS

# summaries
summary(m.3)
summary(b.3)

## Extract parameter estimates
# intercepts (beta[i])
summary(m.3)$solutions
summary(b.3)$fixed

# phylogenetic (co)variances (sig.B)
sqrt(summary(m.2)$Gcovariances)[c(1,4),]; m.2.Gcov <- c(posterior.mode(m.2$VCV[,'traity1:traity2.animal']/sqrt((m.2$VCV[,'traity1:traity1.animal']*m.2$VCV[,'traity2:traity2.animal']))),HPDinterval(m.2$VCV[,'traity1:traity2.animal']/sqrt((m.2$VCV[,'traity1:traity1.animal']*m.2$VCV[,'traity2:traity2.animal'])))[,]); names(m.2.Gcov)[1] <- "y1:y2.animal"; m.2.Gcov
summary(b.2)$random

# residual (co)variances (sig.C)
sqrt(summary(m.2)$Rcovariances)[c(1,4),]; m.2.Rcov <- c(posterior.mode(m.2$VCV[,'traity1:traity2.units']/sqrt((m.2$VCV[,'traity1:traity1.units']*m.2$VCV[,'traity2:traity2.units']))),HPDinterval(m.2$VCV[,'traity1:traity2.units']/sqrt((m.2$VCV[,'traity1:traity1.units']*m.2$VCV[,'traity2:traity2.units'])))[,]); names(m.2.Rcov)[1] <- "y1:y2.units"; m.2.Rcov
summary(b.2)$spec_pars;summary(b.2)$rescor_pars




##---------- MULTIVARIATE GAUSSIAN, MULTIPLE GROUPING FACTORS -----------##

# grouping factors much easier to specify in brms: 'by=' compared with 'at.level' coding in MCMCglmm


## MCMCglmm
p4=list(G = list(G1=list(V = diag(2), nu = 1.002),
                 G2=list(V = diag(2), nu = 1.002)),
        R = list(R1=list(V = diag(2), nu = 1.002),
                 R2=list(V = diag(2), nu = 1.002)))

m.4<-MCMCglmm(cbind(y1, y2) ~ trait:y5-1,
                random = ~us(at.level(y5,'L1'):trait):animal + us(at.level(y5,'L2'):trait):animal,
                rcov = ~us(at.level(y5,'L1'):trait):units + us(at.level(y5,'L2'):trait):units,
                family = c("gaussian","gaussian"), 
                pedigree=t.toy, data = d.toy, prior=p4, 
                nitt=210000, burnin=10000, thin=200,
                pr=TRUE,verbose = FALSE)
summary(m.4)
# plot(m.4$VCV)
autocorr(m.4$VCV)[,1:3,1]

mean(m.4$VCV[,'at.level(y5, "L1"):traity1:at.level(y5, "L1"):traity2.animal']) # negative 

# PHYLOGENETIC CORRELATION #
# L1
posterior.mode(m.4$VCV[,'at.level(y5, "L1"):traity1:at.level(y5, "L1"):traity2.animal']/sqrt((m.4$VCV[,'at.level(y5, "L1"):traity1:at.level(y5, "L1"):traity1.animal']*m.4$VCV[,'at.level(y5, "L1"):traity2:at.level(y5, "L1"):traity2.animal'])))
round(HPDinterval(m.4$VCV[,'at.level(y5, "L1"):traity1:at.level(y5, "L1"):traity2.animal']/sqrt((m.4$VCV[,'at.level(y5, "L1"):traity1:at.level(y5, "L1"):traity1.animal']*m.4$VCV[,'at.level(y5, "L1"):traity2:at.level(y5, "L1"):traity2.animal']))),2)[1:2]
# L2
posterior.mode(m.4$VCV[,'at.level(y5, "L2"):traity1:at.level(y5, "L2"):traity2.animal']/sqrt((m.4$VCV[,'at.level(y5, "L2"):traity1:at.level(y5, "L2"):traity1.animal']*m.4$VCV[,'at.level(y5, "L2"):traity2:at.level(y5, "L2"):traity2.animal'])))
round(HPDinterval(m.4$VCV[,'at.level(y5, "L2"):traity1:at.level(y5, "L2"):traity2.animal']/sqrt((m.4$VCV[,'at.level(y5, "L2"):traity1:at.level(y5, "L2"):traity1.animal']*m.4$VCV[,'at.level(y5, "L2"):traity2:at.level(y5, "L2"):traity2.animal']))))


## BRMS
b.4 <- brm(
  mvbind(y1, y2) ~ (1|p|gr(by=y5, animal, cov = A)), 
  data = d.toy, 
  family = gaussian(), 
  data2 = list(A = A.mat),
  cores=4,
  chains=4, iter = 6000, thin = 3
)
summary(b.4)
b.4 %>% plot
b.4 %>% pp_check(resp = "y1",nsamples = 100)
b.4 %>% pp_check(resp = "y2",nsamples = 100)

summary(b.4)[["random"]][]
summary(b.4)[["spec_pars"]][]
summary(b.4)[["rescor_pars"]][]







##------------- MULTIVARIATE NON-GAUSSIAN, MULTIPLE GROUPING FACTORS -------------##


## MCMCglmm
# NOTE - NEED TO USE DIFFERENT PRIORS FOR NON-GAUSSIAN RESPONSES
p5 <- list(G = list(G1 = list(V = diag(4), nu = 3.002)),
           R = list(V = diag(4), nu = 3.002))

m.5 <- MCMCglmm(cbind(y1, y2, y3, y4) ~ trait-1,
                random = ~us(trait):animal,
                rcov = ~us(trait):units,
                pedigree=t.toy,
                family = c("gaussian","gaussian","poisson", "categorical"), 
                data = d.toy, prior=p5, 
                nitt=210000, burnin=10000, thin=200,
                pr=TRUE,verbose = FALSE)



## BRMS
# remember to create obs column 1:n(row) for obs, i.e. residual error
bf_y1 <- bf(y1 ~ 1 + (1|p|gr(by=y5,animal, cov = A)) + (1|q|obs), sigma = 0.01) + gaussian()
bf_y2 <- bf(y2 ~ 1 + (1|p|gr(by=y5,animal, cov = A)) + (1|q|obs), sigma = 0.01) + gaussian()
bf_y3 <- bf(y3 ~ 1 + (1|p|gr(by=y5,animal, cov = A)) + (1|q|obs)) + poisson()
bf_y4 <- bf(y4 ~ 1 + (1|p|gr(by=y5,animal, cov = A)) + (1|q|obs)) + bernoulli()

b.5 <- brm(
  bf_y1 + bf_y2 + bf_y3 + bf_y4 + set_rescor(FALSE),
  data = d.toy, 
  family = gaussian(), 
  data2 = list(A = A.mat),
  cores=4,
  chains=4, iter = 2000, thin = 1
)
b.5



#### END ####



## --- FUTURE WORK --- ##

# If separate traits in a MV analysis show different magnitudes of phylogenetic signal
# a single A matrix (i.e. assuming BM) will not adequately capture the expected within-trait
# covariance for all traits. Is it possible then to apply different (i.e. transformed) A matrices 
# to each element/block of B to produce a custom mixed Gaussian kronecker product? More importantly,
# can we get MCMCglmm/BRMS to apply different phylo VCV matrices to different traits? (currently 
# passed as a single VCV through A.mat, with kronecker presumably calculated internally)

kronecker(B,A)[1:100,1:100] # variance in y1
kronecker(B,A)[1:100,100:200] # covariance y1y2
kronecker(B,A)[100:200,100:200] # variance in y2
