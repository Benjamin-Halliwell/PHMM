library("MCMCglmm");library("brms");library("tidyverse");library("ape");
library("phytools");library("MASS");library("plyr");library("bindata");library('phangorn')


##--------------------------- ~ ------------------------------##

#### DEFINE SIMULATION PARAMETERS ####

# simulate tree
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

# plot gaussians
plot(d.toy$g1, d.toy$g2)
# table binomials
table(d.toy$b1, d.toy$b2)
table(d.toy$b1, d.toy$b2)/250 # prop

### take a look a probabilities on the inverse-link scale
# plogis(beta[1])
# plogis(beta[1] + a1 + e1) 

### END ####

##------------------- BIVARIATE GAUSSIAN ---------------------##

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
                nitt=350000, burnin=100000, thin=300, 
                pr=TRUE,verbose = FALSE) 
summary(m.1)

## DIAGNOSTICS
# view chain
plot(m.1$VCV)
# calculate autocorrelation among samples
autocorr(m.1$VCV)


### BRMS

b.1 <- brm(
  mvbind(g1, g2) ~ (1|p|gr(animal, cov = A)), 
  data = d.toy,
  data2 = list(A = A.mat),
  family = gaussian(),
  cores = 2,
  chains = 2, iter = 4000, thin = 2
)
summary(b.1)
pairs(b.1)
b.1 %>% plot
b.1 %>% pp_check(resp = "g1",ndraws = 100)
b.1 %>% pp_check(resp = "g2",ndraws = 100)
conditional_effects(b.1)

# saveRDS(b.1, "PMM_brms_2_VAR.rds")



## COMPARISON BETWEEN METHODS
# N.B. Intercept estimates are directly comparable between MCMCglmm and BRMS but the variance components reported
# by MCMCglmm need to be sqrt() to compare to the SDs reported by BRMS (as well as the generating parameters for our sims).

# summary outputs
summary(m.1) # MCMCglmm
summary(b.1) # brms

# intercept
summary(m.1)$solutions
summary(b.1)[["fixed"]]

# phylogenetic variance (sig.B)
sqrt(abs(summary(m.1a)$Gcovariances)) # N.B. remember sign of covariance
summary(b.1)[["random"]]

# residual variance (sig.C)
sqrt(summary(m.1)$Rcovariances)
summary(b.1)[["spec_pars"]]



### END ####

##------------- BIVARIATE (GAUSSIAN, BERNOULLI) --------------##

### FIT MODELS ####

### MCMCglmm

p <- list(B = list(mu=c(0,0), V=diag(c(1,1+pi^2/3))),
          R = list(R1=list(V = diag(2), nu=2, fix=2)),
          G = list(G1=list(V = diag(2), nu=c(0,1000), alpha.mu = c(0,0), alpha.V = diag(c(1,1)))))
m.1 <- MCMCglmm(cbind(g1, b2) ~ trait-1, 
                random = ~us(trait):animal, 
                rcov = ~us(trait):units, 
                pedigree=t.toy,
                family = c("gaussian","categorical"), 
                nodes="ALL", data = d.toy, prior=p,
                nitt=350000, burnin=100000, thin=300, 
                pr=TRUE,verbose = FALSE) 
summary(m.1)

## DIAGNOSTICS
# view chain
plot(m.1$VCV)
# calculate autocorrelation among samples
autocorr(m.1$VCV)


### BRMS

bf_g1 <- bf(g1 ~ 1 + (1|a|gr(animal, cov = A)) + (1|b|obs), sigma = 0.001) + gaussian() # fix residual error at very low level and model through obs instead
bf_b2 <- bf(b2 ~ 1 + (1|a|gr(animal, cov = A)) + (1|b|obs)) + bernoulli()

b.1 <- brm(
  bf_g1 + bf_b2 + set_rescor(FALSE), 
  data = d.toy,
  data2 = list(A = A.mat),
  family = gaussian(),
  cores = 4,
  chains = 4, iter = 6000, thin = 3
)
summary(b.1)
pairs(b.1)
b.1 %>% plot
b.1 %>% pp_check(resp = "g1",ndraws = 100)
b.1 %>% pp_check(resp = "g2",ndraws = 100)
conditional_effects(b.1)

# saveRDS(b.1, "PMM_brms_2_VAR.rds")


## COMPARISON BETWEEN METHODS
# N.B. Intercept estimates are directly comparable between MCMCglmm and BRMS but the variance components reported
# by MCMCglmm need to be sqrt() to compare to the SDs reported by BRMS (as well as the generating parameters for our sims).

# summary outputs
summary(m.1) # MCMCglmm
summary(b.1) # brms

# intercept
summary(m.1)$solutions
summary(b.1)[["fixed"]]

# phylogenetic variance (sig.B)
sqrt(abs(summary(m.1a)$Gcovariances)) # N.B. remember sign of covariance
summary(b.1)[["random"]]

# residual variance (sig.C)
sqrt(summary(m.1)$Rcovariances)
summary(b.1)[["spec_pars"]]


### END ####

##------------- BIVARIATE (BERNOULLI, BERNOULLI) -------------##

### FIT MODELS ####

### MCMCglmm 

# For binary responses, the residual variance is not identifiable. Therefore, the diagonal elements of rcov need to be fixed
# in the prior specification (ref: https://stat.ethz.ch/pipermail/r-sig-mixed-models/2015q3/023966.html). Ideally, we would employ ~corg()
# to fix residual variances while still allowing estimation of residual covariances. Unfortunatley, ~corg() does not seem to behave well for family = categorical.
# Therefore, our options are 1) fix rcov to an identity matrix with family = categorical or 2) use ~corg() with family = threshold


## MV PRIORS
p1a <- list(G = list(G1=list(V=diag(2), nu=1000, alpha.mu=c(0,0), alpha.V=diag(c(1,1)))), # Chi^2 - parameter expanded G priors increase ESS of intercept and phylogenetic variance
          R = list(V=diag(2), nu=0, fix=1))

# p1a <- list(G = list(G1=list(V=diag(2), nu=1000, alpha.mu=c(0,0), alpha.V=diag(c(1,1)))), # Chi^2 prior - parameter expanded G priors increase ESS of intercept and phylogenetic variance
#             R = list(V=diag(c(0.001,0.001)), nu=0, fix=1)) # fix residual variance at some very low level?

p1b <- list(G = list(G1=list(V=diag(2), nu=2)),
          R = list(V=diag(2), nu=0, fix=1))

p1c <- list(G = list(G1=list(V=diag(2), nu=3, alpha.mu=c(0,0), alpha.V=diag(c(1000,1000)))), # scaled Fisher - parameter expanded G priors increase ESS of intercept and phylogenetic variance
          R = list(V=diag(2), nu=0, fix=1))


m.1a<-MCMCglmm(cbind(b1, b2) ~ trait-1,
              random = ~us(trait):animal, 
              rcov = ~us(trait):units,
              family = c("categorical","categorical"),
              pedigree = t.toy,
              data = d.toy, 
              prior = p1a,
              nitt = 350000, burnin = 100000, thin = 300,
              pr = T, verbose = F)
summary(m.1a)


m.1b<-MCMCglmm(cbind(b1, b2) ~ trait-1,
               random = ~us(trait):animal, 
               rcov = ~us(trait):units,
               family = c("categorical","categorical"),
               pedigree = t.toy,
               data = d.toy, 
               prior = p1b,
               nitt = 350000, burnin = 100000, thin = 300,
               pr = T, verbose = F)
summary(m.1b)

m.1c<-MCMCglmm(cbind(b1, b2) ~ trait-1,
               random = ~us(trait):animal, 
               rcov = ~us(trait):units,
               family = c("categorical","categorical"),
               pedigree = t.toy,
               data = d.toy, 
               prior = p1c,
               nitt = 350000, burnin = 100000, thin = 300,
               pr = T, verbose = F)
summary(m.1c)



## THRESHOLD MODEL 
# For bivariate problems, Hadfield recommends use of the 'threshold' family with probit link. Refs:
# https://stat.ethz.ch/pipermail/r-sig-mixed-models/2015q3/023966.html
# https://stat.ethz.ch/pipermail/r-sig-mixed-models/2014q1/021875.html

# PRIORS
p2a <- list(G = list(G1=list(V=diag(2), nu=1000, alpha.mu=c(0,0), alpha.V=diag(c(1,1)))), # parameter expanded G priors increase ESS of intercept and phylogenetic variance
             R = list(V=diag(2), nu=0))

m.2a<-MCMCglmm(cbind(b1, b2) ~ trait-1,
              random = ~us(trait):animal, 
              rcov = ~corg(trait):units, # corg() fixes the residual variances to 1 but estimates the residual correlation. Only possible for threshold?
              family = c("threshold","threshold"),
              pedigree = t.toy,
              data = d.toy, 
              prior = p2a,
              nitt = 35000, burnin = 10000, thin = 30,
              pr = T, verbose = F)
summary(m.2a)


## DIAGNOSTICS
# view chain
plot(m.1$VCV)
# calculate autocorrelation among samples. All values above lag 0 should be close to 0
autocorr(m.1$VCV)



## RESCALE BASED ON RESIDUAL VARIANCE ? ##
# For binomial, family variance only makes sense in context of the residual variance.
# Use intra-class correlation coefficient to rescale the phylogenetic variance 
# according to the assumed residual variance (fixed at 1)
head(m.1$VCV)

ICC.1 <- m.1$VCV[, 1]/(m.1$VCV[, 1] + 1 + pi^2/3) # 1 = fixed additive overdispersion
mean(ICC.1);HPDinterval(ICC.1)

ICC.2 <- m.1$VCV[, 4]/(m.1$VCV[, 4] + 1 + pi^2/3)
mean(ICC.2);HPDinterval(ICC.2)


# rescale estimates by the estimated residual variance in order to obtain the
# posterior distributions of the parameters under the assumption that the 
# actual residual variance is equal to some other value.
c2 <- ((16 * sqrt(3))/(15 * pi))^2 # constant for logit link, c = 1 for probit
Sol.s <- m.1$Sol/sqrt(1 + c2 * m.1$VCV[, 5])
head(colnames(m.1$Sol))
mean(Sol.s[,"traitb1.1"]);HPDinterval(Sol.s[,"traitb1.1"])
mean(Sol.s[,"traitb2.1"]);HPDinterval(Sol.s[,"traitb2.1"])


#---------------------------------------------------------------------#

### BRMS


## PRIORS

# see which priors brms is using by default for our desired model
priors_b.1<-get_prior(mvbind(b1,b2) ~ 1 + (1|a|gr(animal, cov = A)) + (1|obs), # (1|obs) fits independent obs level group effects for each response
                      data = d.toy,
                      data2 = list(A = A.mat),
                      family = bernoulli())  
priors_b.1


# + (1|b|obs)) manually added allows for residual (co)variances to be estimated on the link scale, as
# residual variance (additive over dispersion) is not fit by default for a bernoulli response in brms

b.1 <- brm(
  bf(mvbind(b1,b2) ~ 1 + (1|a|gr(animal, cov = A)) + (1|obs)) + bernoulli(), # + bernoulli() redundant if family is set to bernoulli() (must specify separate bf() if prob dist of responses differs)
  data = d.toy,
  prior = set_prior("constant(1)", class = "sd", group = "obs", resp = c("b1","b2")),
  family = bernoulli(),
  control = list(adapt_delta = 0.9),
  data2 = list(A = A.mat),
  cores=4,
  chains=4, iter = 100, thin = 1
)
b.1


b.2 <- brm(
  bf(mvbind(b1,b2) ~ 1 + (1|a|gr(animal, cov = A))),
  data = d.toy, 
  family = bernoulli(),
  control = list(adapt_delta = 0.9),
  data2 = list(A = A.mat),
  cores=4,
  chains=4, iter = 1000, thin = 1
)
b.2

# b.3 %>% saveRDS("b.3.rds")
b.1 %>% pairs

## DIAGNOSTICS
b.1 %>% plot
b.1 %>% pp_check(nsamples = 100)
pairs(b.1)


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
sqrt(abs(summary(m.1a)$Gcovariances)) # N.B. remember sign of covariance
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

##--------------------------- ~ ------------------------------##

## NOTES ####

## UV PRIORS
# p <- list(G = list(G1=list(V=1, nu=0.002)), # Inv Wish for random intercepts (MCMCglmm default)
#           R = list(V = 1, fix = 1))
# p <- list(G = list(G1=list(V=1, nu=1000, alpha.mu=0, alpha.V=1)), # parameter expanded G priors increase ESS of intercept and phylogenetic variance
#           R = list(V = 1, fix = 1))
# p <- list(G = list(G1=list(V=1, nu=1, alpha.mu=0, alpha.V=1)), # reducing nu reduces ESS and gives inflated estimates of phylogenetic variance
#           R = list(V = 1, fix = 1))
# p <-  list(B = list(mu = c(0, 0), V = diag(2) * (1 + pi^2/3)), # prior to use when incomplete separation is a problem for fixed effects
#                    R = list(V = 1, fix = 1))
# p <- list(G = list(G1=list(V=1, nu=1, alpha.mu=0, alpha.V=1000)), # parameter expanded prior recommended for binary data in MCMCMglmm course notes reduces ESS and gives inflated estimates of phylogenetic variance
#           R = list(V = 1,  fix = 1))


# p <- list(R = list(R1=list(V = diag(1), nu=0.002)),
#           G = list(G1=list(V = diag(1), nu=0.002)))
# m.1 <- MCMCglmm(g1 ~ g2, 
#                 random = ~animal, 
#                 rcov = ~units, 
#                 pedigree=t.toy,
#                 family = c("gaussian"), 
#                 nodes="ALL", data = d.toy, prior=p,
#                 nitt=35000, burnin=10000, thin=30, 
#                 pr=T, saveX = T, saveZ = T,
#                 verbose = F) 
# summary(m.1)
# 
# head(d.toy)
# # conditional mean
# W.1<-cBind(m.1$X, m.1$Z) # note X and Z are sparse so use cBind
# prediction.1<-W.1%*%posterior.mode(m.1$Sol)
# plot(d.toy$g1,prediction.1)
# plot(order(as.factor(d.toy$animal)), prediction.1@x)
# 
# # fixed effects only
# prediction.2<-m.1$X%*%posterior.mode(m.1$Sol[,1:2])
# plot(d.toy$g1,prediction.2);abline(0,1)
# plot(order(as.factor(d.toy$animal)), prediction.2@x)
# plot(order(as.factor(d.toy$animal)), d.toy$g2)
# 
# # cond mean has higher variance due to consideration of group level effects
# plot(prediction.2,prediction.1)
# var(prediction.1@x);var(prediction.2@x)
# 
# prediction.2@x
# 
# # group level effects
# plot(order(as.factor(d.toy$animal)),prediction.1-prediction.2);abline(0,0)
# 
# 
# ## using predict()
# #
# dummy <- data.frame(g2 = seq(0,1, length.out = 10),
#                     animal = paste0("s",rep(1:250, each = 10)),
#                     y=0)
# 
# # conditional on group level effects
# g1.c <- predict(m.1, dummy, marginal = NULL)
# # marginalizing over group level effects
# g1.m <- predict(m.1, dummy)
# 
# df.pred <- data.frame(animal=dummy$animal, g2=dummy$g2, g1.c, g1.m)
# 
# lattice::xyplot(g1.c + g1.m ~ g2 | animal, data = df.pred)[1:10]
# 
# head(df.pred)
# head(d.toy)
summary(m.1)

## END ####

#--------------------------------------------------------------#

#### DATA HANDLING ####

# READ IN DATA AND TREE
d <- read.csv('Liz_Master_FINAL.csv')
d <- d[d$taxon=='Sauria'| d$taxon=='Sphenodontia',] # subset to lizards. include Sphenodon as outgroup.
d <- d[complete.cases(d$SG),] # subset to complete cases for SG
tree <- read.tree('Liz_comp_FINAL.tre')
is.ultrametric(tree)
tree.u <- force.ultrametric(tree)# make tree ultrametric

# DATA PREP
d$animal <- d$species # set up reserved term 'animal'
d$parity <- as.factor(d$parity)
d$parity <- revalue(d$parity, c('0' = 'Ov', '1' = 'Viv'))
d$Lat <- abs(d$Lat)
d$region <- as.factor(ifelse(d$Lat<23.5, 'Tropics', (as.character(d$Hemisphere))))
d$region2 <- as.factor(ifelse(d$region=='Tropics', 'Tropics', 'Temperate'))

# TRANSFORMATIONS (do we also need to scale?)
d$SVL <- log(d$SVL)
d$mass <- log(d$mass)
d$juvSVL <- log(d$juvSVL)
d$juv_mass <- log(d$juv_mass)
d$GmeanCS <- round(d$GmeanCS) # grand mean of clutch size observations rounded to nearest integer
d$meanMature <- log(d$meanMature)

# NAs TREATED AS ABSENCE FOR HABITAT ASSOCIATION (i.e., NOT FOUND IN THAT HABITAT TYPE)
d$arboreal[is.na(d$arboreal)] <- 0
d$saxicolous[is.na(d$saxicolous)] <- 0
d$terrestrial[is.na(d$terrestrial)] <- 0
d$fossorial[is.na(d$fossorial)] <- 0


#### END ####

#### LIZ MODEL CHECK ####

## AUTO CORRELATION CORRECTED SE FOR CLIM VARIABLES?

# SAVED RUNS
load(file="m.1");load(file="m.2")

summary(d)
d <- d[complete.cases(d$pMeanPred),] # some missing data in pMeanPred

p1 <- list(B = list(mu=rep(0,10), V = diag(c(1+pi^2/3,rep(1,9)))),
           R = list(R1=list(V = diag(1), fix = 1)),
           G = list(G1=list(V = diag(1), nu=1000, alpha.mu=0,alpha.V=1)))

m.1 <- MCMCglmm(SG ~ 1 + pMean + tMean + pMeanPred + tMeanPred + 
                  fossorial + saxicolous + terrestrial + arboreal +
                  parity,
                random = ~animal,
                rcov = ~units,
                pedigree=tree.u,
                family = c("categorical"),
                nodes="ALL", data=d, prior=p1,
                nitt=600000, burnin=100000, thin=500,
                pr=T, saveX = T, saveZ = T,
                verbose = F)
summary(m.1)
# save(m.1, file="m.1")

d.b <- d[complete.cases(d$saxicolous),]
p1b <- list(B = list(mu=rep(0,9), V = diag(c(1+pi^2/3,rep(1,8)))),
           R = list(R1=list(V = diag(1), fix = 1)),
           G = list(G1=list(V = diag(1), nu=1000, alpha.mu=0,alpha.V=1)))

m.1b <- MCMCglmm(SG ~ 1 + pMean + tMean + pMeanPred + tMeanPred + 
                  fossorial + saxicolous + terrestrial + arboreal,
                random = ~animal,
                rcov = ~units,
                pedigree=tree.u,
                family = c("categorical"),
                nodes="ALL", data=d.b, prior=p1b,
                nitt=60000, burnin=10000, thin=50,
                pr=T, saveX = T, saveZ = T,
                verbose = F)
summary(m.1b)
# save(m.1b, file="m.1b")

# save(m.1b, file="m.1b")

#
p2 <- list(B=list(mu=rep(0,4), V=diag(c(1,1+pi^2/3,1+pi^2/3,1+pi^2/3))),
        R = list(R1=list(V = diag(4), nu=4.002, fix = 2)),
        G = list(G1=list(V = diag(4), nu=1000, alpha.mu = c(0,0,0,0), alpha.V = diag(c(1,1,1,1)))))

p2.2 <- list(B=list(mu=rep(0,4), V=diag(c(1,1+pi^2/3,1+pi^2/3,1+pi^2/3))),
           R = list(R1=list(V = diag(4), nu=4.002, fix = 2)),
           G = list(G1=list(V = diag(4), nu=3, alpha.mu = c(0,0,0,0), alpha.V = diag(c(1,1000,1000,1000)))))


m.2 <- MCMCglmm(cbind(tMean, SG, saxicolous, parity) ~ trait-1,
                random = ~us(trait):animal,
                rcov = ~idh(trait):units,
                pedigree=tree.u,
                family = c("gaussian","categorical","categorical","categorical"),
                nodes="ALL", data=d, prior=p2,
                nitt=600000, burnin=100000, thin=500,
                pr=T, saveX = T, saveZ = T,
                verbose = F)
summary(m.2)


m.2.2 <- MCMCglmm(cbind(tMean, SG, saxicolous, parity) ~ trait-1,
                random = ~us(trait):animal,
                rcov = ~idh(trait):units,
                pedigree=tree.u,
                family = c("gaussian","categorical","categorical","categorical"),
                nodes="ALL", data=d, prior=p2.2,
                nitt=600000, burnin=100000, thin=500,
                pr=T, saveX = T, saveZ = T,
                verbose = F)

result <- m.2$VCV[,"traitSG.1:traitparity.2.animal"]/
  sqrt(m.2$VCV[,"traitSG.1:traitSG.1.animal"]*
         m.2$VCV[,"traitparity.2:traitparity.2.animal"])
posterior.mode(result);HPDinterval(result)

#
p3 <- list(B=list(mu=rep(0,3), V=diag(c(1,1+pi^2/3,1+pi^2/3,1,1+pi^2/3,1+pi^2/3))),
        R = list(R1=list(V = diag(3), nu=3.002, fix = 2),
                R2=list(V = diag(3), nu=3.002, fix = 2)),
        G = list(G1=list(V = diag(3), nu=1000, alpha.mu = c(0,0,0), alpha.V = diag(c(1,1,1))),
                G2=list(V = diag(3), nu=1000, alpha.mu = c(0,0,0), alpha.V = diag(c(1,1,1)))))

m.3 <- MCMCglmm(cbind(tMean, SG, saxicolous) ~ trait:parity-1
                random = ~us(at.level(parity,"Ov"):trait):animal + us(at.level(parity,"Viv")trait):animal,
                rcov = ~idh(at.level(parity,'Ov'):trait):units + idh(at.level(parity,'Viv'):trait):units,
                pedigree=tree.u,
                family = c("gaussian","categorical","categorical"),
                nodes="ALL", data=d, prior=p3,
                nitt=600000, burnin=100000, thin=500,
                pr=T, saveX = T, saveZ = T,
                verbose = F)






###########################################################.
### Calculate partial R2 (Nakagawa and Schielzeth 2013) ###
###########################################################.

### NON-GAUSSIAN RESPONSE (BINARY) ###

# ## INTERCEPT ONLY MODEL
# 
# fit <- m.0
# 
# # calculate posterior mean of total phenotypic variance (adding additive over dispersion (1) and distribution-specific variance (pi^2/3)) and 95% CRI
# V = as.mcmc(fit$VCV[,1] + 1 + pi^2/3)
# mean(V);HPDinterval(V)
# 
# # calculate proportion of variance (point estimates and 95% CRI) attributed to each variance component
# posterior.mode(fit$VCV/V); HPDinterval(fit$VCV/V)
# head(fit$VCV)


## MODEL INCLUDING FIXED EFFECTS

summary(m.1)
summary(m.2)


fit <- m.1

# summary(fit)
vmVarF<- numeric(length(fit$Sol[,1])) # number = length(fit$Sol[,1])
for(i in 1:length(fit$Sol[,1])){
  Var<-var(as.vector(fit$Sol[i,1:nrow(summary(fit)$solutions)] %*% t(fit$X))) # Sol[i, 1:k]
  vmVarF[i]<-Var}

# calculate marginal and conditional R^2 adding distribution specific variance (for binary, pi^2/3)
# head(fit$VCV)
R2m <- vmVarF/(vmVarF+fit$VCV[,1]+fit$VCV[,2]+pi^2/3) # variance explained by fixed effects divided by total variance
R2c <- (vmVarF+fit$VCV[,1])/(vmVarF+fit$VCV[,1]+fit$VCV[,2]+pi^2/3) # variance explained by fixed and random effects combined

## RESULTS
plot(R2c);posterior.mode(R2c); HPDinterval(R2c) # prop of total phenotypic variance explained by model
plot(R2m);posterior.mode(R2m); HPDinterval(R2m) # prop explained by fixed effects
posterior.mode(R2c-R2m) # prop explained by group level (random) effects


#################.
## PREDICTION ###
#################.

mod <- m.1
dat <- d

# full conditional prediction
W.1<-cbind(mod$X, mod$Z) # note X and Z are sparse so use cbind
pred.c<-W.1%*%posterior.mode(mod$Sol) # multiply design matrix by posterior coefficient estimates
# prediction from fixed effects only
pred.f<-mod$X%*%posterior.mode(mod$Sol[,1:nrow(summary(mod)$solutions)])
# put predictions on the probability scale
dat$pSG.c <- plogis(pred.c@x)
dat$pSG.f <- plogis(pred.f@x)

# plot(dat$SG,pred.c)
# plot(order(as.factor(dat$animal)), pred.c@x) # viviparous clades seen as vertical striations
# plot(dat$SG,pred.f)
# plot(order(as.factor(dat$animal)), pred.f@x)

# compare the association between SG and predictors between 
# conditional and fixed effect predictions

# predicted probabilities of SG on observed states
ggplot(dat, aes(x = as.factor(SG), y = pSG.f, colour = parity)) +
  geom_boxplot(position = position_dodge(width = 0.5))  +
  theme_classic() + scale_color_manual(values = c("blue", "red"))
ggplot(dat, aes(x = as.factor(SG), y = pSG.c, colour = parity)) +
  geom_boxplot(position = position_dodge(width = 0.5)) +
  theme_classic() + scale_color_manual(values = c("blue", "red"))

# on tMean
ggplot(dat, aes(x = tMean, y = pSG.f, colour = parity:as.factor(SG))) +
  geom_point( size = 2) + theme_classic() + scale_color_manual(values = c("light blue","blue", "pink", "red"))
ggplot(dat, aes(x = tMean, y = pSG.c, colour = parity:as.factor(SG))) +
  geom_point(size = 2) + theme_classic() + scale_color_manual(values = c("light blue","blue", "pink", "red"))

# on pMean
ggplot(dat, aes(x = pMean, y = pSG.f, colour = parity:as.factor(SG))) +
  geom_point( size = 2) + theme_classic() + scale_color_manual(values = c("light blue","blue", "pink", "red"))
ggplot(dat, aes(x = pMean, y = pSG.c, colour = parity:as.factor(SG))) +
  geom_point(size = 2) + theme_classic() + scale_color_manual(values = c("light blue","blue", "pink", "red"))
+ ylim(0,1)





# cond mean has higher variance due to consideration of group level effects
plot(dat$pSG.c,dat$pSG.f)
var(dat$pSG.c);var(dat$pSG.f)


# PICKING "COLD" THRESHOLD FOR BINARY T VARIABLE
ggplot(d, aes(x = tMean, y = pMean, colour = parity:as.factor(SG))) +
  geom_point(size = 2) + theme_classic() + scale_color_manual(values = c("light blue","blue", "pink", "red"))


### END ####
