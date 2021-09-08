library("MCMCglmm");library("brms");library("tidyverse");library("ape");
library("phytools");library("MASS");library("plyr");library("bindata")


load(file="t.toy")
load(file="d.toy")

### SIMULATE DATA AND TREE ####

## simulate tree
t.toy <- geiger::sim.bdtree(b=1, d=0, stop="taxa", n=100, extinct=FALSE) # may be best not to use pure birth trees (Adams and Collyer 2018)
t.toy <- multi2di(t.toy)
# save(t.toy, file="t.toy")

# Create VCV matrix from tree
A.mat <- ape::vcv.phylo(t.toy, corr=T) # BM

# define A.mat under the assumption of OU

## define simulation parameters
# intercepts for each trait
beta = c(1,2,1,-2,-1)
# number of traits
k = 5

# B/b is phylogenetic trait-level VCV matrix
# ? Do the elements of sig.B represent the strength of signal in each trait ?
# ? how would we define interactions between variables, e.g. b14_rho is higher when y5 = 1?
sig.B <- c(b11 = 1, b22 = 1, b33 = 1, b44 = 2, b55 = 2) # standard deviation of diagonal components
b12_rho = 0.75; b13_rho = 0.5; b14_rho = 0.5; b15_rho = 0.5 # off-diagonal correlation coefficients
b23_rho = 0.5; b24_rho = 0.5; b25_rho = 0.5; 
b34_rho = 0.5; b35_rho = 0.5; b45_rho = 1

Bcor <- matrix(c(c(1,b12_rho,b13_rho,b14_rho,b15_rho), # correlation matrix
                 c(b12_rho,1,b23_rho,b24_rho,b25_rho),
                 c(b13_rho,b23_rho,1,b34_rho,b35_rho),
                 c(b14_rho,b24_rho,b34_rho,1,b45_rho),
                 c(b15_rho,b25_rho,b35_rho,b45_rho,1)),k,k, byrow = T) 

B <- matrix(kronecker(sig.B, sig.B),k,k)*Bcor # VCV


# C/c is the residual trait-level VCV matrix
sig.C <- c(c11 = 1, c22 = 1, c33 = 1, c44 = 1, c55 = 1) # standard deviation of diagonal components 
c12_rho = 0.5; c13_rho = 0.5; c14_rho = 0.5; c15_rho = 0.5 # off-diagonal correlation coefficients
c23_rho = 0.5; c24_rho = 0.5; c25_rho = 0.5; 
c34_rho = 0.5; c35_rho = 0.5; c45_rho = 0.5

Ccor <- matrix(c(c(1,c12_rho,c13_rho,c14_rho,c15_rho), # correlation matrix
                 c(c12_rho,1,c23_rho,c24_rho,c25_rho),
                 c(c13_rho,c23_rho,1,c34_rho,c35_rho),
                 c(c14_rho,c24_rho,c34_rho,1,c45_rho),
                 c(c15_rho,c25_rho,c35_rho,c45_rho,1)),k,k, byrow = T) 

C <- matrix(kronecker(sig.C, sig.C),k,k)*Ccor # VCV

# supplied phylogenetic VCV matrix
A = A.mat

# number of species
n = nrow(A)

# identitiy matrix
I = diag(n)

# simulate phylogenetic random effects
a <- mvrnorm(1,rep(0,n*k),kronecker(B,A)) # animal effect enters via A. 
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

# define response (apply inverse link)
y1 <- beta[1] + a1 + e1
y2 <- beta[2] + a2 + e2
y3 <- rpois(n,exp(beta[3] + a3 + e3))
y4 <- rbinom(n,1,plogis(beta[4] + a4 + e4)) # plogis(beta[4] + a4 + e4)?
y5 <- rbinom(n,1,plogis(beta[5] + a5 + e5))

# generate df
species <- t.toy$tip.label
d.toy <- data.frame(species,y1,y2,y3,y4,y5) # as second binary (parity)
d.toy$animal <- d.toy$species
d.toy$y5 <- as.factor(y5) # second Bernoulli as factor for separate VCVs via at.level/by coding
d.toy$y5 <- revalue(d.toy$y5, c("0"="L1", "1"="L2"))
d.toy$obs <- 1:nrow(d.toy)
# save(d.toy, file="d.toy")
head(d.toy)
summary(d.toy)

table(d.toy$y4, d.toy$y5)
table(d.toy$y4)
table(d.toy$y5)

### END ####


### PLOTTING ####

# plot data
plot(y1~y2, d.toy)
plot(y1~y3, d.toy)
# phylo dist of y1
names(y1) <- t.toy$tip.label
plotTree.barplot(t.toy, y1)
# phylo dist of y3
names(y3) <- t.toy$tip.label
par(mfrow=c(1,1), mar=c(2,2,2,2))
dotTree(t.toy, y3)


# trait level cor is hidden by phylogenetic cov
par(mfrow=c(1,1), mar=c(2,2,2,2))
plot(y1,y2)
# view trait level cor for fixed animal
traitT <- mvrnorm(100,rep(0,n*nrow(B)),kronecker(B^2,A))
plot(traitT[,n+1],traitT[,1]) # 100 samples of both t1 and t2 for a given animal reveals correlation

### END ####


#### PHYLOGENETIC SIGNAL ####

phylosig(t.toy, x = d.toy$y1, method= "lambda", test=T, nsim=10000)
phylosig(t.toy, x = d.toy$y1, method= "K", test=T, nsim=10000)

phylosig(t.toy, x = d.toy$y2, method= "lambda", test=T, nsim=10000)
phylosig(t.toy, x = d.toy$y2, method= "K", test=T, nsim=10000)


#### END ####


### MODELS ###

##------------- UNIVARIATE -------------##

## MCMCglmm 
p <- list(G = list(G1 = list(V = 1, nu = 0.002)), 
          R = list(V = 1, nu = 0.002))

m.1<-MCMCglmm(y1 ~ 1,
                random = ~animal,
                pedigree=t.toy, 
                family = c("gaussian"), 
                data = d.toy, prior = p,
                nitt=210000, burnin=10000, thin=200,
                pr=TRUE,verbose = FALSE)
summary(m.1)
par(mar=(c(2,2,2,2)))
plot(m.1$Sol)
plot(m.1$VCV)
# calculate autocorrelation among samples
# - seeking all values above lag 0 to be close to 0
autocorr(m.1$VCV)
# obtain estimates of the additive genetic/phylogenetic and residual variances
posterior.mode(m.1$VCV)
# credible intervals
HPDinterval(m.1$VCV)

  
## BRMS
mb.1 <- brm(
  y1 ~ 1 + (1|gr(animal, cov = A)), 
  data = d.toy, 
  family = gaussian(), 
  data2 = list(A = A.mat),
  cores=4,
  chains=4, iter = 6000, thin = 3
)
mb.1
mb.1 %>% plot
mb.1 %>% pp_check(nsamples = 100)
pairs(mb.1)
names(as.data.frame(mb.1))
# compare estimates between methods
summary(m.1)
summary(mb.1)

# intercept estimate is directly comparable but variance components 
# need to be sqrt to compare MCMCglmm outputs to BRMS (and generating parameters)

# phylo
paste0(round(sqrt(mean(m.1$VCV[,"animal"])),2)," (",round(sqrt(HPDinterval(m.1$VCV[,"animal"]))[1],2),
       ", ",round(sqrt(HPDinterval(m.1$VCV[,"animal"]))[2],2),")")
# residual
paste0(round(sqrt(mean(m.1$VCV[,"units"])),2)," (",round(sqrt(HPDinterval(m.1$VCV[,"units"]))[1],2),
       ", ",round(sqrt(HPDinterval(m.1$VCV[,"units"]))[2],2),")")

summary(mb.1)[["random"]][]
summary(mb.1)[["spec_pars"]][]

## NOTES
#
# Both methods provide good point estimates
# and similar CIs
#


##--------- MULTIVARIATE -----------##

## MCMCglmm
#
# nu > 1 for proper prior
# result sensitive to value of V
#

p2 <- list(G = list(G1 = list(V = diag(1,2), nu = 1.002)), 
          R = list(V = diag(1,2), nu = 1.002))

m.m.1<-MCMCglmm(cbind(y1, y2) ~ trait-1,
                random = ~us(trait):animal,
                rcov = ~us(trait):units,
                pedigree=t.toy,
                family = c("gaussian","gaussian"), 
                data = d.toy, prior=p2, 
                nitt=210000, burnin=10000, thin=200,
                pr=TRUE,verbose = FALSE)
summary(m.m.1)
plot(m.m.1$VCV)
autocorr(m.m.1$VCV)[,,1]

# PHYLOGENETIC HERITABILITY #
posterior.mode(m.m.1$VCV[,'traity1:traity1.animal']/(m.m.1$VCV[,'traity1:traity1.animal']+m.m.1$VCV[,'traity1:traity1.units']))
posterior.mode(m.m.1$VCV[,'traity2:traity2.animal']/(m.m.1$VCV[,'traity2:traity2.animal']+m.m.1$VCV[,'traity2:traity2.units']))
# PHYLOGENETIC CORRELATION #
posterior.mode(m.m.1$VCV[,'traity1:traity2.animal']/sqrt((m.m.1$VCV[,'traity1:traity1.animal']*m.m.1$VCV[,'traity2:traity2.animal'])))


## BRMS
#
#
#

m.mb.1 <- brm(
  mvbind(y1, y2) ~ (1|p|gr(animal, cov = A)), 
  data = d.toy, 
  family = gaussian(), 
  data2 = list(A = A.mat),
  cores=4,
  chains=4, iter = 6000, thin = 3
)
summary(m.mb.1)
m.mb.1 %>% plot
m.mb.1 %>% pp_check(resp = "y1",nsamples = 100)
m.mb.1 %>% pp_check(resp = "y2",nsamples = 100)


## COMPARE MODEL OUTPUTS
mean(m.m.1$Sol[,"traity1"]);mean(m.m.1$Sol[,"traity2"])
summary(m.mb.1)[["fixed"]][]
# phylo
sqrt(mean(m.m.1$VCV[,"traity1:traity1.animal"]))
sqrt(mean(m.m.1$VCV[,"traity2:traity2.animal"]))
sqrt(mean(m.m.1$VCV[,"traity1:traity2.animal"])) # underestimates - reports cov. Should this be sqrt()?
summary(m.mb.1)[["random"]][] # underestimates - reports cor
# res
sqrt(mean(m.m.1$VCV[,"traity1:traity1.units"]))
sqrt(mean(m.m.1$VCV[,"traity2:traity2.units"]))
summary(m.mb.1)[["spec_pars"]][]
# res covar
sqrt(mean(m.m.1$VCV[,"traity1:traity2.units"])) # overestimates
summary(m.mb.1)[["rescor_pars"]][] # overestimates


## NOTES
#
# fairly close agreement between models
# - BRMS underestimates phylo cov
# - MCMCglmm overestimates residual cov
#


##---------- MULTIVARIATE, MULTIPLE GROUPING FACTORS -----------##

# MCMCglmm
#
p3=list(G = list(G1=list(V = diag(2), nu = 1.002),
                 G2=list(V = diag(2), nu = 1.002)),
        R = list(R1=list(V = diag(2), nu=1.002),
                 R2=list(V = diag(2), nu=1.002)))

m.m.2<-MCMCglmm(cbind(y1, y2) ~ trait:y5-1,
                random = ~us(at.level(y5,'L1'):trait):animal + us(at.level(y5,'L2'):trait):animal,
                rcov = ~us(at.level(y5,'L1'):trait):units + us(at.level(y5,'L2'):trait):units,
                family = c("gaussian","gaussian"), 
                pedigree=t.toy, data = d.toy, prior=p3, 
                nitt=210000, burnin=10000, thin=200,
                pr=TRUE,verbose = FALSE)
summary(m.m.2)
# plot(m.m.2$VCV)
autocorr(m.m.2$VCV)[,1:3,1]

mean(m.m.2$VCV[,'at.level(y5, "L1"):traity1:at.level(y5, "L1"):traity2.animal']) # negative 

# PHYLOGENETIC CORRELATION #
# L1
posterior.mode(m.m.2$VCV[,'at.level(y5, "L1"):traity1:at.level(y5, "L1"):traity2.animal']/sqrt((m.m.2$VCV[,'at.level(y5, "L1"):traity1:at.level(y5, "L1"):traity1.animal']*m.m.2$VCV[,'at.level(y5, "L1"):traity2:at.level(y5, "L1"):traity2.animal'])))
round(HPDinterval(m.m.2$VCV[,'at.level(y5, "L1"):traity1:at.level(y5, "L1"):traity2.animal']/sqrt((m.m.2$VCV[,'at.level(y5, "L1"):traity1:at.level(y5, "L1"):traity1.animal']*m.m.2$VCV[,'at.level(y5, "L1"):traity2:at.level(y5, "L1"):traity2.animal']))),2)[1:2]
# L2
posterior.mode(m.m.2$VCV[,'at.level(y5, "L2"):traity1:at.level(y5, "L2"):traity2.animal']/sqrt((m.m.2$VCV[,'at.level(y5, "L2"):traity1:at.level(y5, "L2"):traity1.animal']*m.m.2$VCV[,'at.level(y5, "L2"):traity2:at.level(y5, "L2"):traity2.animal'])))
round(HPDinterval(m.m.2$VCV[,'at.level(y5, "L2"):traity1:at.level(y5, "L2"):traity2.animal']/sqrt((m.m.2$VCV[,'at.level(y5, "L2"):traity1:at.level(y5, "L2"):traity1.animal']*m.m.2$VCV[,'at.level(y5, "L2"):traity2:at.level(y5, "L2"):traity2.animal']))))


## BRMS
#
m.mb.2 <- brm(
  mvbind(y1, y2) ~ (1|p|gr(by=y5, animal, cov = A)), 
  data = d.toy, 
  family = gaussian(), 
  data2 = list(A = A.mat),
  cores=4,
  chains=4, iter = 6000, thin = 3
)
summary(m.mb.2)
m.mb.2 %>% plot
m.mb.2 %>% pp_check(resp = "y1",nsamples = 100)
m.mb.2 %>% pp_check(resp = "y2",nsamples = 100)

summary(m.mb.2)[["random"]][]
summary(m.mb.2)[["spec_pars"]][]
summary(m.mb.2)[["rescor_pars"]][]

#--------------------------- MULTIVARIATE WITH MIXED DISTRIBUTIONS, MULTIPLE GROUPING FACTORS ----------------------------------#


# remember to create obs column 1:n(row) for obs, i.e. residual error
bf_y1 <- bf(y1 ~ 1 + (1|p|gr(by=y5,animal, cov = A)) + (1|q|obs), sigma = 0) + gaussian()
bf_y2 <- bf(y2 ~ 1 + (1|p|gr(by=y5,animal, cov = A)) + (1|q|obs), sigma = 0) + gaussian()
bf_y3 <- bf(y3 ~ 1 + (1|p|gr(by=y5,animal, cov = A)) + (1|q|obs)) + poisson()
bf_y4 <- bf(y4 ~ 1 + (1|p|gr(by=y5,animal, cov = A)) + (1|q|obs)) + bernoulli()

m.mb.3 <- brm(
  bf_y1 + bf_y2 + bf_y3 + bf_y4 + set_rescor(FALSE),
  data = d.toy, 
  family = gaussian(), 
  data2 = list(A = A.mat),
  cores=4,
  chains=4, iter = 2000, thin = 1
)

m.mb.3

brms::set_rescor()
bf_y1 <- bf(y1 ~ 1  + (1|q|obs), sigma = 0.01) + gaussian()
bf_y2 <- bf(y2 ~ 1  + (1|q|obs), sigma = 0.01) + gaussian()
bf_y3 <- bf(y3 ~ 1  + (1|q|obs)) + poisson()
bf_y4 <- bf(y4 ~ 1  + (1|q|obs)) + bernoulli()

m.mb.3 <- brm(
  bf_y1 + bf_y2 + bf_y3 + bf_y4 + set_rescor(FALSE),
  data = d.toy, 
  family = gaussian(), 
  data2 = list(A = A.mat),
  cores=4,
  chains=4, iter = 2000, thin = 1
)

m.mb.3


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
