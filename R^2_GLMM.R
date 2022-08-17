######################################################################
### Calculate partial R2 (Nakagawa and Schielzeth 2013) - MCMCglmm ###
######################################################################

library(arm)

## GAUSSIAN RESPONSE 

fit <- UV.tMean.1
summary(fit)

# create vector to hold sigma_f values
vmVarF<- numeric(length(fit$Sol[,1])) # number = length(fit$Sol[,1])

# predict fitted values based on fixed effects alone by multiplying the 
# design matrix of the fixed effects with the vector of fixed effect estimates
# then calculate variance of these fitted values. repeat for each MCMC sample
for(i in 1:length(fit$Sol[,1])){
  Var<-var(as.vector(fit$Sol[i,1:nrow(summary(fit)$solutions)] %*% t(fit$X))) # Sol[i, 1:k]
  vmVarF[i]<-Var}

# stats
summary(vmVarF)
hist(vmVarF) # distribution of sigma_f (variance of fixed effects) across samples
m <- mean(vmVarF)
s <- sd(vmVarF)

# calculate marginal and conditional R^2
head(fit$VCV)
R2m <- vmVarF/(vmVarF+fit$VCV[,1]+fit$VCV[,2]) # variance explained by fixed effects divided by total variance
R2c <- (vmVarF+fit$VCV[,1])/(vmVarF+fit$VCV[,1]+fit$VCV[,2]) # variance explained by fixed and random effects combined

posterior.mode(R2m); HPDinterval(R2m) # variance explained by fixed effects
posterior.mode(R2c); HPDinterval(R2c) # variance explained by the entire model (fixed + random effects)

# calculate PCV (proportion change in variance) between null (intercept only) and full model
# for each level of the model hierarchy
fit0 <- UV.tMean.0 # null model

PCVran <- 1 - (fit$VCV[,1] / fit0$VCV[,1])
PCVres <- 1 - (fit$VCV[,2] / fit0$VCV[,2])

# positive PCV value indicates that variance explained by the focal level has decreased
# with the addition of the specified fixed effects, negative PCV indicates increase in variance explained.
posterior.mode(PCVran); HPDinterval(PCVran)
posterior.mode(PCVres); HPDinterval(PCVres)


## NON-GAUSSIAN RESPONSE (BINARY)

fit <- UV.SG.1
summary(fit)
vmVarF<- numeric(length(fit$Sol[,1])) # number = length(fit$Sol[,1])
for(i in 1:length(fit$Sol[,1])){
  Var<-var(as.vector(fit$Sol[i,1:nrow(summary(fit)$solutions)] %*% t(fit$X))) # Sol[i, 1:k]
  vmVarF[i]<-Var}

# calculate marginal and conditional R^2 adding distribution specific variance (for binary, pi^2/3)
head(fit$VCV)
R2m <- vmVarF/(vmVarF+fit$VCV[,1]+fit$VCV[,2]+pi^2/3) # variance explained by fixed effects divided by total variance
R2c <- (vmVarF+fit$VCV[,1])/(vmVarF+fit$VCV[,1]+fit$VCV[,2]+pi^2/3) # variance explained by fixed and random effects combined

posterior.mode(R2m); HPDinterval(R2m)
posterior.mode(R2c); HPDinterval(R2c)



# calculate PCV adding distribution specific variance
fit0 <- UV.SG.0 # null model

PCVran <- 1 - (fit$VCV[,1] / fit0$VCV[,1])
# PCVres <- 1 - ((fit$VCV[,2]+pi^2/3) / (fit0$VCV[,2]+pi^2/3)) # PVCres not applicable to binary response

# positive PCV value indicates that variance explained by the focal level has decreased
# with the addition of the specified fixed effects, negative PCV indicates increase in variance explained.
posterior.mode(PCVran); HPDinterval(PCVran)



##################################################################
### Calculate partial R2 (Nakagawa and Schielzeth 2013) - BRMS ### - NOT RIGHT YET!
##################################################################


# Partition Variance - NOT RIGHT!
fit <- b.1
# as_draws(b.1)[1][[1]]$b_Intercept # use as-draws to extract posterior samples?

# extract posterior samples
samp <- posterior_samples(fit) # all effects
samp.f <- posterior_samples(fit, "^b") # fixed effects
samp.r <- posterior_samples(fit, "^r") # random effect

# two additional parameters
setdiff(names(samp),c(names(samp.f),names(samp.r)))

samp.sd.r <- posterior_samples(fit, "^sd") # sd of random effect
samp.sd.r <- as.numeric(samp.sd.r[,1]) # make vector

# create sparse design matrix of the fixed effects
i <- rep(1:nrow(fit$data), 3); j <- rep(c(1,2,3), each=nrow(fit$data)); x <- c(rep(1,nrow(fit$data)), fit$data[,2], fit$data[,3])
X <- sparseMatrix(i, j, x=x)

# predict fitted values based on fixed effects alone
vmVarF<- numeric(length(samp.f))
for(i in 1:nrow(samp.f)){
  Var<-var(as.vector(as.numeric(samp.f[i,]) %*% t(X)))
  vmVarF[i]<-Var
}

# calculate marginal and conditional R^2
R2m <- as.mcmc(vmVarF/(vmVarF+samp.sd.r+pi^2/3)) # variance explained by fixed effects divided by total variance
R2c <- as.mcmc((vmVarF+samp.sd.r)/(vmVarF+samp.sd.r+pi^2/3)) # variance explained by fixed and random effects combined

posterior.mode(R2m); HPDinterval(R2m) # variance explained by fixed effects
posterior.mode(R2c); HPDinterval(R2c) # variance explained by the entire model (fixed + random effects)

