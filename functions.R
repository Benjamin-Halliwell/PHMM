
# simulate bivariate trait data under BM (is it really BM?)
sim.biv.BM <- function (B, C, tree) {

# convert tree to correlation matrix
A.mat <- vcv.phylo(tree, corr = T) 

# ensure A.mat is positive definite (need to extend to ensure matrix is non-singular?)
while (dim(table(eigen(A.mat)$values > 0))==2) { # should have all eigenvalues > 0
  A.mat <-  A.mat + diag(1e-6,nrow(A.mat)) } # if there are zero values, adjust by adding a small constant to the diagonal

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

return(d)

}

# fit PMM and PGLS to bivariate data
# TO-DO
# 1. optim error for MLE of lambda on some data/trees
sim.biv.fits <- function (rep) {
  
  # extract data and tree nested within higher list element (CHANGE DEPENDING ON LIST STRUCTURE)
  d <- rep[[1]]
  tree <- rep[[2]]
  
  # convert tree to correlation matrix for brms
  A.mat <- vcv.phylo(tree, corr = T)
  # ensure A.mat is positive definite (need to extend to ensure matrix is non-singular?)
  while (dim(table(eigen(A.mat)$values > 0))==2) { # should have all eigenvalues > 0
    A.mat <-  A.mat + diag(1e-6,nrow(A.mat)) } # if there are zero values, adjust by adding a small constant to the diagonal
  
  # fit MV-PMM
  fit.pmm <- brm(
    mvbind(y1, y2) ~ (1|p|gr(animal, cov = A)),
    data = d,
    data2 = list(A = A.mat),
    family = gaussian(),
    cores = 4,
    chains = 4, iter = 1000, thin = 1)
  
  # fit PGLS
  comp <- comparative.data(tree, d, animal, vcv=TRUE, na.omit = F)
  comp$vcv <- comp$vcv + diag(1e-6,nrow(comp$vcv)) # PGLS spits error that matrix is singular. this trick OK here?
  fit.pgls <- pgls(y1 ~ y2, data = comp, lambda = 1) # really want lambda = "ML" but convergence issues for some data/trees
  
  fits <- list(fit.pgls, fit.pmm)
  
  # fits <- fit.pgls
  
  return(fits)
  
  
}
