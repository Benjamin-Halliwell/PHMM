
# simulate bivariate trait data under BM (is it really BM?)
sim_BM_trait <- function (B, C, tree, seed, m=2, beta = c(0,0)) {

# convert tree to correlation matrix
A.mat <- vcv.phylo(tree, corr = T) 

## CHANGED - preventing zero length terminal branches in sim_Price() solves this issue
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

return(d)

}

fit_brms <- function(A,trait, cores = 2, chains = 2, iter = 2000, future = F, fit = NA) {
  
  fit = brm(
    bf(mvbind(y1, y2) ~ (1|p|gr(animal, cov = A))) + set_rescor(TRUE),
    data = trait,
    data2 = list(A = A),
    family = gaussian(),
    fit = fit,
    future = future,
    chains = chains, 
    thin = 4,
    iter = iter,
    file_refit = "always")
  
  pars <- c("sd_animal__y1_Intercept",
            "sd_animal__y2_Intercept",
            "cor_animal__y1_Intercept__y2_Intercept",
            "sigma_y1",
            "sigma_y2",
            "rescor__y1__y2")
  
  # extract posterior sims for selected parameters
  post = fit %>% as.data.frame %>% 
    as_tibble() %>% 
    select(rho_phy = cor_animal__y1_Intercept__y2_Intercept,
           s2_phy_1 = sd_animal__y1_Intercept,
           s2_phy_2 = sd_animal__y2_Intercept,
           s2_res_1 = sigma_y1,
           s2_res_2 = sigma_y2,
           rho_res = rescor__y1__y2)
  
  # compure rhat statistic
  rhat_est = rhat(fit)[pars]
  names(rhat_est) <- names(post)
  list(post = post, rhat_est = rhat_est)
}

fit_pgls <- function(tree,trait) {
  ## CHANGED - phylolm a more succinct function call
  # comp <- comparative.data(tree, trait, animal, vcv=TRUE, na.omit = F)
  # comp$vcv <- comp$vcv + diag(1e-6,nrow(comp$vcv))
  # pgls(y1 ~ y2, data = comp, lambda = 1) # change to ML 
  rownames(trait) <- tree$tip.label
  phylolm(y1 ~ y2, phy=tree, model = "lambda", data=trait)
}

calc_A <- function(tree, eps = 1e-6){
  A <- vcv.phylo(tree, corr = T) 
  A + diag(eps,nrow(A))
}

# simulate bivariate data under the adaptive radiation model of Price 1997
sim_Price <- function (N, Sigma, trait = F) {
    
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
    
    ## CHANGED - grow tips once more to prevent zero length terminal sisters
    growth <- rep(0, nrow(atree$edge))
    growth[atree$edge[,2]<=num.tips] <-rep(rexp(1), num.tips)
    atree$edge.length <- atree$edge.length + growth
    
    res <- list(niche.space, atree)
    if(trait) return(niche.space)
    return(atree)
    
}

get_tree <- function(N,pr_vcv,type,seed){
  seed_save <- .Random.seed
  set.seed(seed)
  phy <- rtree(2); phy$edge.length <- c(1,1)
  a <- sim_Price(N,pr_vcv)
  b <- sim_Price(N,pr_vcv)
  a$edge.length<-a$edge.length/max(nodeHeights(a)[,2])*0.9
  b$edge.length<-b$edge.length/max(nodeHeights(b)[,2])*0.9
  len <- ifelse(type=="long",1.1,0.1)
  phy.1 <- phy;phy.1$edge.length<-phy.1$edge.length/max(nodeHeights(phy.1)[,2])*len
  phy1 <- bind.tree(phy.1, a, where = 1, position = 0);phy1 <- bind.tree(phy1, b, where = 1, position = 0)
  phy1$tip.label <- paste0("t", 1:(2*N))
  .Random.seed <- seed_save
  phy1
}

sim_Price_trait <- function(N,pr_vcv,type,seed){
  set.seed(seed)
  a <- sim_Price(N,pr_vcv,T)
  b <- sim_Price(N,pr_vcv,T)
  colnames(a) <- colnames(b) <- c("y1","y2")
  data.frame(animal = paste0("t",1:(2*N)), rbind(a,b))
}

get_trait <- function(N,evo,B, C, pr_vcv,tree,seed){
  if(evo == "Price") sim_Price_trait(N,pr_vcv,type,seed)
  else sim_BM_trait(B,C,tree,seed)
}

make_vcv <- function(sig2_1, sig2_2, rho){
  matrix(c(sig2_1, sqrt(sig2_1)*sqrt(sig2_2)*rho,sqrt(sig2_1)*sqrt(sig2_2)*rho,sig2_2),2,2)
}

get_price_vcv <- function(s1_p = NA,s2_p = NA,rho_p= NA,s1_r= NA,s2_r= NA,rho_r= NA, rho_fixed = NA) matrix(c(1,rho_fixed,rho_fixed,1),2,2) 
