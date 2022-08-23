
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

fit_brms <- function(A, trait, brms_model, cores = 1, chains = 2, iter = 3000, future = F) {

  fit = update(brms_model, newdata = trait, data2 = list(A = A),
               cores = cores, chains = chains, iter = iter, future = future)
  
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
  
  # compute rhat statistic
  rhat_est = rhat(fit)[pars]
  names(rhat_est) <- names(post)
  time_elapsed = rstan::get_elapsed_time(fit$fit) %>% apply(1,sum) %>% max
  rstan::nlist(post,rhat_est,time_elapsed)
}

fit_pgls <- function(tree,trait) {

  trait$animal <- tree$tip.label
  comp <- comparative.data(tree, trait, animal, vcv=TRUE, na.omit = F)
  # comp$vcv <- comp$vcv + diag(1e-6,nrow(comp$vcv))
  res <- list(lambda_ML = pgls(y1 ~ y2, data = comp, lambda = "ML"),
              lambda_1 = pgls(y1 ~ y2, data = comp, lambda = 1),
              lambda_0 = pgls(y1 ~ y2, data = comp, lambda = 1e-06))
  
  ## PHYLOLM()
  # rownames(trait) <- tree$tip.label
  # res <- list(phylolm(y1 ~ y2, phy=tree, model = "lambda", data=trait),
  # phylolm(y1 ~ y2, phy=tree, model = "BM", data=trait),
  # phylolm(y1 ~ y2, phy=tree, model = "lambda", lower.bound = 1, upper.bound = 1, data=trait))
  # 
  return(res)
  
}

calc_A <- function(tree, eps = 1e-6){
  A <- vcv.phylo(tree, corr = T) 
  A + diag(eps,nrow(A))
}


# #test of sim_price() function
# r.scalar = 2.125
# Sigma <- matrix(c(1,0.5,0.5,1),2,2); Sigma
# d <- sim_Price(30, Sigma, trait = T, r.scalar=r.scalar)
# plot(d)
# 
# # ellipse from runif_on_ellipsoid() and ellipse() do not match
# plot(d[,1],d[,2], xlim=c(-5,5), ylim=c(-5,5));ellipse(0, Sigma)
# # sample border of ellipse instead to confirm matches with ellipse(0, Sigma)
# # ellipse sampled in sim_price too small, need to adjust radius in runif_in_ellipsoid() to match diag(Sigma) but how?
# Sigma.inv <- Sigma*-1;diag(Sigma.inv) <- 1;d2 <- uniformly::runif_on_ellipsoid(100, Sigma.inv, Sigma[1,1]*r.scalar)
# plot(d2[,1],d2[,2], xlim=c(-5,5), ylim=c(-5,5));ellipse(0, Sigma)

# simulate bivariate data under the adaptive radiation model of Price 1997
# set limits on min.dist based on radius of ellipse? upper limit around 0.35 with radius = 1
sim_Price <- function (N, Sigma, trait = F, min.dist = 0.1, r.scalar = 2.125) {
    
    # sample uniformly from an ellipse defined by Sigma
    # radius of ellipse defined as diagonal element of Sigma but will
    # need to change if we want to  accommodate unequal variances between traits
    Sigma <- Sigma*-1;diag(Sigma) <- 1 # very strange, but runif_in_ellipsoid() produces an ellipse with the wrong sign! hack to fix but must be a better way
    
    # sample first species niche
    niche.space <- uniformly::runif_in_ellipsoid(1, Sigma, Sigma[1,1]*r.scalar) # radius != Sigma[1,1] but will do for now
    
    # add one new niche at a time subject to the distance condition min.dist
    while (nrow(niche.space) < N){
      new.niche <- uniformly::runif_in_ellipsoid(1, Sigma, Sigma[1,1]*r.scalar)
      temp.niche.space <- rbind(niche.space, new.niche)
      temp.niche.dist <- as.matrix(dist(temp.niche.space))
      if (min(temp.niche.dist[upper.tri(temp.niche.dist)]) > min.dist) {
        niche.space <- temp.niche.space
      } else NULL
    }
    
    # calculate distance matrix on final niche matrix
    niche.dist <- as.matrix(dist(niche.space))
    
    # # WHY DOESNT THIS WORK!?
    # # ensure a minimum euclidean distance between candidate phenotypes in niche space (n.dist)
    # while (min(niche.dist[upper.tri(niche.dist)]) < n.dist){
    #   niche.space <- uniformly::runif_in_ellipsoid(N, Sigma, Sigma[1,1])
    #   niche.dist <- as.matrix(dist(niche.space))
    #   print(min(niche.dist[upper.tri(niche.dist)]) < n.dist)
    # }
    
    ## PREVIOUS SAMPLING METHOD
    # # randomly generate 2D points from a multivariate normal using Sigma as VCV
    # # matrix and set mean to 0
    # niche.space <- mvrnorm(n = N, rep(0, 2), Sigma)
    # niche.dist <- as.matrix(dist(niche.space)) # euclidean distance between N in niche space
    
    
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
    
    # grow tips once more to prevent zero length terminal sisters
    growth <- rep(0, nrow(atree$edge))
    growth[atree$edge[,2]<=num.tips] <-rep(rexp(1), num.tips)
    atree$edge.length <- atree$edge.length + growth
    
    res <- list(niche.space, atree)
    if(trait) return(niche.space)
    return(atree)
    
}

# seed = 123
# pr_vcv = matrix(c(1,0.5,0.5,1),2,2)
# tree = get_tree(N,pr_vcv,seed)
# plot(tree)

get_tree <- function(N,seed){ # pr_vcv removed from arguments for sim_Price derived tree
  seed_save <- .Random.seed
  set.seed(seed)
  # phy <- rtree(2); phy$edge.length <- c(1,1)
  # phy1 <- sim_Price(N,pr_vcv)
  # b <- sim_Price(N,pr_vcv)
  # a$edge.length<-a$edge.length/max(nodeHeights(a)[,2])*0.9
  # b$edge.length<-b$edge.length/max(nodeHeights(b)[,2])*0.9
  # len <- ifelse(type=="long",1.1,0.1)
  # phy.1 <- phy;phy.1$edge.length<-phy.1$edge.length/max(nodeHeights(phy.1)[,2])*len
  # phy1 <- bind.tree(phy.1, a, where = 1, position = 0);phy1 <- bind.tree(phy1, b, where = 1, position = 0)
  # phy1$tip.label <- paste0("t", 1:(2*N))
  # .Random.seed <- seed_save
  
  phy1 <- geiger::sim.bdtree(b=1, d=0, stop=c("taxa"), n=N, seed=0, extinct=F)
  # phy1 <- ape::rcoal(n = N)
  phy1$edge.length <- phy1$edge.length/max(phytools::nodeHeights(phy1))
  phy1
  

}

sim_Price_trait <- function(N,pr_vcv,seed){
  set.seed(seed)
  a <- sim_Price(N,pr_vcv,T)
  # b <- sim_Price(N,pr_vcv,T)
  colnames(a) <- c("y1","y2")
  data.frame(animal = paste0("t",1:N), a)
}

get_trait <- function(N,evo,B, C, pr_vcv,tree,seed){
  if(evo == "Price") sim_Price_trait(N,pr_vcv,seed)
  else sim_BM_trait(B,C,tree,seed)
}

make_vcv <- function(sig2_1, sig2_2, rho){
  matrix(c(sig2_1, sqrt(sig2_1)*sqrt(sig2_2)*rho,sqrt(sig2_1)*sqrt(sig2_2)*rho,sig2_2),2,2)
}

get_price_vcv <- function(s1_p = NA,s2_p = NA,rho_p= NA,s1_r= NA,s2_r= NA,rho_r= NA, rho_fixed = NA) matrix(c(1,rho_fixed,rho_fixed,1),2,2) 
