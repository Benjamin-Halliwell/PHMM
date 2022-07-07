library(tidyverse)

rm(list = ls())

N <- 5

#set.seed(944307) # sample(1e6,1)
tree <- ape::rtree(N); plot(tree)
A <- vcv.phylo(tree, corr = T); A
I <- diag(N); I

rho_phy <- 0.1
rho_res <- 0.5
sigma_phy_1

Sigma_phy <- matrix(c(1,rho_phy,rho_phy,1),2,2); Sigma_phy
Sigma_res <- matrix(c(1,rho_res,rho_res,1),2,2); Sigma_res
Sigma <- kronecker(Sigma_phy,A) + kronecker(Sigma_res,I); Sigma

Sigma_partition <- function(i,j) Sigma_phy[i,j]*A + Sigma_res[i,j]*I
Sigma_12 <- Sigma_partition(1,2)
Sigma_22 <- Sigma_partition(2,2)

Sigma_12
solve(Sigma_22)
beta <- Sigma_12 %*% solve(Sigma_22); beta

