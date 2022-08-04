library(tidyverse)
options(scipen = 999)
#rm(list = ls())

N <- 6

#set.seed(944307) # sample(1e6,1)
tree <- ape::rtree(N); plot(tree)
A <-ape::vcv.phylo(tree, corr = T); A
I <- diag(N); I


s2_phy_1 <- 1
s2_phy_2 <- 1
s2_res_1 <- 1
s2_res_2 <- 1
rho_phy <- 0.5*sqrt(s2_phy_1*s2_phy_2)
rho_res <- 0.5*sqrt(s2_res_1*s2_res_2)


Sigma_phy <- matrix(c(s2_phy_1,rho_phy,rho_phy,s2_phy_2),2,2); Sigma_phy
Sigma_res <- matrix(c(s2_res_1,rho_res,rho_res,s2_res_2),2,2); Sigma_res
Sigma <- kronecker(Sigma_phy,A) + kronecker(Sigma_res,I); Sigma

kronecker(Sigma_phy,A)

Sigma_partition <- function(i,j) Sigma_phy[i,j]*A + Sigma_res[i,j]*I
Sigma_12 <- Sigma_partition(1,2)
Sigma_22 <- Sigma_partition(2,2)

Sigma_12
solve(Sigma_22)

beta <- Sigma_12 %*% solve(Sigma_22)

c(s2_phy_1 = s2_phy_1,
  s2_phy_2 = s2_phy_2,
  s2_res_1 = s2_res_1,
  s2_res_2 = s2_res_2,
  rho_phy = rho_phy,
  rho_res = rho_res)

beta %>% round(3)

diag(beta) %>% mean
