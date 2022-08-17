library(tidyverse)
library(expm)

rm(list = ls())

# number of traits
N <- 6

#set.seed(944307) # sample(1e6,1)
tree <- ape::rtree(N); plot(tree)
A <-ape::vcv.phylo(tree, corr = T); A
I <- diag(N); I

# set variance values
s2_phy_1 <- 1
s2_res_1 <- 1.2
s2_phy_2 <- 1
s2_res_2 <- 0.5

# set correlation values and compute covariances
rho_phy <- 0.5
cov_phy <- rho_phy*sqrt(s2_phy_1*s2_phy_2);cov_phy
rho_res <- 0.6
cov_res <- rho_res*sqrt(s2_res_1*s2_res_2);cov_res

# signal
lambda_2 =  s2_phy_2/(s2_res_2 +s2_phy_2); lambda_2
lambda_1 =  s2_phy_1/(s2_res_1 +s2_phy_1); lambda_1

# if T, change rho_phy to force scalar solution
if(T){ 
  rho_phy = (lambda_2*sqrt(s2_res_1*s2_res_2)*rho_res)/((1-lambda_2)*sqrt(s2_phy_1*s2_phy_2));rho_phy
  cov_phy <- rho_phy*sqrt(s2_phy_1*s2_phy_2)
  (1-lambda_2)*sqrt(s2_phy_1*s2_phy_2)*rho_phy - lambda_2*sqrt(s2_res_1*s2_res_2)*rho_res # check = 0
}

# build 2x2 trait VCV matrices
Sigma_phy <- matrix(c(s2_phy_1,cov_phy,cov_phy,s2_phy_2),2,2)
Sigma_res <- matrix(c(s2_res_1,cov_res,cov_res,s2_res_2),2,2)

# compute partitions
Sigma_partition <- function(i,j) Sigma_phy[i,j]*A + Sigma_res[i,j]*I
Sigma_12 <- Sigma_partition(1,2)
Sigma_22 <- Sigma_partition(2,2)

# compute coefficient
beta <- Sigma_12 %*% solve(Sigma_22)
beta %>% round(5)

(cov_phy + cov_res)/(s2_res_2+s2_phy_2)

# to do: compute beta_pgls using known (i.e., true) VCV components 



# check series expansion
if(F){
  S <- function(a,b) a*A + b*I
  
  # nth order series term
  dfdx_n <- function(X,a,b,n,k=1){
    I <- diag(nrow(X))
    coef = ((-1)^n)*(b^n)*((a+b*k)^(-n-1))
    if(n==0) return(coef*I)
    else return(coef*((A - I) %^% n))
  }
  # sum of first N series terms
  f_inv_aprx <- function(X,a,b,N,k=1) 0:N %>% map(~dfdx_n(X,a,b,.x)) %>% reduce(`+`)
  
  A-I
  (A-I) %^% 10
  
  dfdx_n(A,1.6,1,8)
  
  
  # example
  a <- 1; b <- 0.5
  solve(S(a,b))
  f_inv_aprx(A,a,b,N = 2)
  
  solve(S(a,b)) - f_inv_aprx(A,a,b,N = 4)
}





