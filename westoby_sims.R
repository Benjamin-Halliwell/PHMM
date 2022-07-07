library(brms);library(tidyverse);library(ape);library(caper);
library(phytools);library(MASS) #library(bindata)
library(phangorn);library(mvtnorm);library(mixtools)
library(future)
rm(list = ls())

map <- purrr::map

source("functions.R")

N = 5 # number of species
Sigma <- matrix(c(1,0.75,0.75,1),2,2) # parameters for mvnorm
reps = 2
clades = 2 # fixed (generate 2 clades to stitch together---long and short for each rep)


evo = c("BM1","BM2","BM3","BM4","Price")
type = c("long","short")


# Notes
# 1) add short-long variation to Price i.e., start with close/wide separation of initial species

parameters <- expand_grid(evo,type) %>% 
  mutate(N = N,
         model = 1:n(),
         s2_phy_1 = c(rep(0.75,8),NA,NA),
               s2_phy_2 = rep(c(0,0.75,0.75,0.75,NA),each =2),
               rho_phy = rep(c(0,0,0.5,0.5,NA),each=2),
               s2_res_1 = c(rep(0.25,8),NA,NA),
               s2_res_2 = c(rep(0.25,8),NA,NA),
               rho_res = rep(c(0,0,0,0.5,NA),each=2),
               sim = 1,
               seed = 123)


sim_data <- parameters %>% 
  rowwise %>% 
  mutate(pr_vcv = list(get_price_vcv()),
         tree = list(get_tree(N,pr_vcv,type,seed)),
         trait = list(get_trait(N,evo,B = make_vcv(s2_phy_1,s2_phy_2,rho_phy), 
                                C = make_vcv(s2_res_1,s2_res_2,rho_res), 
                                pr_vcv,tree,seed)),
         A = list(calc_A(tree)))


plan(multisession(workers = 30))
fits_brms <- furrr::future_map2(sim_data$A, sim_data$trait, fit_brms, file = "m.brms")
fits_pgls <- furrr::future_map2(sim_data$A, sim_data$trait, fit_pgls)


sims <- sim_data %>% 
  rowwise %>% 
  mutate(m.brms = list(fits_brms[[model]]),
         m.pgls = list(fits_pgls[[model]]))

