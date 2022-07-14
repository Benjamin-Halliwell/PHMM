#--------------------------------------------------------------------
#
# Main script to run simulations
# Authors: L. Yates, B. Halliwell
#
#--------------------------------------------------------------------

# Notes
# 1) add short-long variation to Price i.e., 
#    start with close/wide separation of initial species


library(tidyverse)
library(brms)
library(ape)
library(caper)
library(phytools)
library(MASS) #library(bindata)
library(phangorn)
library(mvtnorm)
library(mixtools)
library(future)
library(phylolm)
rm(list = ls())

theme_set(theme_classic())
map <- purrr::map
select <- dplyr::select

source("00_functions.R")

# set run parameters
N = 75 # number of taxa
n_sims = 100
rho_fixed = 0.75
random_seed <- 3587 # sample(1e4,1)
run_date <- "2022_07_13"
save_dir <- paste0("99_sim_results/",run_date)
save_run <- T
fit_models <- T
if(!dir.exists(save_dir)) dir.create(save_dir)

# load this into global environment to save recompilation
brms_model <- readRDS("m.brms.rds")

# set evolutionary models
evo = c("BM1","BM2","BM3","BM4","Price")
type = c("long","short")

set.seed(random_seed)

# set simulation parameters
parameters <- expand_grid(evo,type) %>% 
  mutate(N = N,
         model = 1:n(),
         s2_phy_1 = c(rep(0.75,8),NA,NA),
               s2_phy_2 = rep(c(0,0.75,0.75,0.75,NA),each =2),
               rho_phy = rep(c(0,0,0.5,0.5,NA),each=2),
               s2_res_1 = c(rep(0.25,8),NA,NA),
               s2_res_2 = c(rep(0.25,8),NA,NA),
               rho_res = rep(c(0,0,0,0.5,NA),each=2))

# generate simulated trees and traits
sim_data <- 
  expand_grid(sim = 1:n_sims,parameters) %>% 
  mutate(seed = sim + random_seed) %>% 
  rowwise %>% 
  mutate(pr_vcv = list(get_price_vcv(rho_fixed = 0.75)),
         tree = list(get_tree(N,pr_vcv,type,seed)),
         trait = list(get_trait(N,evo,B = make_vcv(s2_phy_1,s2_phy_2,rho_phy), 
                               C = make_vcv(s2_res_1,s2_res_2,rho_res), 
                               pr_vcv,tree,seed)),
         A = list(calc_A(tree)))

if(save_run) saveRDS(sim_data, paste0(save_dir,"/sim_data.rds"))


# fit models
if(fit_models){
  plan(multisession(workers = 45))
  
  fits_brms <- furrr::future_map2(sim_data$A, sim_data$trait, fit_brms, brms_model = brms_model,.progress = T)
  if(save_run) saveRDS(fits_brms, paste0(save_dir,"/fits_brms.rds"))
  
  fits_pgls <- furrr::future_map2(sim_data$tree, sim_data$trait, fit_pgls)
  if(save_run) saveRDS(fits_pgls, paste0(save_dir,"/fits_pgls.rds"))
} else {
  fits_brms <- readRDS(paste0(save_dir,"/fits_brms.rds"))
  fits_pgls <- readRDS(paste0(save_dir,"/fits_pgls.rds"))
}


# add model fits summaries to main tibble
sims <- sim_data %>% 
  ungroup %>% 
  mutate(brms_samples = map(fits_brms,"post"),
         brms_rhat = map(fits_brms,"rhat_est"),
         brms_time = map_dbl(fits_brms,"time_elapsed"),
         pgls_fit = fits_pgls)

if(save_run) saveRDS(sims, paste0(save_dir,"/sims.rds"))


