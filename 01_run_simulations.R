#--------------------------------------------------------------------
#
# Main script to run simulations
# Authors: L. Yates, B. Halliwell
#
#--------------------------------------------------------------------

# Notes
# 1) add short-long variation to Price i.e., 
#    start with close/wide separation of initial species

library(dplyr)
library(tidyr)
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
library(ggplot2)
rm(list = ls())

theme_set(theme_classic())
map <- purrr::map
select <- dplyr::select

source("00_functions.R")

# set run parameters
N = 20 # number of taxa
n_sims = 2
rho_fixed = 0.5
random_seed <- 3587 # sample(1e4,1)
run_date <- "2022_08_15"
save_dir <- paste0("99_sim_results/",run_date)
save_run <- F
fit_models <- F
if(!dir.exists(save_dir)) dir.create(save_dir)

# load this into global environment to save recompilation
brms_model <- readRDS("m.brms.rds")

# set evolutionary models
evo = c("BM1","BM2","BM3","BM4","BM5")
#type = c("long","short")

set.seed(random_seed)

# set simulation parameters
parameters <- 
  tribble(
  ~ evo, ~N, ~model, ~s2_phy_1, ~s2_phy_2, ~rho_phy, ~s2_res_1, ~s2_res_2, ~rho_res,
  "BM1",N,1,0.5,0,0,0.5,1,0.7,
  "BM2",N,2,0.5,0.5,0,0.5,0.5,0.7,
  "BM3",N,3,0.5,0.5,0.7,0.5,0.5,0,
  "BM4",N,4,0.5,0.5,0.7,0.5,0.5,0.7,
  "BM5",N,5,0.5,0.5,0.7,0.5,0.5,-0.7
  #"Price",N,5,NA,NA,NA,NA,NA,NA
  )
parameters

# generate simulated trees and traits
sim_data <- 
  expand_grid(sim = 1:n_sims,parameters) %>% 
  mutate(seed = sim + random_seed) %>% 
  rowwise %>% 
  mutate(pr_vcv = list(get_price_vcv(rho_fixed = 0.5)),
         tree = list(get_tree(N,pr_vcv,seed)),
         trait = list(get_trait(N,evo,B = make_vcv(s2_phy_1,s2_phy_2,rho_phy), 
                               C = make_vcv(s2_res_1,s2_res_2,rho_res), 
                               pr_vcv,tree,seed)),
         A = list(calc_A(tree)))


if(save_run) saveRDS(sim_data, paste0(save_dir,"/sim_data.rds"))


# fit models
if(fit_models){
  plan(multisession(workers = 4))
  
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
         pgls_ML = fits_pgls %>% map("lambda_ML"),
         pgls_0 = fits_pgls %>% map("lambda_0"),
         pgls_1 = fits_pgls %>% map("lambda_1"))

if(save_run) saveRDS(sims, paste0(save_dir,"/sims.rds"))


