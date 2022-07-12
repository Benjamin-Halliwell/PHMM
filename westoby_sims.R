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
library(brms);library(ape);library(caper);
library(phytools);library(MASS) #library(bindata)
library(phangorn);library(mvtnorm);library(mixtools)
library(future)
rm(list = ls())

theme_set(theme_classic())
map <- purrr::map
select <- dplyr::select

source("functions.R")


# set parameters
N = 50 # number of taxa
n_sims = 50 
rho_fixed = 0.75
random_seed <- 3587 # sample(1e4,1)
run_date <- "2022_07_13"
save_dir <- paste0("results/",run_date)
save_run <- F
fit_models <- F
if(!dir.exists(save_dir)) dir.create(save_dir)

# set evolutionary models
evo = c("BM1","BM2","BM3","BM4","Price")
type = c("long","short")


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
  plan(multisession(workers = 30))
  
  fits_brms <- furrr::future_map2(sim_data$A, sim_data$trait, fit_brms, fit = "m.brms")
  if(save_run) saveRDS(fits_brms, paste0(save_dir,"/fits_brms.rds"))
  
  fits_pgls <- furrr::future_map2(sim_data$tree, sim_data$trait, fit_pgls)
  if(save_run) saveRDS(fits_pgls, paste0(save_dir,"/fits_pgls.rds"))
} else {
  fits_brms <- readRDS(paste0(save_dir,"/fits_brms.rds"))
  fits_pgls <- readRDS(paste0(save_dir,"/fits_pgls.rds"))
}

# add model fits summaries to main tibble
sims <- sim_data %>% ungroup %>% mutate(m.brms = fits_brms, m.pgls = fits_pgls)
if(save_run) saveRDS(sim_data, paste0(save_dir,"/sims.rds"))


sims
sims %>% 
  filter(evo == "BM1", type == "short") %>% 
  pull(m.brms) %>% 
  map(~ .x %>% as.data.frame %>% 
        as_tibble() %>% 
        select(rho_phy_est = starts_with("cor_")) %>% 
        sample_n(200)) %>% 
  bind_rows() %>% 
 pull(1) %>% 
 density %>% plot




c("sd_animal__y1_Intercept","sd_animal__y2_Intercept",
  "cor_animal__y1_Intercept__y2_Intercept","sigma_y1",
  "sigma_y2","rescor__y1__y2")  


fit %>% as.data.frame %>% 
  as_tibble() %>% names

parameters

extract_sims <- function(fits){
  fits %>% 
    map(~ .x %>% as.data.frame %>% 
          as_tibble() %>% 
          select(rho_phy = cor_animal__y1_Intercept__y2_Intercept,
                 s2_phy_1 = sd_animal__y1_Intercept,
                 s2_phy_2 = sd_animal__y2_Intercept,
                 s2_res_1 = sigma_y1,
                 s2_res_2 = sigma_y2,
                 rho_res = rescor__y1__y2
                 ) %>% 
          sample_n(200)) %>% bind_rows %>% list
}

sims_est <- sims %>% 
  group_by(evo,type) %>% 
  summarise(est = m.brms %>% extract_sims) %>% 
  left_join(parameters, by = c("evo","type")) 

sims_est %>% 
  filter(evo != "Price") #%>% 
  rowwise() %>% 
  mutate(panel_plot = list(make_plot(est,s2_phy_1,s2_phy_2,s2_res_1,s2_res_2,rho_res)))


make_plot <- function(est,rho_phy,s2_phy_1,s2_phy_2,s2_res_1,s2_res_2,rho_res){
  xvals = rstan::nlist(rho_phy,s2_phy_1,s2_phy_2,s2_res_1,s2_res_2,rho_res)
  est %>% 
    imap(~ .x %>% density %>% 
           {tibble(x = .$x, y = .$y)} %>% 
           ggplot(aes(x,y)) + geom_line() + labs(subtitle = .y) + 
           geom_vline(xintercept = xvals[.y])) %>% 
    ggpubr::ggarrange(plotlist = .)
}


plots <- lapply(sims_est$est, function(x) x %>% 
  imap(~ .x %>% density %>% 
        {tibble(x = .$x, y = .$y)} %>% 
        ggplot(aes(x,y)) + geom_line() + labs(subtitle = .y)) %>% 
  ggpubr::ggarrange(plotlist = .))

plots[[2]]
