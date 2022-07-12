#--------------------------------------------------------------------
#
# Script to summarise model fits and generate plots
# Authors: L. Yates, B. Halliwell
#
#--------------------------------------------------------------------

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

rm(list = ls())

theme_set(theme_classic())
map <- purrr::map
select <- dplyr::select

source("00_functions.R")

# set parameters
run_date <- "2022_07_13"
save_dir <- paste0("99_sim_results/",run_date)

# load sims
sims <- readRDS(sim_data, paste0(save_dir,"/sims.rds"))




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
