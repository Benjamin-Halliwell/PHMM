#--------------------------------------------------------------------
#
# Script to summarise model fits and generate plots
# Authors: L. Yates, B. Halliwell
#
#--------------------------------------------------------------------

library(tidyverse)
library(ggpubr)
library(phylolm)

rm(list = ls())

theme_set(theme_classic())
map <- purrr::map
select <- dplyr::select

# plot functions ---------------------------------------
plot_density <- function(post_dens, par_value = NULL, par_name = ""){
  ggplot(post_dens, aes(x,y)) +  
    geom_line() + 
    geom_vline(xintercept = par_value) +
    labs(subtitle = par_name)
}

calc_density <- function(x,adjust = 1,par_name = ""){
  if(str_starts(par_name,fixed("rho"))){
    dens <- density(x, adjust = 1, from = -1, to = 1) 
  } else{
    dens <- density(x, adjust = 1)
  }
  dens %>% {tibble(x = .$x, y = .$y)}
}

plot_group <- function(plot_list, evo, stat = "") {
  ggarrange(plotlist = plot_list) %>% 
    annotate_figure(top = paste(evo,stat))
}

#-------------------------------------------------------



# set parameters
run_date <- "2022_07_26"
save_dir <- paste0("99_sim_results/",run_date)

# load sims
sims <- readRDS(paste0(save_dir,"/sims.rds"))

n_taxa = sims$N[1]
n_sims = sims$sim %>% unique %>% length

##-----
## PGLS 
##-----


pgls_est <- sims %>% 
  select(evo,pgls_fit) %>% 
  rowwise() %>% 
  mutate(pgls_coefs = list(summary(pgls_fit)$coefficients),
         beta = pgls_coefs["y2","Estimate"]) %>% 
  select(evo,beta) %>% 
  group_by(evo) %>% 
  summarise(ll = quantile(beta, probs = 0.05),
            l = quantile(beta, probs = 0.25),
            m = quantile(beta, probs = 0.5),
            u = quantile(beta, probs = 0.75),
            uu = quantile(beta, probs = 0.95)) %>% 
  mutate(par_name = "beta")



##-----
## BRMS 
##-----


brms_est <- sims %>% 
  select(!seed:A,-brms_time, -brms_rhat) %>% 
 pivot_longer(s2_phy_1:rho_res, names_to = "par_name", values_to = "par_value") %>% 
  filter(par_name == c("rho_phy","rho_res")) %>% 
  rowwise() %>% 
  mutate(est = brms_samples %>% pull(par_name) %>% median) %>% 
  group_by(evo, par_name) %>% 
  summarise(ll = quantile(est, probs = 0.05),
            l = quantile(est, probs = 0.25),
            m = quantile(est, probs = 0.5),
            u = quantile(est, probs = 0.75),
            uu = quantile(est, probs = 0.95))



# main plot

parameters <- 
  tribble(
    ~ evo, ~model, ~s2_phy_1, ~s2_phy_2, ~rho_phy, ~s2_res_1, ~s2_res_2, ~rho_res,
    "BM1", 1, 1 , 0 , 0, 1, 1, 0.5,
    "BM2", 2,1,1,0,1,1,0.5,
    "BM3",3,1,1,0.5,1,1,0,
    "BM4",4,1,1,-0.5,1,1,0.5
  )

parameters

true_vals <-  true_vals %>% 
  select(evo, rho_phy,rho_res) %>% 
  pivot_longer(-1, names_to = "par_name", values_to = "est")


bind_rows(brms_est, pgls_est) %>% 
  ggplot(aes(y = par_name)) +
  geom_vline(xintercept = 0, col = "grey70", lty = "longdash") +
  geom_linerange(aes(xmin = ll, xmax = uu), col = blues9[6]) +
  geom_linerange(aes(xmin = l, xmax = u), size = 1, col = blues9[8]) +
  geom_point(aes(x = m), col = blues9[9]) +
  geom_point(aes(x = est), data = true_vals, shape = "|", size = 6) +
  facet_wrap(~ evo, ncol = 1, strip.position="left") +
  theme(strip.placement = "outside",
        axis.line.y = element_blank())

  
  


  
  
  
  
  