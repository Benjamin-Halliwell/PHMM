#--------------------------------------------------------------------
#
# Script to summarise model fits and generate plots
# Authors: L. Yates, B. Halliwell
#
#--------------------------------------------------------------------

library(tidyverse)
library(ggpubr)

rm(list = ls())

theme_set(theme_classic())
map <- purrr::map
select <- dplyr::select

# plot functions ---------------------------------------

plot_density <- function(post_dens, par_value, par_name){
  ggplot(post_dens, aes(x,y)) +  
    geom_line() + 
    geom_vline(xintercept = par_value) +
    labs(subtitle = par_name)
}

calc_density <- function(x,adjust = 1,par_name = NULL){
  density(x, adjust = 1) %>% {tibble(x = .$x, y = .$y)}
}

plot_group <- function(plot_list, evo, type, stat = "") {
  ggarrange(plotlist = plot_list) %>% 
    annotate_figure(top = paste(evo, type,stat))
}

#-------------------------------------------------------



# set parameters
run_date <- "2022_07_13"
save_dir <- paste0("99_sim_results/",run_date)

# load sims
sims <- readRDS(paste0(save_dir,"/sims.rds"))

# Transform data into long form
sims2 <- 
  sims %>% 
  select(!seed:A,-brms_time, -brms_rhat)  %>% 
  pivot_longer(s2_phy_1:rho_res, names_to = "par_name", values_to = "par_value") %>% 
  rowwise() %>% 
  mutate(brms_samples = list(brms_samples %>% pull(par_name)),
         post_median = list(quantile(brms_samples, probs = c(0.5)))) %>% 
  group_by(evo,type,par_name, par_value) %>% 
  summarise(post_all = list(as_vector(brms_samples)),
            post_median = list(as_vector(post_median)))


# compute summaries
sims_summary <- 
  sims2 %>% 
  rowwise() %>% 
  mutate(post_dens_all = list(calc_density(post_all)),
         post_dens_median = list(calc_density(post_median)),
         plot_dens_median = list(plot_density(post_dens_median, par_value, par_name)),
         plot_dens_all = list(plot_density(post_dens_all, par_value, par_name)))
  

# combine plots
sims_plots <- 
  sims_summary %>% 
  group_by(evo,type) %>% 
  summarise(plot_median = list(plot_group(plot_dens_median, evo, type, stat = "(median)")),
            plot_all = list(plot_group(plot_dens_all, evo, type, stat = "(all)")))

# save plots
n_taxa = sims$N[1]
n_sims = sims$sim %>% unique %>% length
ggarrange(plotlist = sims_plots %>% pull(plot_median), ncol = 1) %>% 
  annotate_figure(top = paste("N =",n_taxa,"#sims =",n_sims,"\n")) %>% 
  ggsave(paste0(save_dir,"/plots_median.pdf"),plot = ., width = 7, height = 40)
ggarrange(plotlist = sims_plots %>% pull(plot_all), ncol = 1) %>% 
  annotate_figure(top = paste("N =",n_taxa,"#sims =",n_sims,"\n")) %>% 
  ggsave(paste0(save_dir,"/plots_all.pdf"),plot = ., width = 7, height = 40)
