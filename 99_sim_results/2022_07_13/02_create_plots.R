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

n_taxa = sims$N[1]
n_sims = sims$sim %>% unique %>% length

##-----
## PGLS 
##-----

sims_summaries_pgls <- 
  sims %>% 
  select(evo,type,pgls_fit) %>% 
  rowwise() %>% 
  mutate(pgls_coefs = list(summary(pgls_fit)$coefficients),
         beta = pgls_coefs["y2","Estimate"],
         se_beta = pgls_coefs["y2","StdErr"],
         p_beta = pgls_coefs["y2","p.value"],
         sigma2 = pgls_fit$sigma2,
         lambda = pgls_fit$optpar,
         intercept = pgls_coefs["(Intercept)","Estimate"],) %>% 
  select(-starts_with("pgls")) %>% 
  pivot_longer(beta:intercept,names_to = "par_name", values_to = "par_value") %>% 
  group_by(evo,type,par_name) %>% 
  summarise(par_values = list(as_vector(par_value))) %>% 
  rowwise() %>% 
  mutate(dens = list(calc_density(par_values)),
         plot = list(plot_density(dens, par_name = par_name)))
  
# combine plots
sims_pgls_plots <- 
  sims_summaries_pgls %>% 
  group_by(evo,type) %>% 
  summarise(plot = list(plot_group(plot, evo, type)))

sims_pgls_plots$plot[[1]]

# save plots
ggarrange(plotlist = sims_pgls_plots %>% pull(plot), ncol = 1) %>% 
  annotate_figure(top = paste("PGLS --- N =",n_taxa,"#sims =",n_sims,"\n")) %>% 
  ggsave(paste0(save_dir,"/plots_pgls.pdf"),plot = ., width = 7, height = 40)


##-----
## BRMS 
##-----


# Transform data into long form
sims_brms <- 
  sims %>% 
  select(!seed:A,-brms_time, -brms_rhat)  %>% 
  pivot_longer(s2_phy_1:rho_res, names_to = "par_name", values_to = "par_value") %>% 
  rowwise() %>% 
  mutate(brms_samples = list(brms_samples %>% 
                               mutate(across(starts_with("s2"), ~ .x^2)) %>% 
                               pull(par_name)),
         post_median = list(quantile(brms_samples, probs = c(0.5)))) %>% 
  group_by(evo,type,par_name, par_value) %>% 
  summarise(post_all = list(as_vector(brms_samples)),
            post_median = list(as_vector(post_median)))

# compute summaries
sims_brms_summary <- 
  sims_brms %>% 
  rowwise() %>% 
  mutate(post_dens_all = list(calc_density(post_all)),
         post_dens_median = list(calc_density(post_median)),
         plot_dens_median = list(plot_density(post_dens_median, par_value, par_name)),
         plot_dens_all = list(plot_density(post_dens_all, par_value, par_name)))
  

# combine plots
sims_brms_plots <- 
  sims_brms_summary %>% 
  group_by(evo,type) %>% 
  summarise(plot_median = list(plot_group(plot_dens_median, evo, type, stat = "(median)")),
            plot_all = list(plot_group(plot_dens_all, evo, type, stat = "(all)")))

# save plots
ggarrange(plotlist = sims_brms_plots %>% pull(plot_median), ncol = 1) %>% 
  annotate_figure(top = paste("N =",n_taxa,"#sims =",n_sims,"\n")) %>% 
  ggsave(paste0(save_dir,"/plots_brms_median.pdf"),plot = ., width = 7, height = 40)
ggarrange(plotlist = sims_brms_plots %>% pull(plot_all), ncol = 1) %>% 
  annotate_figure(top = paste("N =",n_taxa,"#sims =",n_sims,"\n")) %>% 
  ggsave(paste0(save_dir,"/plots__brms_all.pdf"),plot = ., width = 7, height = 40)
