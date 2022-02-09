library(brms)
library(bayesplot)
library(ggpubr)

rm(list = ls())
fit <- readRDS("PMM_brms_2022_02_08.rds")

# time to fit (seconds): slowest chain of four
rstan::get_elapsed_time(fit$fit) %>% apply(1,sum) %>% max

# names of model parameters, excluding random effects and log-probability
pars <- fit %>% as_tibble %>% select(-starts_with("r_"),-"lp__") %>% names
r_pars <- fit %>% as_tibble %>% select(starts_with("r_")) %>% names

#------------------------
# convergence diagnostics
#------------------------
mcmc_trace(fit, pars = r_pars)
mcmc_rhat(rhat(fit, pars = r_pars))
mcmc_neff(neff_ratio(fit, pars = pars))

rstan::check_hmc_diagnostics(fit$fit)
rhat(fit)
rstan::check_hmc_diagnostics(fit$fit)

#----------------------------
# posterior predictive checks
#----------------------------

# prediction conditional on species-level random-intercept estimates
ggarrange(
  pp_check(fit, resp= "y1", type = "intervals") + labs(subtitle = "Response: y1"),
  pp_check(fit, resp= "y2", type = "intervals") + labs(subtitle = "Response: y2"),
  nrow = 2, legend = "none") %>% 
  annotate_figure(top = "PP-check: conditional on species-level random-intercept estimates")

# fixed effects only (i.e., ignore random effects by setting them to zero)
ggarrange(
  pp_check(fit, resp= "y1", type = "intervals", re_formula = NA) + labs(subtitle = "Response: y1"),
  pp_check(fit, resp= "y2", type = "intervals", re_formula = NA) + labs(subtitle = "Response: y2"),
  nrow = 2, legend = "none") %>% 
  annotate_figure(top = "PP-check: fixed effects only (i.e., ignore random effects by setting them to zero)")

# prediction marginalised over Gaussian distribution of random effects
pp1 <- pp_check(fit, resp= "y1", type = "intervals", 
         allow_new_levels = T, 
         sample_new_levels = "gaussian",
         newdata = fit$data %>% mutate(animal = NA))

pp2 <- pp_check(fit, resp= "y2", type = "intervals", 
                allow_new_levels = T, 
                sample_new_levels = "gaussian",
                newdata = fit$data %>% mutate(animal = NA))

ggarrange(
  pp1 + labs(subtitle = "Response: y1"),
  pp2 + labs(subtitle = "Response: y2"),
  nrow = 2, legend = "none") %>% 
  annotate_figure(top = "PP-check: marginalised over Gaussian distribution of random effects")

#------------------
# residual analysis
#------------------
pe1 <- predictive_error(fit, resp = "y1") %>% as.data.frame()
pe2 <- predictive_error(fit, resp = "y2") %>% as.data.frame()
# qq pnorm plots
pe1 %>% map_dbl(mean) %>% qqnorm()
pe2 %>% map_dbl(mean) %>% qqnorm()


#------------------------------------
# recovery of simulation parameters
#------------------------------------
x <- fit %>% as_tibble %>% select(-starts_with("r_"),-"lp__")
true_values <- c(b_y1_Intercept = 1,
                 b_y2_Intercept = 2,
                 sd_animal__y1_Intercept = 1,
                 sd_animal__y2_Intercept = 2, # sqrt needed here?
                 cor_animal__y1_Intercept__y2_Intercept = 0.75,
                 sigma_y1 = 2, # sqrt needed here?
                 sigma_y2 = 1,
                 rescor__y1__y2 = 0.25
                 )
names(true_values) == names(x)
identical(names(true_values),names(x))
mcmc_recover_hist(x, true_values)

