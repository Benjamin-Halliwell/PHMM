library(tidyverse)
library(dplyr)
library(tidyr)
library(brms)
library(bayesplot)
library(ggpubr)
library(purrr)


rm(list = ls())
#fit <- readRDS("PMM_brms_2022_02_08.rds")
fit <- readRDS("b.1.rds")

select <- dplyr::select

# time to fit (seconds): slowest chain of four
rstan::get_elapsed_time(fit$fit) %>% apply(1,sum) %>% max

# names of model parameters, excluding random effects and log-probability
pars <- fit %>% as_tibble %>% select(-starts_with("r_"),-"lp__") %>% names
r_pars <- fit %>% as_tibble %>% select(starts_with("r_")) %>% names


#------------------------
# convergence diagnostics
#------------------------
mcmc_trace(fit, pars = pars)
mcmc_rhat(rhat(fit, pars = r_pars))
mcmc_neff(neff_ratio(fit, pars = r_pars))

rstan::check_hmc_diagnostics(fit$fit)

#----------------------------
# posterior predictive checks
#----------------------------

# prediction conditional on species-level random-intercept estimates
ppA1 <- pp_check(fit, resp= "logSLA", type = "intervals") + aes(x=reorder(factor(x), m)) + labs(x="taxon", subtitle = "SLA") + theme(axis.text.x = element_blank(), axis.ticks = element_blank())
ppA2 <- pp_check(fit, resp= "N", type = "intervals") + aes(x=reorder(factor(x), m)) + labs(x="taxon", subtitle = "N") + theme(axis.text.x = element_blank(), axis.ticks = element_blank())
ppA3 <- pp_check(fit, resp= "d13C", type = "intervals") + aes(x=reorder(factor(x), m)) + labs(x="taxon", subtitle = "d13C")  + theme(axis.text.x = element_blank(), axis.ticks = element_blank())

ggarrange(ppA1, ppA2, ppA3, nrow = 3, legend = "none") %>% 
  annotate_figure(top = "PP-check: conditional on species-level random-intercept estimates")

get_coverage <- function(data) data %>%
  transmute(cov90 = (hh >= y_obs) & (ll <= y_obs),
            cov50 = (h >= y_obs) & (l <= y_obs)) %>% 
  map_dbl(sum)

get_coverage(ppA1$data)/length(fit$data$animal)
get_coverage(ppA2$data)/length(fit$data$animal)
get_coverage(ppA3$data)/length(fit$data$animal)

## Re-order by pred mean
# fixed effects only (i.e., ignore random effects by setting them to zero)
ppB1 <- pp_check(fit, resp= "logSLA", type = "intervals", re_formula = NA) + labs(subtitle = "Response: logSLA")
ppB2 <- pp_check(fit, resp= "N", type = "intervals", re_formula = NA) + labs(subtitle = "Response: N")
ppB3 <- pp_check(fit, resp= "d13C", type = "intervals", re_formula = NA) + labs(subtitle = "Response: d13C")

ggarrange(ppB1, ppB2, ppB3, nrow = 3, legend = "none") %>% 
  annotate_figure(top = "PP-check: fixed effects only (i.e., ignore random effects by setting them to zero)")

get_coverage(ppB1$data)/length(fit$data$animal)
get_coverage(ppB2$data)/length(fit$data$animal)
get_coverage(ppB3$data)/length(fit$data$animal)

# prediction marginalised over Gaussian distribution of random effects
pp1 <- pp_check(fit, resp= "logSLA", type = "intervals", 
         allow_new_levels = T, 
         sample_new_levels = "gaussian",
         newdata = fit$data %>% mutate(animal = NA))

pp2 <- pp_check(fit, resp= "N", type = "intervals", 
                allow_new_levels = T, 
                sample_new_levels = "gaussian",
                newdata = fit$data %>% mutate(animal = NA))

pp3 <- pp_check(fit, resp= "d13C", type = "intervals", 
                allow_new_levels = T, 
                sample_new_levels = "gaussian",
                newdata = fit$data %>% mutate(animal = NA))

dat1 <- ppA1$data
dat2 <- ppA2$data
dat3 <- ppA3$data
pp1 <- pp1 + aes(x=reorder(factor(x), dat1$m)) + 
  labs(subtitle = "SLA", x ="taxon") + theme(axis.text.x = element_blank(), 
                                          axis.ticks = element_blank())
pp2 <- pp2 + aes(x=reorder(factor(x), dat2$m)) + 
  labs(subtitle = "N", x = "taxon") + theme(axis.text.x = element_blank(),
                                          axis.ticks = element_blank())
pp3 <- pp3 + aes(x=reorder(factor(x), dat3$m)) + 
  labs(subtitle =  bquote(delta^13*C), x="taxon") + theme(axis.text.x = element_blank(),
                                          axis.ticks = element_blank())

ggarrange(
  pp1,pp2,pp3,
  nrow = 3, legend = "none") #%>% 
  #annotate_figure(top = "PP-check: marginalised over Gaussian distribution of random effects")

pp1$data %>% get_coverage()/length(fit$data$animal)
pp2$data %>% get_coverage()/length(fit$data$animal)
pp3$data %>% get_coverage()/length(fit$data$animal)


bquote(SLA~"~"~delta^13*C)

pp1 <- pp_check(fit, resp= "logSLA", type = "intervals", allow_new_levels = T, sample_new_levels = "gaussian",newdata = fit$data %>% mutate(animal = NA)) +
  aes(x=reorder(factor(x), dat1$m)) + labs(x="taxon",y="log SLA") + theme(axis.title.x = element_blank(), axis.text.x = element_blank(), axis.ticks = element_blank(), axis.title.y = element_text(angle=0, vjust=0.5))
pp2 <- pp_check(fit, resp= "N", type = "intervals", allow_new_levels = T, sample_new_levels = "gaussian",newdata = fit$data %>% mutate(animal = NA)) +
  ylim(4,29) + aes(x=reorder(factor(x), dat2$m)) + labs(x="taxon",y="N") + theme(axis.title.x = element_blank(), axis.text.x = element_blank(), axis.ticks = element_blank(), axis.title.y = element_text(angle=0, vjust=0.5))
pp3 <- pp_check(fit, resp= "d13C", type = "intervals", allow_new_levels = T, sample_new_levels = "gaussian",newdata = fit$data %>% mutate(animal = NA)) +
  aes(x=reorder(factor(x), dat3$m)) + labs(x="taxon", y=bquote(delta^13*C)) + theme(axis.title.x = element_blank(), axis.text.x = element_blank(), axis.ticks = element_blank(), axis.title.y = element_text(angle=0, vjust=0.5))

ppA1 <- pp_check(fit, resp= "logSLA", type = "intervals") + aes(x=reorder(factor(x), m)) + labs(x="taxon") + theme(axis.title.x = element_blank(), axis.text.x = element_blank(), axis.ticks = element_blank())
ppA2 <- pp_check(fit, resp= "N", type = "intervals") + ylim(4,29) + aes(x=reorder(factor(x), m)) + labs(x="taxon") + theme(axis.title.x = element_blank(), axis.text.x = element_blank(), axis.ticks = element_blank())
ppA3 <- pp_check(fit, resp= "d13C", type = "intervals") + aes(x=reorder(factor(x), m)) + labs(x="taxon")  + theme(axis.title.x = element_blank(), axis.text.x = element_blank(), axis.ticks = element_blank())

ppA1 <- ppA1 + labs(subtitle = paste("94.8% at 90th quantile")) + theme(plot.subtitle = element_text(size = 8, vjust=-7, hjust=0.05))
ppA2 <- ppA2 + labs(subtitle = paste("85.9% at 90th quantile")) + theme(plot.subtitle = element_text(size = 8, vjust=-7, hjust=0.05))
ppA3 <- ppA3 + labs(subtitle = paste("85.6% at 90th quantile")) + theme(plot.subtitle = element_text(size = 8, vjust=-7, hjust=0.05))

pp1 <- pp1 + labs(y = "", subtitle = paste("98% at 90th quantile")) + theme(plot.subtitle = element_text(size = 8, vjust=-7, hjust=0.05))
pp2 <- pp2 + labs(y = "", subtitle = paste("97.4% at 90th quantile")) + theme(plot.subtitle = element_text(size = 8, vjust=-7, hjust=0.05))
pp3 <- pp3 + labs(y = "", subtitle = paste("98.7% at 90th quantile")) + theme(plot.subtitle = element_text(size = 8, vjust=-7, hjust=0.05))

ann1 <- ggplot() + geom_text(aes(x=0, y=0, label = "SLA"), 
            parse = TRUE, size = 5, hjust = 0.5) + theme_void()
ann2 <- ggplot() + geom_text(aes(x=0, y=0, label = "N"), 
                             parse = TRUE, size = 5, hjust = 0) + theme_void()
ann3 <- ggplot() + geom_text(aes(x=0, y=0, label = "delta^13*C"), 
                             parse = TRUE, size = 5, hjust = 0.5) + theme_void()

fig <- ggarrange(ann1,pp1,ppA1,
                 ann2,pp2,ppA2,
                 ann3,pp3,ppA3, nrow=3, ncol=3,
                 widths = c(0.4,1,1),
                 legend = "none", align = "hv")

annotate_figure(fig,
                top = text_grob("Marginal                                Conditional",
                face = "bold", size = 15, rot=0, hjust = 0.3),
                bottom = text_grob("taxon                                                    taxon",
                                   size = 12, vjust=-0.25, hjust = 0.275))



#----------------
# loo predictions
#----------------

# new approach
m.new <- fit %>% update(save_pars = save_pars(all = TRUE))
future::plan(multisession(workers = 10))
m.loo <- loo(m.new, moment_match = T, cores = 10, save_psis = T, reloo = T, future = T)
loo_predict(m, psis = m.loo$psis_object, type = "quantile", probs = c(0.025, 0.25, 0.5, 0.75, 0.975))


# draw posterior samples
ppred <- posterior_predict(fit)
ppred1 <- ppred[,,1] # logSLA
ppred2 <- ppred[,,2] # N
ppred3 <- ppred[,,3] # d13C


# compute importance sampling weights
log_ratios1 <- -1*log_lik(fit, resp = "logSLA")
log_ratios2 <- -1*log_lik(fit, resp = "N")
log_ratios3 <- -1*log_lik(fit, resp = "d13C")
r_eff1 <- loo::relative_eff(exp(-log_ratios1), chain_id = rep(1:4, each = 1000))
r_eff2 <- loo::relative_eff(exp(-log_ratios2), chain_id = rep(1:4, each = 1000))
r_eff3 <- loo::relative_eff(exp(-log_ratios3), chain_id = rep(1:4, each = 1000))
psis1 <- loo::psis(log_ratios1, cores = 10, r_eff = r_eff1)
psis2 <- loo::psis(log_ratios2, cores = 10, r_eff = r_eff2)
psis3 <- loo::psis(log_ratios3, cores = 10, r_eff = r_eff3)

loo::E_loo
# generate loo quantile predictions via PSIS
loo_pred1 <- loo::E_loo(ppred1,psis_object = psis1, log_ratios = log_ratios1, type = "quantile", 
                        probs = c(0.05,0.25,0.5,0.75,0.95))
loo_pred2 <- loo::E_loo(ppred2,psis_object = psis2, log_ratios = log_ratios2, type = "quantile", 
                        probs = c(0.05,0.25,0.5,0.75,0.95))
loo_pred3 <- loo::E_loo(ppred3,psis_object = psis3, log_ratios = log_ratios3, type = "quantile", 
                        probs = c(0.05,0.25,0.5,0.75,0.95))

# collate loo estimates
loo_pred_probs1 <- loo_pred1$value %>% t %>% as.data.frame() %>% as_tibble()
loo_pred_probs2 <- loo_pred2$value %>% t %>% as.data.frame() %>% as_tibble()
loo_pred_probs3 <- loo_pred3$value %>% t %>% as.data.frame() %>% as_tibble()
colnames(loo_pred_probs1) <- colnames(loo_pred_probs2) <- colnames(loo_pred_probs3) <- c("ll","l","m","h","hh")

# compute coverage
cov1 <- loo_pred_probs1 %>% mutate(x = 1:n(), k = loo_pred1$pareto_k, y_obs = fit$data$log_SLA) %>% 
  get_coverage() * 100/length(fit$data$animal); cov1 <- round(cov1)
cov2 <- loo_pred_probs2 %>% mutate(x = 1:n(), k = loo_pred2$pareto_k, y_obs = fit$data$N) %>% 
  get_coverage() * 100/length(fit$data$animal); cov2 <- round(cov2)
cov3 <- loo_pred_probs3 %>% mutate(x = 1:n(), k = loo_pred3$pareto_k, y_obs = fit$data$d13C) %>% 
  get_coverage() * 100/length(fit$data$animal); cov3 <- round(cov3)

# plot function
plot_loo_pred <- function(x,y_obs,cov, ylab = "y") x %>%  
  mutate(y_obs = y_obs) %>% 
  arrange(m) %>% 
  mutate(x = 1:n()) %>% 
  ggplot(aes(x)) +
  geom_point(aes(y=m), size =3, col = blues9[4]) +
  geom_linerange(aes(ymin = ll, ymax = hh), col = blues9[4]) +
  geom_linerange(aes(ymin = l, ymax = h), size = 1, col = blues9[4]) +
  geom_point(aes(y = y_obs)) +
  theme_classic() +
  labs(y = ylab, x = "taxon", subtitle = paste("Coverage: ", cov[1], "% at 90th quantile, ",cov[2],"% at 50th quantile"),
       title = "LOO-predictive checks")

# main plot
ggarrange(plot_loo_pred(loo_pred_probs1, fit$data$log_SLA, cov1, "logSLA"),
          plot_loo_pred(loo_pred_probs2, fit$data$N, cov2, "N"),
          plot_loo_pred(loo_pred_probs3, fit$data$d13C, cov3, "d13C"),
          ncol = 1)

#ggsave("loo_post_pred_2022_02_10.png", width = 8, height = 6)



#------------------
# residual analysis
#------------------
pe1 <- predictive_error(fit, resp = "logSLA") %>% as.data.frame() %>% map_dbl(mean)
pe2 <- predictive_error(fit, resp = "N") %>% as.data.frame() %>% map_dbl(mean)
pe3 <- predictive_error(fit, resp = "d13C") %>% as.data.frame() %>% map_dbl(mean)

tibble(r1 = pe1, r2 = pe2, r3 = pe3) %>% 
  ggplot() + 
  stat_qq(aes(sample = r1, col = "logSLA")) +
  stat_qq(aes(sample = r2, col = "N")) +
  stat_qq(aes(sample = r3, col = "d13C")) +
  scale_colour_manual(values = c(blues9[4], blues9[6], blues9[8])) +
  labs(col = "Response", subtitle = "QQ plots") +
  theme_classic()


# #------------------------------------
# # recovery of simulation parameters
# #------------------------------------
# x <- fit %>% as_tibble %>% select(-starts_with("r_"),-"lp__")
# true_values1 <- c(b_y1_Intercept = 1,
#                  b_y2_Intercept = 2,
#                  sd_animal__y1_Intercept = 1,
#                  sd_animal__y2_Intercept = 2, # sqrt needed here?
#                  cor_animal__y1_Intercept__y2_Intercept = 0.75,
#                  sigma_y1 = 2, # sqrt needed here?
#                  sigma_y2 = 1,
#                  rescor__y1__y2 = 0.25
#                  )
# 
# true_values2 <- c(b_y1_Intercept = 1,
#                  b_y2_Intercept = 2,
#                  sd_animal__y1_Intercept = 2,
#                  sd_animal__y2_Intercept = 4, # sqrt needed here?
#                  cor_animal__y1_Intercept__y2_Intercept = 0.75,
#                  sigma_y1 = 1, # sqrt needed here?
#                  sigma_y2 = 2,
#                  rescor__y1__y2 = 0.25
# )
# identical(names(true_values),names(x))
# mcmc_recover_hist(x, true_values2)
# 
# ## Condititonal prediction etc.
# fit
# data <- fit$data
# A <- fit$data2$A
# 
# # leave out first obs (manually)
# data_22 <- data[-1,]
# A_22 <- A[-1,-1]
# 
# fit_22 <- update(fit, newdata = data_22, data2 = list(A = A_22), cores = 4, recompile = T)
# saveRDS(fit_22, "fit_22.rds")
# fit_22
# 
# A_22_inv <- solve(A_22)
# 
# kronecker(A,diag(nrow = 2))
# 
# A_eig <- eigen(A)
# 
# Lambda <- A_eig$values %>% diag
# Lamba_inv <- (1/A_eig$values) %>% diag
# 
# Q <- A_eig$vectors
# 
# (Q) %*% Lambda %*% t(Q)  # eigendecomposition
# (Q) %*% Lambda_Inv %*% t(Q)
# 
# A <- matrix(c(1,3,3,5),2,2); A
# Sigma_phy <- matrix(c(1,0.5,0.5,1),2,2); Sigma_phy
# kronecker(t(A[1,]), Sigma_phy)
# kronecker(Sigma_phy,A)
