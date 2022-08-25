sims <- sim_data %>% 
  ungroup %>% 
  mutate(pgls_ML = fits_pgls %>% map("lambda_ML"),
         pgls_0 = fits_pgls %>% map("lambda_0"),
         pgls_1 = fits_pgls %>% map("lambda_1"))
sims



##-----
## PGLS 
##-----


pgls_ML_est <- sims %>% 
  select(evo,pgls_ML) %>% 
  rowwise() %>% 
  mutate(pgls_coefs = list(summary(pgls_ML)$coefficients),
         beta_ML = pgls_coefs["y2","Estimate"]) %>% 
  select(evo,beta_ML) %>% 
  group_by(evo) %>% 
  summarise(ll = quantile(beta_ML, probs = 0.05),
            l = quantile(beta_ML, probs = 0.25),
            m = quantile(beta_ML, probs = 0.5),
            u = quantile(beta_ML, probs = 0.75),
            uu = quantile(beta_ML, probs = 0.95)) %>% 
  mutate(par_name = "beta_ML")

pgls_1_est <- sims %>% 
  select(evo,pgls_1) %>% 
  rowwise() %>% 
  mutate(pgls_coefs = list(summary(pgls_1)$coefficients),
         beta_1 = pgls_coefs["y2","Estimate"]) %>% 
  select(evo,beta_1) %>% 
  group_by(evo) %>% 
  summarise(ll = quantile(beta_1, probs = 0.05),
            l = quantile(beta_1, probs = 0.25),
            m = quantile(beta_1, probs = 0.5),
            u = quantile(beta_1, probs = 0.75),
            uu = quantile(beta_1, probs = 0.95)) %>% 
  mutate(par_name = "beta_1")

pgls_0_est <- sims %>% 
  select(evo,pgls_0) %>% 
  rowwise() %>% 
  mutate(pgls_coefs = list(summary(pgls_0)$coefficients),
         beta_0 = pgls_coefs["y2","Estimate"]) %>% 
  select(evo,beta_0) %>% 
  group_by(evo) %>% 
  summarise(ll = quantile(beta_0, probs = 0.05),
            l = quantile(beta_0, probs = 0.25),
            m = quantile(beta_0, probs = 0.5),
            u = quantile(beta_0, probs = 0.75),
            uu = quantile(beta_0, probs = 0.95)) %>% 
  mutate(par_name = "beta_0")


true_vals <-  parameters %>% 
  select(evo, rho_phy,rho_res) %>% 
  pivot_longer(-1, names_to = "par_name", values_to = "est")


bind_rows(pgls_ML_est, pgls_1_est, pgls_0_est) %>% 
  ggplot(aes(y = par_name)) +
  geom_vline(xintercept = 0, col = "grey70", lty = "longdash") +
  geom_linerange(aes(xmin = ll, xmax = uu), col = blues9[6]) +
  geom_linerange(aes(xmin = l, xmax = u), size = 1, col = blues9[8]) +
  geom_point(aes(x = m), col = blues9[9]) +
  geom_point(aes(x = est,  col = par_name), data = true_vals, shape = "|", size = 6) +
  facet_wrap(~ evo, ncol = 1, strip.position="left") +
  theme(strip.placement = "outside",
        axis.line.y = element_blank())

parameters
