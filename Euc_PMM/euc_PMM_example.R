library(ips);library(ggpubr);library(tidyr);library(data.table);library(caper);library(RRphylo)
library(stringr);library(stringdist);library(phangorn);library(corrplot);library(grid)
library(RColorBrewer);library(phytools);library(dplyr);library(plyr);library(rr2);library(motmot)
library(phylolm);library(MCMCglmm);library(diversitree);library(ggtree);library(devtools)
library(BiocManager);library(tidytree);library(cowplot);library(ggpubr);library(patchwork)
library(surface);library(ouch);library(geiger);library(igraph);library(Rphylopars);library(rlist)
require(dplyr);library(brms);library(future)


### DATA AND TREES ####

# d.master <- read.csv("d.master.csv")
# d.master <- d.master[d.master$num_occ>2,]
d.master <- read.csv("d.master_V4.csv")


## TREES
# tree <- read.tree(file = 'Eucalypts_ML1_dated_r8s.tre')
tree.ML1 <- read.tree(file="tree.ML1.tre")
tree.ML2 <- read.tree(file="tree.ML2.tre")
tree.Bayes <- read.tree(file="tree.Bayes.tre")

## Sequence data from Thornhill
# seq <- read.fas("Eucalypts_final_PD.fas")

# resolve polytomies and make clock-like
tree.ML1 <- multi2di(tree.ML1, random = F)
tree.ML2 <- multi2di(tree.ML2, random = F)
tree.Bayes <- multi2di(tree.Bayes, random = F)
tree.ML1$node.label <- (length(tree.ML1$tip.label)+1):(tree.ML1$Nnode+length(tree.ML1$tip.label))
tree.ML2$node.label <- (length(tree.ML1$tip.label)+1):(tree.ML1$Nnode+length(tree.ML1$tip.label))
tree.Bayes$node.label <- (length(tree.ML1$tip.label)+1):(tree.ML1$Nnode+length(tree.ML1$tip.label))
tree.ML1.u <- force.ultrametric(tree.ML1)
tree.ML2.u <- force.ultrametric(tree.ML2)
tree.Bayes.u <- force.ultrametric(tree.Bayes)

# subset data to Bayes phy
d.bayes <- d.master[d.master$taxon %in% tree.Bayes$tip.label,]

### END ####


### PMM ####

# trees
phy1 <- tree.ML1.u
phy2 <- tree.ML2.u
phy3 <- tree.Bayes.u

# set plot view
par(mar=(c(2,2,2,2)))

# resolve zero branch lengths
for (i in 1:3){phy1$edge.length[phy1$edge.length<=0.00001] <- 0.0001
phy1 <- force.ultrametric(phy1)}

for (i in 1:3){phy2$edge.length[phy2$edge.length<=0.00001] <- 0.0001
phy2 <- force.ultrametric(phy2)}

for (i in 1:3){phy3$edge.length[phy3$edge.length<=0.00001] <- 0.0001
phy3 <- force.ultrametric(phy3)}


## SUBSET TO SYMPHYOMYRTUS AND EUCALYPTUS ##
d.subgen <- d.master[d.master$subgenus=="Symphyomyrtus"|d.master$subgenus=="Eucalyptus",]
row.names(d.subgen) <-  1:nrow(d.subgen)
d.subgen$subgenus <- as.factor(d.subgen$subgenus)


tree.subgen <- keep.tip(phy1, d.subgen$taxon)
for (i in 1:3){tree.subgen$edge.length[tree.subgen$edge.length<=0.00001] <- 0.005
tree.subgen <- force.ultrametric(tree.subgen)}



## COMPLETE CASES FOR SE
# d.sub <- d.subgen[!is.na(d.subgen$SLA.log)&                   
#                     !is.na(d.subgen$leaf_N_per_dry_mass)&
#                     !is.na(d.subgen$leaf_delta13C)&
#                     !is.na(d.subgen$se.specific_leaf_area)&
#                     !is.na(d.subgen$se.leaf_N_per_dry_mass)&
#                     !is.na(d.subgen$se.leaf_delta13C),]

## COMPLETE CASES FOR TRAITS
d.sub <- d.subgen[!is.na(d.subgen$SLA)&
                    !is.na(d.subgen$leaf_N_per_dry_mass)&
                    !is.na(d.subgen$leaf_delta13C),]

# subset data and other trees to Bayes tree
d.sub <- d.sub[d.sub$taxon %in% phy3$tip.label,]
phy1 <- keep.tip(phy1, d.sub$taxon)
phy2 <- keep.tip(phy2, d.sub$taxon)
phy3 <- keep.tip(phy3, d.sub$taxon)
names(d.sub)[names(d.sub)=="se.specific_leaf_area"] <- "se.SLA"
d.sub$obs <- 1:nrow(d.sub)

# ensure all SE are positive by adding min SE observed for trait to any that equal 0
d.sub[!is.na(d.sub$se.SLA) & d.sub$se.SLA==0,]$se.SLA <- min(d.sub[!is.na(d.sub$se.SLA) & d.sub$se.SLA!=0,]$se.SLA)
d.sub[!is.na(d.sub$se.leaf_N_per_dry_mass) & d.sub$se.leaf_N_per_dry_mass==0,]$se.leaf_N_per_dry_mass <- min(d.sub[!is.na(d.sub$se.leaf_N_per_dry_mass) & d.sub$se.leaf_N_per_dry_mass!=0,]$se.leaf_N_per_dry_mass)
d.sub[!is.na(d.sub$se.leaf_delta13C) & d.sub$se.leaf_delta13C==0,]$se.leaf_delta13C <- min(d.sub[!is.na(d.sub$se.leaf_delta13C) & d.sub$se.leaf_delta13C!=0,]$se.leaf_delta13C)

# for any NAs, make SE the 90th pctl. of SE for that trait
d.sub[is.na(d.sub$se.SLA),]$se.SLA <- quantile(d.sub[!is.na(d.sub$se.SLA),]$se.SLA, 0.9)
d.sub[is.na(d.sub$se.leaf_N_per_dry_mass),]$se.leaf_N_per_dry_mass <- quantile(d.sub[!is.na(d.sub$se.leaf_N_per_dry_mass),]$se.leaf_N_per_dry_mass, 0.9)
d.sub[is.na(d.sub$se.leaf_delta13C),]$se.leaf_delta13C <- quantile(d.sub[!is.na(d.sub$se.leaf_delta13C),]$se.leaf_delta13C, 0.9)


# log-transform SLA and compute first-order standard error estimates (i.e. relative error) via delta method (se.SLA/SLA)
d.sub <- d.sub %>% as_tibble %>% mutate(log_SLA = log(SLA), se_log_SLA = se.SLA/SLA,
                                        N = leaf_N_per_dry_mass, se_N = se.leaf_N_per_dry_mass,
                                        d13C = leaf_delta13C, se_d13C = se.leaf_delta13C,
                                        N_area = leaf_N_per_dry_mass/SLA)

# subset to relevant columns
d.sub2 <- d.sub %>% select(animal, log_SLA, se_log_SLA,N, se_N,d13C,se_d13C, N_area)


# ## SAVED FITS
b.0 <- readRDS("b.0.rds")
b.1 <- readRDS("b.1.rds")
b.2 <- readRDS("b.2.rds")
b.3 <- readRDS("b.3.rds")

## brms
A.mat <- ape::vcv.phylo(phy1, corr = T)


## MOD 0 - raw trait scale
bf_y1 <- bf(SLA | resp_se(se.SLA, sigma = TRUE) ~ 1 + (1|b|gr(animal, cov = A)))
bf_y2 <- bf(leaf_N_per_dry_mass | resp_se(se.leaf_N_per_dry_mass, sigma = TRUE) ~ 1 + (1|b|gr(animal, cov = A)))
bf_y3 <- bf(leaf_delta13C | resp_se(se.leaf_delta13C, sigma = TRUE) ~ 1 + (1|b|gr(animal, cov = A)))

b.0 <- brm(bf_y1 + bf_y2 + bf_y3 + set_rescor(TRUE),
           data = d.sub, 
           family = gaussian(), 
           data2 = list(A = A.mat),
           control=list(adapt_delta = 0.85, max_treedepth = 12),
           cores = 4, chains = 4, iter = 6000, thin = 3)
saveRDS(b.0, file = "b.0.rds")
# b.0
# plot(b.0) # posteriors for phy cor terms skewed and stuck around 0 or 1

# pp checks show SLA needs to be transform to avoid boundary issues
# p1 <- pp_check(b.0, resp = "SLA", ndraws = 100) + ggtitle("SLA") # prediction poor for SLA
# p2 <- pp_check(b.0, resp = "leafNperdrymass", ndraws = 100) + ggtitle("N")
# p3 <- pp_check(b.0, resp = "leafdelta13C", ndraws = 100) + ggtitle("d13C")
# ggarrange(p1,p2,p3, nrow = 1, ncol = 3)


# fit with log transformed SLA instead

## phy1
bf_y1 <- bf(log_SLA | resp_se(se_log_SLA, sigma = TRUE) ~ 1 + (1|b|gr(animal, cov = A))) 
bf_y2 <- bf(N | resp_se(se_N, sigma = TRUE) ~ 1 + (1|b|gr(animal, cov = A))) 
bf_y3 <- bf(d13C | resp_se(se_d13C, sigma = TRUE) ~ 1 + (1|b|gr(animal, cov = A))) 

b.1 <- brm(bf_y1 + bf_y2 + bf_y3 + set_rescor(TRUE),
            data = d.sub2, 
            family = gaussian(), 
            data2 = list(A = A.mat),
            control=list(adapt_delta = 0.85, max_treedepth = 12),
            cores = 4, chains = 4, iter = 6000, thin = 3)
saveRDS(b.1, file = "b.1.rds")
# readRDS(file = "b.1.rds")

# # pp checks now look good. COntinue fit for other trees
# p1 <- pp_check(b.1, resp = "logSLA", ndraws = 100) + ggtitle("SLA") # prediction poor for SLA
# p2 <- pp_check(b.1, resp = "N", ndraws = 100) + ggtitle("N")
# p3 <- pp_check(b.1, resp = "d13C", ndraws = 100) + ggtitle("d13C")
# ggarrange(p1,p2,p3, nrow = 1, ncol = 3)


## phy2
A.mat2 <- ape::vcv.phylo(phy2, corr = T)
b.2 <- brm(bf_y1 + bf_y2 + bf_y3 + set_rescor(TRUE),
             data = d.sub2, 
             family = gaussian(), 
             data2 = list(A = A.mat2),
             control=list(adapt_delta = 0.85, max_treedepth = 12),
             cores = 4, chains = 4, iter = 6000, thin = 3)
saveRDS(b.2, file = "b.2.rds")


## phy3
A.mat3 <- ape::vcv.phylo(phy3, corr = T)
b.3 <- brm(bf_y1 + bf_y2 + bf_y3 + set_rescor(TRUE),
             data = d.sub2, 
             family = gaussian(), 
             data2 = list(A = A.mat3),
             control=list(adapt_delta = 0.85, max_treedepth = 12),
             cores = 4, chains = 4, iter = 6000, thin = 3)
saveRDS(b.3, file = "b.3.rds")




# calculate h2 
mod <- b.1
mean(unlist(purrr::map(mod$fit@sim$samples, ~.x$sd_animal__logSLA_Intercept))/(unlist(purrr::map(mod$fit@sim$samples, ~.x$sd_animal__logSLA_Intercept))+unlist(purrr::map(mod$fit@sim$samples, ~.x$sigma_logSLA))))
quantile(unlist(purrr::map(mod$fit@sim$samples, ~.x$sd_animal__logSLA_Intercept))/(unlist(purrr::map(mod$fit@sim$samples, ~.x$sd_animal__logSLA_Intercept))+unlist(purrr::map(mod$fit@sim$samples, ~.x$sigma_logSLA))), probs = c(0.025, 0.975))


## COMBINE MODELS ##

## results similar across 3 available trees so combine
# b.mod <- combine_models(b.1, b.2, b.3)
# saveRDS(b.mod, file = "b.mod.rds")
b.mod <- readRDS(file = "b.mod.rds")


# calculate h2 (b.mod for combined chains across trees)
mod <- b.1
mean(unlist(purrr::map(mod$fit@sim$samples, ~.x$sd_animal__logSLA_Intercept))/(unlist(purrr::map(mod$fit@sim$samples, ~.x$sd_animal__logSLA_Intercept))+unlist(purrr::map(mod$fit@sim$samples, ~.x$sigma_logSLA))))
quantile(unlist(purrr::map(mod$fit@sim$samples, ~.x$sd_animal__logSLA_Intercept))/(unlist(purrr::map(mod$fit@sim$samples, ~.x$sd_animal__logSLA_Intercept))+unlist(purrr::map(mod$fit@sim$samples, ~.x$sigma_logSLA))), probs = c(0.025, 0.975))

mean(unlist(purrr::map(mod$fit@sim$samples, ~.x$sd_animal__N_Intercept))/(unlist(purrr::map(mod$fit@sim$samples, ~.x$sd_animal__N_Intercept))+unlist(purrr::map(mod$fit@sim$samples, ~.x$sigma_N))))
quantile(unlist(purrr::map(mod$fit@sim$samples, ~.x$sd_animal__N_Intercept))/(unlist(purrr::map(mod$fit@sim$samples, ~.x$sd_animal__N_Intercept))+unlist(purrr::map(mod$fit@sim$samples, ~.x$sigma_N))), probs = c(0.025, 0.975))

mean(unlist(purrr::map(mod$fit@sim$samples, ~.x$sd_animal__d13C_Intercept))/(unlist(purrr::map(mod$fit@sim$samples, ~.x$sd_animal__d13C_Intercept))+unlist(purrr::map(mod$fit@sim$samples, ~.x$sigma_d13C))))
quantile(unlist(purrr::map(mod$fit@sim$samples, ~.x$sd_animal__d13C_Intercept))/(unlist(purrr::map(mod$fit@sim$samples, ~.x$sd_animal__d13C_Intercept))+unlist(purrr::map(mod$fit@sim$samples, ~.x$sigma_d13C))), probs = c(0.025, 0.975))


# worth fitting a skew normal to SLA?
p1 <- pp_check(b.mod, resp = "logSLA", ndraws = 100) + ggtitle("SLA") # prediction poor for SLA
p2 <- pp_check(b.mod, resp = "N", ndraws = 100) + ggtitle("N")
p3 <- pp_check(b.mod, resp = "d13C", ndraws = 100) + ggtitle("d13C")
ggarrange(p1,p2,p3, nrow = 1, ncol = 3)



## PLOTTING ####

library(tidybayes)

b.mod$fit %>%
  gather_draws(cor_animal__logSLA_Intercept__N_Intercept,
               cor_animal__logSLA_Intercept__d13C_Intercept,
               cor_animal__N_Intercept__d13C_Intercept,
               rescor__logSLA__N,
               rescor__logSLA__d13C,
               rescor__N__d13C) %>%
  median_hdi(.width = c(0.95, 0.5)) %>% 
  mutate(cor = rep(rep(c("SLA ~ d13C","SLA ~ N_mass","N_mass ~ d13C"),2),2)) %>% 
  mutate(level = rep(rep(c("phy","ind"),each=3),2)) %>%
  ggplot(aes(y = cor, x = .value, xmin = .lower, xmax = .upper, col = level)) +
  geom_pointinterval(position = position_dodge(width = 0.25)) +
  geom_vline(xintercept = 0, col = "grey70", lty = "longdash") +
  theme_classic() + 
  theme(axis.text = element_text(size = 12),
        axis.title = element_text(size = 14),
        legend.text = element_text(size = 12),
        legend.title = element_text(size = 14)) +
  xlab("correlation estimate") + ylab("") + xlim(-1,1) +
  scale_color_manual(values = c("red", "#0072B2")) + 
  scale_y_discrete(labels=c("SLA ~ N_mass" = bquote(SLA~"~"~N[mass]),
                            "SLA ~ d13C" = bquote(SLA~"~"~delta^13*C),
                            "N_mass ~ d13C" = bquote(N[mass]~"~"~delta^13*C)))

b.N %>%
  gather_draws(cor_animal__logSLA_Intercept__Narea_Intercept,
               cor_animal__logSLA_Intercept__d13C_Intercept,
               cor_animal__Narea_Intercept__d13C_Intercept,
               rescor__logSLA__Narea,
               rescor__logSLA__d13C,
               rescor__Narea__d13C) %>%
  median_hdi(.width = c(0.95, 0.5)) %>% 
  mutate(cor = rep(rep(c("SLA ~ d13C","SLA ~ N_area","N_area ~ d13C"),2),2)) %>% 
  mutate(level = rep(rep(c("phy","res"),each=3),2)) %>%
  ggplot(aes(y = cor, x = .value, xmin = .lower, xmax = .upper, col = level)) +
  geom_pointinterval(position = position_dodge(width = 0.25)) +
  geom_vline(xintercept = 0, col = "grey70", lty = "longdash") +
  theme_classic() + 
  theme(axis.text = element_text(size = 12),
        axis.title = element_text(size = 14),
        legend.text = element_text(size = 12),
        legend.title = element_text(size = 14)) +
  xlab("correlation estimate") + ylab("") + xlim(-1,1) +
  scale_color_manual(values = c("black", blues9[6])) + 
  scale_y_discrete(labels=c("SLA ~ N_area" = bquote(SLA~"~"~N[area]),
                            "SLA ~ d13C" = bquote(SLA~"~"~delta^13*C),
                            "N_area ~ d13C" = bquote(N[area]~"~"~delta^13*C)))



## MULTI PANEL PLOT


reps <- 3 # fixed for now
N <- 30
rho_b = 0.9
rho_c = 0
sig_b = 1
sig_c = 1
C <- matrix(c(sig_c,rho_c,rho_c,sig_c),2,2) # residual covariance matrix
B <- matrix(c(sig_b,rho_b,rho_b,sig_b),2,2) # phylogenetic covariance matrix
C.scale <- 1; C <- C*C.scale
B.scale <- 1; B <- B*B.scale
blength <- 0.1


cols <- c("#E69F00","#56B4E9")

p1.1 <- ggplot(d.sub, aes(N,log_SLA)) + 
  geom_point(aes(col = subgenus)) + 
  theme_classic() +
  theme(legend.position = "none",
        plot.margin = margin(0.25,0.2,0.2,0.675, "cm")) +
  scale_color_manual(values=cols) +
  ylab("log SLA")
p1.2 <- ggplot(d.sub, aes(d13C,log_SLA)) + 
  geom_point(aes(col = subgenus)) + 
  theme_classic() +
  theme(legend.position = "none",
        plot.margin = margin(0.25,0.2,0.2,0.675, "cm")) +
  scale_color_manual(values=cols) +
  xlab(bquote(delta^13*C)) +
  ylab("log SLA")
p1.3 <- ggplot(d.sub, aes(N,d13C)) + 
  geom_point(aes(col = subgenus)) + 
  theme_classic() +
  theme(legend.position = "none",
        plot.margin = margin(0.25,0.2,0.2,0.25, "cm")) +
  scale_color_manual(values=cols) +
  ylab(bquote(delta^13*C))



# # 3D scatter plot
# library(plotly)
# fig <- plot_ly(d.sub, x = ~log_SLA, y = ~N, z = ~d13C, color = ~subgenus, colors = c('#BF382A', '#0C4B8E'), size = I(50))
# fig <- fig %>% add_markers()
# fig <- fig %>% layout(scene = list(xaxis = list(title = 'SLA (log)'),
#                                    yaxis = list(title = 'leaf N'),
#                                    zaxis = list(title = 'd13C')))
# fig



# PLOT 2
p2 <- ggtree(phy2)
p2 <- p2 + theme(plot.margin = margin(1,0,1.2,0.1, "cm"))
clade.lab <- unique(d.sub$subgenus)
grp <- list()
for(i in 1:length(clade.lab)){grp[[i]] = d.sub[d.sub$subgenus==clade.lab[1],]$taxon}
names(grp) <-  unique(d.sub$subgenus)

p2 <- groupOTU(p2, grp, 'clade') + aes(color=clade) + 
  theme(legend.position = c(0.4,0.975)) +
  # geom_tiplab(hjust = 0, align=TRUE, size=2) +
  scale_color_manual(values=cols,
                     name = "subgenus", labels = c("Eucalyptus", "Symphyomyrtus"))

# p2 + geom_tiplab(aes(label=node), size=2, offset = 0.05)
# p2.f <-  flip(p2, 92,151)
# p2.f + geom_tiplab(aes(label=node), size=2, offset = 0.05)

# PLOT 3
target <- get_taxa_name(p2)
d2 <- d.sub[match(target, d.sub$taxon),]
d2 <- data.frame(
  label = rep(d2$taxon,3),
  trait = rep(c("SLA","N","d13C"), each=length(phy1$tip.label)),
  value = c(scale(d2$log_SLA),scale(d2$N),scale(d2$d13C)))
d2$label <- factor(d2$label, level = get_taxa_name(p2)) # re-order to match ggtree

d2[d2$value==max(d2$value),]$value <- 5

p3 <- ggplot(d2, aes(x=trait, y=label)) + 
  geom_tile(aes(fill=value)) + 
  # scale_fill_gradientn(limits=c(-2,2),
  #                      colours = c("cyan2", "white", "orange2"),
  #                      values = scales::rescale(c(-0.5, -0.2, 0, 0.2, 0.5))) +
  scale_fill_viridis_c(name = "xyz",
                       # limits = c(-2.5, 2.5),
                       # breaks = c(-2, -1, 0, 1, 2))
  ) + theme_minimal() + xlab(NULL) + ylab(NULL) +
  scale_x_discrete(position = "top",
                   labels=c("d13C" = bquote(delta^13*C))) +
  scale_y_discrete(limits=rev(levels(d2$label))) +
  # scale_x_discrete(expand=c(0,0)) +
  theme(axis.text.y = element_blank(),
        axis.text.x = element_text(size=9),
        # axis.text.x = element_blank(),
        panel.grid = element_blank(),
        plot.margin = margin(0.4,0,1.3,0.1, "cm"),
        # plot.margin = margin(1,0,1.2,0.1, "cm"),
        legend.title = element_blank(),
        legend.key.size = unit(0.5, 'cm'))

pA <- cowplot::plot_grid(p1.1, p1.2, p1.3, axis = "l", nrow=3)
pB <- cowplot::plot_grid(p2, p3, ncol=2, rel_widths = c(1,0.75))

cowplot::plot_grid(pA, pB, nrow=1, rel_widths = c(0.65,1))


## END ####


## PREDICTION

# choose fit
pp <- predict(b.mod)

str(pp)
pSLA <- pp[,1,1]
pN <- pp[,1,2]
pd13C <- pp[,1,3]

d.sub <- d.sub %>% mutate(SLA_res = d.sub$SLA - pp[,1,1])
d.sub <- d.sub %>% mutate(N_res = d.sub$N - pp[,1,2])
d.sub <- d.sub %>% mutate(d13C_res = d.sub$d13C - pp[,1,3])

## while there is no relationship between N_mass and N_area, there is a significant
# positive relationship between the residuals of N_mass (i.e., after accounting for phylogeny) 
# and N_area (strengthened by removing outlier). Therefore, the positive residual correlation between
# N_mass and d13C partly reflects co-variation between N_area and d13C, for which we predict a positive
# relationship reflecting lower isotope discrimination with higher carboxylation capacity.
plot(d.sub$N, d.sub$N_area); cor.test(d.sub$N, d.sub$N_area);abline(lm(N_area ~ N, data = d.sub))
plot(d.sub$N_res, d.sub$N_area); cor.test(d.sub$N_res, d.sub$N_area);abline(lm(N_area ~ N_res, data = d.sub))
plot(d.sub[d.sub$N_res<10,]$N_res, d.sub[d.sub$N_res<10,]$N_area); cor.test(d.sub[d.sub$N_res<10,]$N_res, d.sub[d.sub$N_res<10,]$N_area);abline(lm(N_area ~ N_res, data = d.sub[d.sub$N_res<10,]))


plot(d.sub$N_res, d.sub$log_SLA); cor.test(d.sub$N_res, d.sub$log_SLA);abline(lm(log_SLA ~ N_res, data = d.sub))


plot(log(d.sub$SLA), log(d.sub$N));cor.test(log(d.sub$SLA), log(d.sub$N))
plot(log(d.sub$SLA), log(d.sub$N_area));cor.test(log(d.sub$SLA), log(d.sub$N_area))

## BRMS Multiple ##

# make list of data and trees
t.list <- list(tree.sub1,tree.sub2)
t.list <- lapply(t.list, ape::vcv.phylo)
t.list <- list(list(A=t.list[[1]]),list(A=t.list[[2]]))
d.list <- list(d.sub, d.sub)

# fit over trees and parallel process
plan(multiprocess)
brm_multiple(bf_y1 + bf_y2 + bf_y3 + set_rescor(TRUE),
             data = d.list,
             data2 = t.list,
             family = gaussian(), 
             control=list(adapt_delta = 0.85, max_treedepth = 12),
             cores = 4, chains = 4, iter = 6000, thin = 3)


# ggplot(d.sub, aes(x=leaf_N_per_dry_mass,y=leaf_delta13C,col=subgenus)) + geom_point() + theme_classic() + geom_smooth(method = lm)


### END ####




