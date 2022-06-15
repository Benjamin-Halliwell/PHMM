library("MCMCglmm");library("brms");library("ape");
library("phytools");library("MASS");#library("plyr");
library("bindata");library('phangorn')
library("tidyverse")

rm(list = ls())

select <- dplyr::select
summarise <- dplyr::summarise
map <- purrr::map
theme_set(theme_classic())

fits <- readRDS("fits.rds")

# list nesting = [[tree.rep]][[tree.type]][[model.evo]]
# tree.rep = 1:ntrees
# tree.type = c(early split, balanced)
# model.evo = c(BM1, BM2) (PRICE to be added)

compare_fits <- function(f) f %>% 
  map_dfr(~ .x %>% as.data.frame() %>% as_tibble %>% 
            select(-starts_with(c("r","lp"))) %>% 
            imap_dfr(~ c(mean = mean(.x), err = sd(.x)),.id = "var"), .id = "m") %>% 
  ggplot(aes(y = var)) +
  geom_linerange(aes(xmin = mean-2*err, xmax = mean+2*err), data = ~ .x %>% filter(m==1)) +
  geom_point(aes(x = mean), data = ~ .x %>% filter(m==1)) +
  geom_linerange(aes(xmin = mean-2*err, xmax = mean+2*err), 
                 position = position_nudge(y = 0.1),
                 col = "blue", data = ~ .x %>% filter(m==2)) +
  geom_point(aes(x = mean), position = position_nudge(y = 0.1),
             col = "blue", data = ~ .x %>% filter(m==2)) +
  labs(x = "x")


## compare fits to BM1 and BM2 from the same tree
c(fits[[1]][[1]][1],fits[[1]][[1]][2]) %>% compare_fits + labs(subtitle = "Fits 111 and 112")

d111 <- fits[[1]][[1]][[1]]$data %>% mutate(clade = c(rep("A",nrow(.)/2),rep("B",nrow(.)/2))) %>% as.tibble
d112 <- fits[[1]][[1]][[2]]$data %>% mutate(clade = c(rep("A",nrow(.)/2),rep("B",nrow(.)/2))) %>% as.tibble
d111 %>% group_by(clade) %>% summarise(y1 =mean(y1),y2 = mean(y2)) 
d112 %>% group_by(clade) %>% summarise(y1 =mean(y1),y2 = mean(y2)) 
A111 <- fits[[1]][[1]][[2]]$data2$A
A112 <- fits[[1]][[1]][[2]]$data2$A # How do I reconstruct the tree from A?

## compare fits to BM1 from balanced and early split trees
c(fits[[1]][[1]][1],fits[[1]][[2]][1]) %>% compare_fits + labs(subtitle = "Fits 111 and 121")
summary(fits[[1]][[1]][[1]]);plot(trees[[1]][[1]]) # balanced
summary(fits[[1]][[2]][[1]]);plot(trees[[1]][[2]]) # early split

## compare fits to BM2 from balanced and early split trees
c(fits[[1]][[1]][2],fits[[1]][[2]][2]) %>% compare_fits + labs(subtitle = "Fits 112 and 122")
summary(fits[[1]][[1]][[2]])
summary(fits[[1]][[2]][[2]])




## TO DO / RESOLVE ##

# 1. WHY ARE MEANS FOR CLADES A AND B NOT MORE DIFFERENT IN OUR SIM DATA?

# means of clade A and clade B
mean(d[d$clade=="A",]$y1);mean(d[d$clade=="B",]$y1)
# sig difference between means of clade A and clade B?
t.test(d[d$clade=="A",]$y1,d[d$clade=="B",]$y1)
# compare to difference when BM simulated directly on the trees - sig2 != b_11?
x<-fastBM(trees[[1]][[1]],sig2=1)
mean(x[1:50])-mean(x[51:100])
y<-fastBM(trees[[1]][[2]],sig2=1)
mean(y[1:50])-mean(y[51:100])



##--------------------------- NOTES --------------------------------##

# final tree heights and plot first pair of trees in list
nodeheight(trees[[1]][[1]],1);plot(trees[[1]][[1]])
nodeheight(trees[[1]][[2]],1);plot(trees[[1]][[2]])

# # match heights of clades by adding difference to shorter of the pair
# x <- c(nodeheight(t1,1)-nodeheight(t1,n+1),nodeheight(t1,n+1)-nodeheight(t1,1))
# if (which(x<=0)==1){t1$edge.length[1] <- t1$edge.length[1]+x[which(x>=0)]} else {t1$edge.length[100] <- t1$edge.length[100]+x[which(x>=0)]}
