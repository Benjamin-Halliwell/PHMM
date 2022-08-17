library(tidyverse); library(ggplot2); library(knitr); library(brms)
library(rstan); library(geiger);library(ape);library(ape);library(caper)
library(ggtree)

tree <- read.tree(text='(((1,2)),(3,((4,5),6)));')

tree2 <- rcoal(6)
plot(tree2)

vcv <- read.tree(text='((1:2,2:2,3:2):3,(4:2,5:2,6:2):3);')
tree <- read.tree(text='((1:2.2,2:2.2):2.8,(3:3.25,((4:1.5,5:1.5):1.5,6:3):0.25):1.75);')

vcv(vcv)
vcv(tree)-1

# re-order tips for ggtree
vcv <- read.tree(text='((6:2,5:2,4:2):3,(3:2,2:2,1:2):3);')
tree <- read.tree(text='((6:2.2,5:2.2):2.8,(((2:1.5,1:1.5):1.5,3:3):0.3,4:3.3):1.7);')

vcv(vcv)
vcv(tree)

ggtree(tree) + geom_tiplab()

tree$tip.label <- as.character(c(1,2,5,6,4,3))

p <- ggtree(vcv, size = 1) + theme_tree2() + geom_tiplab(size = 5, offset = 0.1) + #scale_x_continuous(labels=abs(5:0)) +
  theme(axis.line.x = element_line(size = 1),
        axis.ticks.x = element_line(size = 1),
        axis.text.x = element_text(size = 14),
        axis.title.x = element_text(size = 16)) +
        labs(x ="similarity") +
    geom_segment(color="red", x=3, xend=5, y=1, yend=1, size=1) +
    geom_segment(color="red", x=3, xend=5, y=2, yend=2, size=1) +
    geom_segment(color="red", x=3, xend=5, y=3, yend=3, size=1) + 
    geom_segment(color="red", x=3, xend=5, y=4, yend=4, size=1) +
    geom_segment(color="red", x=3, xend=5, y=5, yend=5, size=1) +
    geom_segment(color="red", x=3, xend=5, y=6, yend=6, size=1) 
p

p2 <- ggtree(tree, size=1) + theme_tree2() + geom_tiplab(size = 5, offset = 0.1) + #scale_x_continuous(labels=abs(5:0)) +
  theme(axis.line.x = element_line(size = 1),
        axis.ticks.x = element_line(size = 1),
        axis.text.x = element_text(size = 14),
        axis.title.x = element_text(size = 16)) +
        labs(x ="similarity") +
     geom_segment(color="red", x=4, xend=5, y=1, yend=1, size=1) +
     geom_segment(color="red", x=4, xend=5, y=2, yend=2, size=1) +
     geom_segment(color="red", x=4, xend=5, y=3, yend=3, size=1) + 
     geom_segment(color="red", x=4, xend=5, y=4, yend=4, size=1) +
     geom_segment(color="red", x=4, xend=5, y=5, yend=5, size=1) +
     geom_segment(color="red", x=4, xend=5, y=6, yend=6, size=1)

p
p2

