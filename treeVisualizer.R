library(ape)
library(car)
library(plotrix)
library(ggplot2)
library(doBy)
library(plyr)
library(fmsb)
library(Hmisc)


target <- read.tree("/Users/ugurmeryem/Dropbox/701CMSC/Neigbor-Joining/data/85VASTdomains.tree")
current <- read.tree("/Users/ugurmeryem/Dropbox/701CMSC/Neigbor-Joining/data/our_tree.txt")


## S3 method for class 'phylo'
result<-all.equal(target, current, use.edge.length = TRUE,
          use.tip.label = TRUE, index.return = FALSE,
          tolerance = .Machine$double.eps ^ 0.5,
          scale = NULL);

STOP("XXX");

p<-plot(target, type = "phylogram", use.edge.length = TRUE,
     node.pos = NULL, show.tip.label = TRUE, show.node.label = FALSE,
     edge.color = "black", edge.width = 1, edge.lty = 1, font = 3,
     cex = par("cex"), adj = NULL, srt = 0, no.margin = FALSE,
     root.edge = FALSE, label.offset = 0, underscore = FALSE,
     x.lim = NULL, y.lim = NULL, direction = "rightwards",
     lab4ut = "horizontal", tip.color = "black");

print(p)
dev.print(device=postscript, file= paste("/Users/ugurmeryem/Desktop/tree2.eps"))
dev.off()