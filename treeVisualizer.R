library(ape)
library(car)
library(plotrix)
library(ggplot2)
library(doBy)
library(plyr)
library(fmsb)
library(Hmisc)

treeM = "/Users/ugurmeryem/Dropbox/701CMSC/Neigbor-Joining/data/85VASTdomains.newick";
MyTree <- read.tree(treeM)

plot(MyTree, type = "phylogram", use.edge.length = TRUE,
     node.pos = NULL, show.tip.label = TRUE, show.node.label = FALSE,
     edge.color = "black", edge.width = 1, edge.lty = 1, font = 3,
     cex = par("cex"), adj = NULL, srt = 0, no.margin = FALSE,
     root.edge = FALSE, label.offset = 0, underscore = FALSE,
     x.lim = NULL, y.lim = NULL, direction = "rightwards",
     lab4ut = "horizontal", tip.color = "black");