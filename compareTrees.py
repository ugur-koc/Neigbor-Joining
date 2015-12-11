#!/opt/local/bin/python2.7

# This converts output tree of our NJ implementation into newick format
# and also generates the accuracy score for given two phylogenetic trees
# input: tree 1 in newick format (.newick file) (this tree will be treaded as groung truh)
#        tree 2 in our putput format (.txt file contains our NJ implementation putput)
# output: it will print accuracy score for tree 2 and also tree 2  in newik format
# @auhtor Khanh, Ugur Koc
# @date Dec.8.2015

from __future__ import division
from ete2 import Tree
import networkx as nx

import matplotlib.pyplot as plot
import sys

def dfs(g, u):
    global count
    count += 1
    for v in u.children:
        if (v.name == ""):
            v.name = "New_" + str(len(g))
        g.add_edge(u.name, v.name, weight=v.dist)
        dfs(g, v)

def newickToGraph(tree, format=3):
    g = nx.Graph()
    dfs(g, tree)
    return g

def parseFile(filename, format=3):
   if (format == 2):
      g = myOutputToGraph(filename)
      treeStr = networkxToNewick(g)
   else:
      treeStr = ""
      with open(filename, "r") as f:
         for line in f:
            treeStr += line.rstrip()
   return treeStr

def stringToNewick(treeStr):
   t = Tree(treeStr)
   if (t.name == ""):
      t.name = "New_0"
   return t

def myOutputToGraph(filename):
    g = nx.Graph()
    with open(filename, "r") as f:
        f.readline()
        f.readline()
        for line in f:
            tokens = line.rstrip().split()
            g.add_edge(tokens[0], tokens[1], weight=float(tokens[2]))
    return g

def traverse(g, u, p):
    numChild = 0
    for v in g.neighbors(u):
        if (v == p):
            continue
        numChild += 1
    if (numChild == 0):
        return u + ":" + str(g[u][p]['weight'])
    childStr = []
    for v in g.neighbors(u):
        if (v == p):
            continue
        childStr.append(traverse(g, v, u))
    ans = "(" + ",".join(childStr) + ")"
    if (p == "?"):
        ans += ";"
    else:
        ans += ":" + str(g[u][p]['weight'])
    return ans

def networkxToNewick(g):
    return traverse(g, 'New_0', '?') 

def drawGraph(g, filename):
    edges = [(u, v) for (u, v, d) in g.edges(data=True)]
    pos = nx.spring_layout(g)
    nx.draw_networkx_nodes(g, pos, node_size=700)
    nx.draw_networkx_edges(g, pos, edgelist=edges, width=6)
    nx.draw_networkx_labels(g, pos, font_size=10, font_family='sans-serif')
    plot.axis("off")
    plot.savefig(filename)

if __name__ == "__main__":
    count = 0
    file1 = sys.argv[1]
    file2 = sys.argv[2]

    t1Str = parseFile(file1)
    t1 = stringToNewick(t1Str)
    g1 = newickToGraph(t1)
    
    t2Str = parseFile(file2, 2)
    t2 = stringToNewick(t2Str)
    g2 = newickToGraph(t2)

    temp=t2.compare(t1,use_collateral=False, min_support_source=0.0,
                    min_support_ref=0.0, has_duplications=False, expand_polytomies=False, unrooted=True,
                    max_treeko_splits_to_be_artifact=1000)
    print "source_edges_in_ref:%0.3f,norm_rf:%0.3f" % (temp['source_edges_in_ref'], temp['norm_rf'])
    print t2Str
    #Do not change output format of this file, Ugur is parsing it in the experiments