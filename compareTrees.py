
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

def newickToGraph(filename, format=3):
    treeStr = ""
    with open(filename, "r") as f:
        for line in f:
            treeStr += line.rstrip() 
    print treeStr
    t = Tree(treeStr)
    g = nx.Graph()
    if (t.name == ""):
        t.name = "New_0"
    dfs(g, t)
    return g

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
    file1 = sys.argv[1]
    file2 = sys.argv[2]
    g1 = newickToGraph(file1)
    g2 = myOutputToGraph(file2)
    print nx.is_isomorphic(g1, g2)
