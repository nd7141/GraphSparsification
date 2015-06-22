'''
This is a script to manipulate different data format in order to create standartizied data format for all other purposes.
The standard is considered a graph as an edge list with probabilities,
where every line is an edge (either directed or undirected) u v p, where p is probability on the edge.
'''
from __future__ import division
__author__ = 'sivanov'

import pandas as pd
import networkx as nx

def _check_format(filename, sep=" ", ncols=3):
    df = pd.read_csv(filename, sep=sep)

    # check number of columns
    if len(df.columns) != ncols:
        assert ValueError, 'Number of columns should be %s (you may specify correct separator)' %(ncols)

    # check types of colunms
    types = df.dtypes
    for num in range(ncols):
        if types[num] != 'float':
            assert ValueError, '%s column should be numeric (float), got %s instead' %(num, types[num])

    print 'Graph "%s" has correct format' %filename

def make_directed(filename, output):
    edges = dict()
    with open(filename) as f:
        for line in f:
            d = line.split()
            edges[(d[0], d[1])] = d[2]
            edges[(d[1], d[0])] = d[2]
    with open(output, 'w+') as f:
        for edge in edges:
            f.write('%s %s %s\n' %(edge[0], edge[1], edges[edge]))

def read_graph(filename, directed=False):
    if not directed:
        G = nx.Graph()
    else:
        G = nx.DiGraph()
    with open(filename) as f:
        for line in f:
            d = line.split()
            G.add_edge(int(d[0]), int(d[1]), weight=float(d[2]))
    return G

def write_graph(G, filename, directed=False):
    if not directed:
        E = dict()
        with open(filename, 'w+') as f:
            for node in G:
                for u in G[node]:
                    if (node, u) not in E and (u, node) not in E and u != node:
                        p = G[node][u]['weight']
                        E[(u, node)] = p
                        E[(node, u)] = p
                        f.write("%s %s %s\n" %(node, u, p))

def convert_idx(filename, output):
    old2new = dict()
    count = 0
    with open(filename) as f:
        with open(output, 'w+') as g:
            for line in f:
                d = line.split()
                u = int(d[0])
                v = int(d[1])
                p = float(d[2])
                if u not in old2new:
                    old2new[u] = count
                    count += 1
                if v not in old2new:
                    old2new[v] = count
                    count += 1
                g.write('%s %s %s\n' %(old2new[u], old2new[v], p))

if __name__ == "__main__":
    # df = pd.read_csv('Datasets/hep.txt', sep=' ')

    # -1. Check correctness of data format
    # 0. Obtain probabilities
    # 1. Take file with probabilities
    # 2. Write graph
    # 3. Convert nodes to correct form
    # 4. Make graph directed

    _check_format('Datasets/hep.txt', ncols=2)
    _check_format('Datasets/hep_prob.txt', ncols=3)

    G = read_graph('Datasets/hep_prob.txt')
    print len(G), len(G.edges()), 2*len(G.edges())

    G = read_graph('Datasets/hep/graph_ic.inf', True)
    print len(G), len(G.edges())

    console = []