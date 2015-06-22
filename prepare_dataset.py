'''
This is a script to manipulate different data format in order to create standartizied data format for all other purposes.
The standard is considered a graph as an edge list with probabilities,
where every line is an edge (either directed or undirected) u v p, where p is probability on the edge.
'''
from __future__ import division
__author__ = 'sivanov'

import pandas as pd

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


if __name__ == "__main__":
    # df = pd.read_csv('Datasets/hep.txt', sep=' ')

    _check_format('Datasets/hep.txt', ncols=2)
    _check_format('Datasets/hep_prob.txt', ncols=3)

    make_directed('Datasets/hep_prob.txt', 'Datasets/hep_dir_prob.txt')

    console = []