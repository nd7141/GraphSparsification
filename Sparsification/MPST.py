'''
Implementation of Most Probable Spanning Tree (MPST).
Extract MPST and delete edges from original graph that belong to that MPST.
Continue until all edges are deleted from G.
Use those MPST to select a possible world (PW) and approximate reliability.

We are solving sparsification problem of uncertain graphs for preserving reliability.
'''
from __future__ import division
import networkx as nx
from itertools import cycle
from math import exp, log
import random, time, sys, json, os
from collections import Counter
from itertools import product, izip
import matplotlib.pyplot as plt
import matplotlib
from memory_profiler import memory_usage
from shutil import copyfile, rmtree
from subprocess import call

def _comprehension_flatten(iter_lst):
    return [item for lst in iter_lst for item in lst]

def minimum_spanning_edges(G,weight='weight',data=True):
    from networkx.utils import UnionFind
    if G.is_directed():
        raise nx.NetworkXError(
            "Mimimum spanning tree not defined for directed graphs.")

    subtrees = UnionFind()
    edges = sorted(G.edges(data=True),key=lambda t: t[2].get(weight,1))
    for u,v,d in edges:
        if subtrees[u] != subtrees[v]:
            if data:
                yield (u,v,d)
            else:
                yield (u,v)
            subtrees.union(u,v)

# adopted MST algo without isolated nodes
# solution found here http://networkx.lanl.gov/_modules/networkx/algorithms/mst.html
def minimum_spanning_tree(G,weight='weight'):
    T=nx.Graph(nx.minimum_spanning_edges(G,weight=weight,data=True))
    return T

def get_MPSTs(G):
    '''
    Expected to have -log(p_e) as a weight for every edge e in G.
    '''
    assert type(G) == type(nx.Graph()), "Graph should be undirected"

    P = sum([exp(1)**(-data["weight"]) for (u,v,data) in G.edges(data=True)]) # expected number of edges
    print 'P:', P
    E = G.copy()

    # extract multiple MPST until all edges in G are removed
    MPSTs = []
    while len(E.edges()):
        mpst = minimum_spanning_tree(E)
        print '|mpst|:', len(mpst), '|mpst.edges|:', len(mpst.edges())
        MPSTs.append(mpst)
        E.remove_edges_from(mpst.edges())
    print 'Found %s MPST' %(len(MPSTs))
    return MPSTs

def get_sparsified_mpst(MPSTs, K):
    '''
    Get sparsified (uncertain graph with K edges) graph using MPST.

    MPSTs ordered from largest to smallest.

    Expected to have -log(p_e) as weights in G.
    '''

    # sort edges
    sorted_edges = []
    for mpst in MPSTs:
        if len(sorted_edges) + len(mpst.edges()) < K:
            sorted_edges.extend(mpst.edges(data = True))
        else:
            sorted_edges.extend(sorted(mpst.edges(data=True),
                                       key = lambda (u,v,data): exp(1)**(-data["weight"]),
                                       reverse=True))
            break
    return sorted_edges[:K]

def get_sparsified_top(G, K):
    '''
    Get sparsified (uncertain graph with K edges) graph using top most probable edges.

    Expected to have -log(p_e) as weights in G.
    '''
    all_edges = G.edges(data=True)
    sorted_edges = sorted(all_edges,
                          key = lambda (u,v,data): exp(1)**(-data["weight"]),
                          reverse=True)
    return sorted_edges[:K]

def get_sparsified_random(G, K):
    '''
    Get sparsified (uncertain graph with K edges) graph using random edges.

    Expected to have -log(p_e) as weights in G.
    '''
    all_edges = G.edges(data=True)
    random.shuffle(all_edges)
    random_edges = random.sample(all_edges, K)
    return random_edges[:K]

def get_sparsified_MPST_remove_leaves(G, K, directed=False):
    '''
    Sparsify graph using most probable spanning tree.
    If |MPST| < K, then add most probable edges that are not included.
    If |MPST| > K, then remove edges that are adjacent to leaves
    '''
    G_edges = G.edges(data=True)
    if directed:
        # MPST_edges = branchings.minimum_spanning_arborescence(G, attr='weight').edges(data=True)
        pass
    else:
        MPST_edges = list(nx.minimum_spanning_edges(G,weight='weight',data=True))
    edges = [e for e in G_edges if e not in MPST_edges]
    mp_edges = sorted(edges,
                    key = lambda (u,v,d): exp(1)**(-d["weight"]),
                    reverse = True)
    if len(MPST_edges) <= K:
        MPST_edges.extend(mp_edges[:(K - len(MPST_edges))])
    else:
        # remove edges that are adjacent to leaves (keeping connectivity)
        # if ties remove with lowest probability (keeping probability)
        #TODO check why in case of directed MPST it doesn't work
        MPST = nx.Graph(MPST_edges)
        degrees = dict()
        leaves = set()
        for u in MPST:
            degrees[u] = len(MPST[u])
            if degrees[u] == 1:
                v, d = MPST[u].items()[0]
                leaves.add((u,v,d["weight"]))
        for _ in range(len(MPST_edges) - K):
            u,v,d = min(leaves, key = lambda (u,v,d): exp(1)**(-d))
            MPST.remove_edge(u,v)
            leaves.remove((u,v,d))
            v_edges = MPST[v].items()
            if len(v_edges) == 1:
                w, t = v_edges[0]
                leaves.add((v,w,t["weight"]))
            elif len(v_edges) == 0:
                leaves.remove((v,u,d))
        print len(MPST.edges()), K
        MPST_edges = MPST.edges(data=True)
    return MPST_edges

def get_sparsified_MPST_keep_connected(G, K):
    '''
    Sparsify graph using most probable spanning tree.
    If K < |MPST|, then add most probable edges that are not included.
    If K > |MPST|, then remove edges based on the number of neighbors
    for endpoints of the edges.

    :param G: undirected graph with -log(p_e) on edges
    :param K: number of edges to preserve
    :return: edges with probabilities -log(p_e) of size K
    '''
    G_edges = G.edges(data=True)
    print 'Finding MPST'
    MPST_edges = list(nx.minimum_spanning_edges(G,weight='weight',data=True))
    print 'Found spanning tree'
    if len(MPST_edges) <= K:
        print 'Start sorting remaining edges'
        edges = [e for e in G_edges if e not in MPST_edges]
        mp_edges = sorted(edges,
                    key = lambda (u,v,d): exp(1)**(-d["weight"]),
                    reverse = True)
        print 'Finished sorting edges'
        MPST_edges.extend(mp_edges[:(K - len(MPST_edges))])
    else:
        #remove edges that will not make isolated vertices
        print 'Start removing spare edges'
        MPST = nx.Graph()
        MPST.add_edges_from(MPST_edges)
        # edges where both endpoints have at least 2 edges (remove these edges first)
        power1 = {tuple(sorted((e[0],e[1]))): len(MPST[e[0]]) + len(MPST[e[1]]) - 2 for e in MPST.edges() if len(MPST[e[0]]) != 1 and len(MPST[e[1]]) != 1}
        # edges where at least one vertex has only one edge (remove these edges second)
        power2 = {tuple(sorted((e[0],e[1]))): len(MPST[e[0]]) + len(MPST[e[1]]) - 2 for e in MPST.edges() if len(MPST[e[0]]) == 1 or len(MPST[e[1]]) == 1}

        def _decrease_power(MPST, endpoint, power1, power2):
            for node in MPST[endpoint]:
                edge = tuple(sorted((endpoint, node)))
                if edge in power1:
                    power1[edge] -= 1
                elif edge in power2:
                    power2[edge] -= 1
                else:
                    raise ValueError, "Edge %s should be in MPST" %edge

        while len(MPST.edges()) > K:
            if len(power1) > 0:
                e, val = max(power1.iteritems(), key = lambda(dk, dv): dv)
                # print e, val
                power1.pop(e)
                MPST.remove_edge(e[0],e[1])
                _decrease_power(MPST, e[0], power1, power2)
                _decrease_power(MPST, e[1], power1, power2)
            elif len(power2) > 0:
                e, val = max(power2.iteritems(), key = lambda(dk, dv): dv)
                # print e, val
                power2.pop(e)
                MPST.remove_edge(e[0],e[1])
                _decrease_power(MPST, e[0], power1, power2)
                _decrease_power(MPST, e[1], power1, power2)
            else:
                raise ValueError, "No more edges"

        MPST_edges = MPST.edges(data=True)
    return MPST_edges

def exp_degrees(G):
    '''
    Expected to have -log(p_e) as weights in G.
    :param G:
    :return:
    '''
    in_d = dict(zip(G, [0]*len(G)))
    out_d = dict(zip(G, [0]*len(G)))
    for v in G:
        for u in G[v]:
            out_d[v] += exp(1)**(-G[v][u]['weight'])
            in_d[u] += exp(1)**(-G[v][u]['weight'])
    return in_d, out_d

def pres_dir_degrees_ChungLu(edges, degrees, dir="in"):
    '''
    Preserve directed degrees
    :param edges: directed edges
    :param in_degrees: dictionary of expected in_degrees
    :return: Q -- graph from edges with preserved in_degrees

    Note the problem that sometimes probabilities = 0 or > 1.
    '''
    Q = nx.DiGraph()
    Q.add_edges_from(edges)
    if dir == "in":
        denom = {v: sum([degrees[u] for (u,_) in Q.in_edges(v)]) for v in Q}
    elif dir == "out":
        denom = {v: sum([degrees[u] for (_,u) in Q.out_edges(v)]) for v in Q}

    assert denom > 0
    for e in edges:
        if dir == "in":
            sum_d = denom[e[1]]
        elif dir == "out":
            sum_d = denom[e[0]]
        if degrees[e[0]]*degrees[e[1]] > sum_d:
            print '(',e[0],e[1],') -->', degrees[e[0]]*degrees[e[1]], sum_d
        print '(',e[0],e[1],') -->', degrees[e[0]], degrees[e[1]], sum_d
        p = degrees[e[0]]*degrees[e[1]]/sum_d
        print p
        Q.add_edge(*e, **{"weight": -log(p)})
    return Q

def pres_dir_degrees_equally(edges, degrees, dir="in"):
    '''
    Preserve directed degrees
    :param edges: directed edges
    :param degrees: dictionary of expected in_degrees
    :return: Q -- graph from edges with preserved in_degrees
    '''
    Q = nx.DiGraph()
    Q.add_edges_from(edges)
    for e in edges:
        if dir == "in":
            p = min(degrees[e[1]]/len(Q.in_edges(e[1])), 1)
        elif dir == "out":
            p = min(degrees[e[0]]/len(Q.out_edges(e[0])), 1)
        Q.add_edge(e[0],e[1],weight=-log(p))
    return Q

def get_undirected_prob(Q):
    '''

    :param Q: directed graph
    :return: undirected graph where probabiltiies is some function of previous probabilities
    '''
    G = nx.Graph()
    for e in Q.edges(data=True):
        u,v,log_p1 = e[0],e[1],e[2]['weight']
        log_p2 = Q[v][u]['weight']
        p = exp(1)**(-log_p1) + exp(1)**(-log_p2)
        G.add_edge(u,v,{'weight': -log(p/2)})
    return G

def get_graph_from_file(filename, directed=False):
    SP = nx.Graph()
    if directed:
        SP = nx.DiGraph()
    with open(filename) as f:
        for line in f:
            u, v, p = map(float, line.split())
            if u != v:
                SP.add_edge(int(u), int(v), weight = -log(p))
    return SP

def save_for_LP(f1, f2, G, G_orig, f3 = "edge_order.txt", f4 = "D.txt"):
    '''
    Saves matrix A and vector b to files f1 and f2,
    for linear programming min||Ax - b||, s.t. 0 <= x <= 1 in MATLAB.

    Save edge order to f3,
    sum of probabilities of sparsified edges to f4.

    Expected to have -log(p_e) as weights in G.
    '''

    edge_order = dict()
    for i, e in enumerate(G.edges()):
        edge_order[(e[0],e[1])] = i + 1
    node_order = dict()
    for i, u in enumerate(G):
        node_order[u] = i + 1
    wi = dict()
    for u in G:
        wi[u] = sum([exp(1)**(-G[u][v]["weight"]) for v in G[u]])

    with open(f1, "w+") as f:
        for e in G.edges():
            e1, e2 = e[0], e[1]
            f.write("%s %s %s\n" %(node_order[e1], edge_order[(e1, e2)], 1))
            f.write("%s %s %s\n" %(node_order[e2], edge_order[(e1, e2)], 1))
    with open(f2, "w+") as f:
        for u in G:
            f.write("%s %s %s\n" %(node_order[u], 1, wi[u]))
    with open(f3, "w+") as f:
        for e in G.edges():
            f.write("%s %s %s\n" %(e[0], e[1], edge_order[(e[0],e[1])]))

    nodes = set(G_orig).difference(G)
    print "# nodes: ", len(nodes),
    d = dict() # degree of remaining nodes
    for u in nodes:
        d[u] = 0
        for v in G_orig[u]:
            p = exp(1)**(-G_orig[u][v]["weight"])
            d[u] += p
    surplus = sum(d.values())
    print "surplus", surplus
    with open(f4, "w+") as f:
        f.write('%s\n' %len(G_orig))
        f.write('%s\n' %surplus)

def save_for_LP_dir (Ainf, Aoutf, dinf, doutf, Q, G, e_ordf, D):
    '''
    Saves matrices for linear program max|x| s.t. din > Ain*x, dout > Aout*x
    :param Ainf: filenames
    :param Aoutf:
    :param dinf:
    :param doutf:
    :param Q: sparsified directed graph
    e_ord: edges for which probabilities should be found
    D: number of nodes and surplus for sparsified nodes
    G: original graph
    :return:
    '''

    assert type(Q) == type(nx.DiGraph())
    edges = Q.edges()

    # calculate in_degree, out_degree of u
    win = dict()
    wout = dict()
    for u in G:
        win[u] = 0
        wout[u] = 0
        for (v,_) in G.in_edges(u):
            p = exp(1)**(-G[v][u]['weight'])
            win[u] += p
        for (_,v) in G.out_edges(u):
            p = exp(1)**(-G[u][v]['weight'])
            wout[u] += p

    # save order of edges
    e_ord = dict() # order of edges
    v_ord = dict() # order of nodes
    count = 1
    for e in edges:
        e_ord[(e[0],e[1])] = count
        e_ord[(e[1],e[0])] = count
        count += 1
    count = 1
    for u in Q:
        v_ord[u] = count
        count += 1
    with open(e_ordf, "w+") as f:
        for e in edges:
            f.write("%s %s\n" %(e[0], e[1]))

    # save Ain and Aout
    nodes = set()
    with open(Ainf, "w+") as f1, open(Aoutf, "w+") as f2:
        for e in edges:
            f1.write("%s %s %s\n" %(v_ord[e[1]], e_ord[e], 1))
            f2.write("%s %s %s\n" %(v_ord[e[0]], e_ord[e], 1))
            nodes.add(e[0])
            nodes.add(e[1])

    sp_n = len(nodes)
    sp_m = len(edges)
    print len(nodes), len(edges)

    # save din and dout
    with open(dinf, "w+") as f1, open(doutf, "w+") as f2:
        for i, u in enumerate(nodes):
            if i == 54307:
                print u, v_ord[u], win[u], wout[u]
            f1.write("%s %s %s\n" %(v_ord[u], 1, win[u]))
            f2.write("%s %s %s\n" %(v_ord[u], 1, wout[u]))

    # save number of nodes and surplus
    n = len(G)
    nodes = set(G).difference(Q)
    print "# nodes: ", len(nodes),
    in_surplus = sum([win[u] for u in nodes])
    out_surplus = sum([wout[u] for u in nodes])
    print "in-degree surplus", in_surplus
    print "out-degree surplus", out_surplus
    with open(D, "w+") as f:
        f.write("%s\n" %(n))
        f.write("%s\n" %(sp_n))
        f.write("%s\n" %(sp_m))
        f.write("%s\n" %(in_surplus + out_surplus))

def construct_lp_graph(lp_file, e_file, output_file):
    with open(lp_file) as f, open(e_file) as g, open(output_file, "w+") as h:
        for x, y in izip(f,g):
            d1 = x.split()
            d2 = y.split()
            h.write("%s %s %s\n" %(d2[0], d2[1], d1[0]))

def get_possible_worlds(graph_file, PW_folder, I):
    '''
    Creates possible worlds from the graph file
    :param graph_file: file of edges with probabilities
    :param PW_folder: folder where to write PWs
    :param I: number of possible worlds
    :return:
    '''
    with open(graph_file) as f:
        f_lines = f.readlines()
    for i in range(I):
        with open(PW_folder + "PW%s.txt" %(i+1), "w+") as f:
            for line in f_lines:
                d = map(float, line.split())
                if random.random() < d[2]:
                    f.write(line)

def estimate_memory(f, args, kw={}, interval=.5):
    return max(memory_usage((f, args, kw), interval=interval))

def main(dt='../Datasets/Ep_hep_random1.txt', pt='percentages.txt', st='Result/Q'):
    '''
    Run Sparsification method to obtain
    :return:
    '''
    arguments = sys.argv
    if len(arguments) == 1:
        dataset_filename = dt
        percentages_filename = pt
        sparsified_filename = st
    elif len(arguments) == 4:
        dataset_filename = arguments[1] # the path to dataset filename
        percentages_filename = arguments[2] # filename with specified percentages of the total number of edges
        sparsified_filename = arguments[3] # the path to the output place of sparsified graphs
    else:
        raise ValueError, 'You should provide 3 arguments: dataset, percentages, name of sparsified file'

    percentages = []
    with open(percentages_filename) as f:
        for line in f:
            percentages.append(int(line))

    for perc in percentages:
        print 'Percentage: %s%%' %(perc)
        print '*Reading graph...'
        G = get_graph_from_file(dataset_filename)
        K = int(len(G.edges())*perc/100)
        print '*Obtaining sparsified edges...'
        kept_edges = get_sparsified_MPST_keep_connected(G, K)
        Q = nx.Graph(kept_edges)

        print '*Saving sparsified edges to files for MATLAB processing...'
        directory = 'tmp/'
        if not os.path.exists(directory):
            os.makedirs(directory)
        copyfile(percentages_filename, directory + os.path.basename(percentages_filename))
        save_for_LP(directory+'A' + str(perc) + '.dat', directory+'b' + str(perc) + '.dat', Q, G, directory+'e' + str(perc) + '.dat', directory+'D' + str(perc) + '.dat')

    print '*Assigning probabilities'
    call(['matlab', '-nodisplay', '-r', "run('LP.m')"])

    print '*Constructing sparsified graph...'
    for perc in percentages:
        construct_lp_graph(directory+'x' + str(perc) + '.dat', directory+'e' + str(perc) + '.dat', sparsified_filename + str(perc) + '.txt')

    rmtree(directory)  # remove tmp folder

if __name__ == "__main__":
    time2execute = time.time()

    # Run from command-line
    # example: python MPST.py ../Datasets/Ep_hep_random1.txt percentages.txt Q
    # ../Datasets/Ep_hep_random1.txt -- filename of dataset
    # percentages.txt -- filename with percentages
    # Q -- name of sparsified graphs (they will be written in the form such as Q10.txt, where 10 is percent)
    main()

    print "Finished execution in %s sec" %(time.time() - time2execute)
    console = []