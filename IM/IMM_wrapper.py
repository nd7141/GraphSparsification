'''
http://delivery.acm.org/10.1145/2730000/2723734/p1539-tang.pdf?ip=89.106.174.122&id=2723734&acc=OPEN&key=4D4702B0C3E38B35%2E4D4702B0C3E38B35%2E4D4702B0C3E38B35%2E6D218144511F3437&CFID=667647771&CFTOKEN=97444750&__acm__=1434964493_0e76cda204d26fbea8a507089b868b5f
'''
from __future__ import division
__author__ = 'sivanov'
import os, sys, time, json

def getNM(filename):
    nodes = set()
    edges = 0
    with open(filename) as f:
        for line in f:
            edges += 1
            nodes.update(map(int, line.split()[:2]))
    return len(nodes), edges

if __name__ == "__main__":
    start2exec = time.time()

    if len(sys.argv) not in [4]:
        assert ValueError, 'command: python IMM_wrapper.py dataset_folder eps k)'

    # dataset_folder should contain graph_ic.inf file with a graph u v p

    dataset = sys.argv[1] # path to dataset
    eps = sys.argv[2] # epsilon

    IMM_time = []
    Spread_time = []
    begin_seed = 45
    end_seed = 55
    step_seed = 5

    # create directory if not exist
    directory = "IMM_data/"
    if not os.path.exists(directory):
        os.makedirs(directory)

    # remove the output file if exists
    output = directory + 'spread.txt'
    if os.path.exists(output):
        os.remove(output)

    if not os.path.exists(dataset + 'graph_ic.inf'):
        raise ValueError, 'The folder with data should contain file graph_ic.inf with graph'

    # create file attribute.txt
    N, M = getNM(dataset + 'graph_ic.inf')
    with open(dataset+'/attribute.txt', 'w+') as f:
        f.write('n=%s\n' %N)
        f.write('m=%s\n' %M)

    # os.system('make -C ./getPossibleWorlds/')
    os.system('make -C ./Algorithms/IMM/')

    for k in range(begin_seed, end_seed, step_seed):
        if os.path.exists(directory + 'seeds%s.txt' %(k)):
            os.remove(directory + 'seeds%s.txt' %(k))
        os.system('./Algorithms/IMM/imm_discrete -dataset %s -epsilon %s -k %s  -model IC -spread %s -seeds %s' %(dataset, eps, k, output, directory + 'seeds%s.txt' %(k)))

    finish2exec = time.time() - start2exec

    # print '* To get worlds: %s sec' %finish2worlds
    print '* Total execution time: %s sec' %finish2exec

    console = []