'''

'''
from __future__ import division
__author__ = 'sivanov'
import os, sys, time, json, subprocess

def getNM(filename):
    nodes = set()
    edges = 0
    with open(filename) as f:
        for line in f:
            edges += 1
            nodes.update(map(int, line.split()[:2]))
    return len(nodes), edges

def get_spread(spreads, folder, eps, k, seeds_filename, spread_filename):
    if os.path.exists(seeds_filename):
        os.remove(seeds_filename)
    os.system('./Algorithms/IMM/imm_discrete -dataset %s -epsilon %s -k %s  -model IC -seeds %s' %(folder, eps, k, seeds_filename))
    os.system('./getSpread/runCascade %s %s %s %s %s' %(dataset, N, directory + 'seeds.txt', I, spread_filename))
    spread = map(float, subprocess.check_output(['tail', '-1', spread_filename]).split())[1]
    spreads[k] = spread

if __name__ == "__main__":
    start2exec = time.time()

    if len(sys.argv) != 5:
        raise ValueError, 'command: python reverse_IMM_wrapper.py folder eps I dataset'

    folder = sys.argv[1] # path to the folder with dataset
    eps = sys.argv[2] # epsilon
    I = sys.argv[3] # number of MC simulations
    dataset = sys.argv[4] # unidirected version of dataset for spread calculations

    N, M = getNM(dataset)

    print N, M

    reverse_time = []
    begin_T = 100
    end_T = 1050
    step_T = 50

    # create directory if not exist
    directory = "data/IMM_data/"
    if not os.path.exists(directory):
        os.makedirs(directory)

    seeds_filename = directory+'seeds.txt'
    spread_filename = directory + "spread.txt"

    # remove the output file if exists
    output = directory + 'gnu_categories_reverse.txt'
    if os.path.exists(output):
        os.remove(output)

    os.system('make -C ./Algorithms/IMM/')

    for T in range(begin_T, end_T, step_T):
        # run reverse
        spreads = dict()
        spreads[0] = 0
        High = 1
        get_spread(spreads, folder, eps, High, seeds_filename, spread_filename)
        spread = spreads[High]
        while spread < T:
            High *= 2
            get_spread(spreads, folder, eps, High, seeds_filename, spread_filename)
            spread = spreads[High]
        Low = High//2
        while Low+1 < High:
            k = Low + (High - Low)//2
            get_spread(spreads, folder, eps, k, seeds_filename, spread_filename)
            spread = spreads[k]
            if spread > T:
                High = k
            else:
                Low = k
        with open(output, 'a+') as f:
            if spreads[Low] > T:
                f.write("%s %s\n" %(T, Low))
            else:
                f.write("%s %s\n" %(T, High))

    finish2exec = time.time() - start2exec

    console = []