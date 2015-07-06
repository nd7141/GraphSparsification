'''

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

    if len(sys.argv) != 4:
        raise ValueError, 'command: python reverse_Harvester_wrapper.py dataset R I'

    dataset = sys.argv[1] # path to dataset
    R = sys.argv[2] # number of worlds
    I = sys.argv[3] # number of Monte-Carlo simulations

    N, M = getNM(dataset)

    print N, M

    reverse_time = []
    begin_T = 15250
    end_T = 15550
    step_T = 50

    # create directory if not exist
    directory = "data/Harvester_data/"
    if not os.path.exists(directory):
        os.makedirs(directory)

    # remove the output file if exists
    output = directory + 'astro_categories_reverse.txt'
    if os.path.exists(output):
        os.remove(output)

    # remove the time file if exists
    timing = directory + 'astro_categories_reverse_time.txt'
    if os.path.exists(timing):
        os.remove(timing)

    os.system('make -C ./Reverse/')

    T_range = range(begin_T, end_T, step_T)
    for T in T_range:
        print 'T:', T
        # run reverse
        start2reverse = time.time()
        os.system('./Reverse/main.o %s %s %s %s %s %s' %(dataset, N, R, T, I, output))
        finish2reverse = time.time() - start2reverse
        reverse_time.append(finish2reverse)

    finish2exec = time.time() - start2exec

    # write timing to file
    with open(timing, 'w+') as f:
        for i in range(len(reverse_time)):
            f.write("%s %s\n" %(T_range[i], reverse_time[i]))

    print '* To run Reverse Harvester: %s sec' %json.dumps(reverse_time)
    print '* Total execution time: %s sec' %finish2exec


    console = []