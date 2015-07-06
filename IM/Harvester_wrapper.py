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

    if len(sys.argv) not in [4,5]:
        raise ValueError, 'command: python Harvester_wrapper.py dataset R I (N)'

    dataset = sys.argv[1] # path to dataset
    R = sys.argv[2] # number of worlds
    I = sys.argv[3] # number of Monte-Carlo simulations
    if len(sys.argv) == 5:
        N = sys.argv[4] # number of vertices in dataset
    else:
        N, M = getNM(dataset)

    print N, M

    Harvester_time = []
    Spread_time = []
    begin_seed = 10
    end_seed = 155
    step_seed = 5

    # create directory if not exist
    directory = "data/Harvester_data/"
    if not os.path.exists(directory):
        os.makedirs(directory)

    if not os.path.exists(directory + 'PW/'):
        os.makedirs(directory + 'PW/')

    # remove the output file if exists
    output = directory + 'spread.txt'
    if os.path.exists(output):
        os.remove(output)

    # remove the time file if exists
    timing = directory + 'astro_categories_time.txt'
    if os.path.exists(timing):
        os.remove(timing)

    # os.system('make -C ./getPossibleWorlds/')
    os.system('make -C ./Algorithms/Harvester2/')
    os.system('make -C ./getSpread/')

    k_range = range(begin_seed, end_seed, step_seed)
    for k in k_range:
        # run Harveseter
        start2Harvester = time.time()
        os.system('./Algorithms/Harvester2/Harvester %s %s %s %s %s' %(dataset, N, R, str(k), directory + 'seeds' + str(k) + '.txt'))
        finish2Harvester = time.time() - start2Harvester
        Harvester_time.append(finish2Harvester)

        # calculate spread
        start2Spread = time.time()
        os.system('./getSpread/runCascade %s %s %s %s %s' %(dataset, N, directory + 'seeds' + str(k) + '.txt', I, output))
        finish2Spread = time.time() - start2Spread
        Spread_time.append(finish2Spread)

        finish2exec = time.time() - start2exec

    # write timing to file
    with open(timing, 'w+') as f:
        for i in range(len(Harvester_time)):
            f.write("%s %s\n" %(k_range[i], Harvester_time[i]))

    # print '* To get worlds: %s sec' %finish2worlds
    print '* To run Harvester: %s sec' %json.dumps(Harvester_time)
    print '* To get spread: %s sec' %json.dumps(Spread_time)
    print '* Total execution time: %s sec' %finish2exec

    console = []