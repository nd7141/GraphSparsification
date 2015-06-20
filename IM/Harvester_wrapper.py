'''

'''
from __future__ import division
__author__ = 'sivanov'
import os, sys, time, json

def getN(filename):
    nodes = set()
    with open(filename) as f:
        for line in f:
            nodes.update(map(int, line.split()[:2]))
    return len(nodes)

if __name__ == "__main__":
    start2exec = time.time()

    if len(sys.argv) != 5:
        assert ValueError, 'command: python Harvester_wrapper.py dataset k R I'

    dataset = sys.argv[1] # path to dataset
    k = sys.argv[2] # number of seeds
    R = sys.argv[3] # number of worlds
    I = sys.argv[4] # number of Monte-Carlo simulations

    N = getN(dataset)
    Harvester_time = []
    Spread_time = []
    begin_seed = 10
    end_seed = 55
    step_seed = 5

    # create directory if not exist
    directory = "data/"
    if not os.path.exists(directory):
        os.makedirs(directory)

    if not os.path.exists(directory + 'PW/'):
        os.makedirs(directory + 'PW/')

    # remove the output file if exists
    output = directory + 'spread.txt'
    if os.path.exists(output):
        os.remove(output)

    os.system('make -C ./getPossibleWorlds/')
    os.system('make -C ./Algorithms/Harvester/')
    os.system('make -C ./getSpread/')

    # create possible worlds
    start2worlds = time.time()
    os.system('./getPossibleWorlds/getPossibleWorlds ' + dataset + " " + str(R) + ' ./'  + directory + 'PW/')
    finish2worlds = time.time() - start2worlds

    for k in range(begin_seed, end_seed, step_seed):
        # run Harveseter
        start2Harvester = time.time()
        os.system('./Algorithms/Harvester/Harvester %s %s %s %s %s' %(dataset, N, directory+'PW/', str(k), directory + 'seeds' + str(k) + '.txt'))
        # os.system('./Algorithms/Harvester/Harvester ' + dataset + ' ' + str(N) + ' ' + directory + 'PW/' + ' ' + str(k) + ' ' + directory + 'seeds.txt')
        finish2Harvester = time.time() - start2Harvester
        Harvester_time.append(finish2Harvester)

        # calculate spread
        start2Spread = time.time()
        os.system('./getSpread/runCascade %s %s %s %s %s' %(dataset, N, directory + 'seeds' + str(k) + '.txt', I, output))
        # os.system('./getSpread/runCascade ' + dataset + ' ' + str(N) + ' ' + directory + 'seeds.txt ' + str(I) + ' ' + directory + 'spread' + str(k) + '.txt')
        finish2Spread = time.time() - start2Spread
        Spread_time.append(finish2Spread)

        finish2exec = time.time() - start2exec

    print '* To get worlds: %s sec' %finish2worlds
    print '* To run Harvester: %s sec' %json.dumps(Harvester_time)
    print '* To get spread: %s sec' %json.dumps(Spread_time)
    print '* Total execution time: %s sec' %finish2exec

    console = []