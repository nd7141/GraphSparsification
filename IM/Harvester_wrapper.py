'''

'''
from __future__ import division
__author__ = 'sivanov'
import os, sys, time

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

    # create directory if not exist
    directory = "data/"
    if not os.path.exists(directory):
        os.makedirs(directory)

    if not os.path.exists(directory + 'PW/'):
        os.makedirs(directory + 'PW/')


    start2worlds = time.time()
    # create possible worlds
    os.system('make -C ./getPossibleWorlds/')
    os.system('./getPossibleWorlds/getPossibleWorlds ' + dataset + " " + str(R) + ' ./'  + directory + 'PW/')
    finish2worlds = time.time() - start2worlds

    start2Harvester = time.time()
    # run Harveseter
    os.system('make -C ./Algorithms/Harvester/')
    os.system('./Algorithms/Harvester/Harvester ' + dataset + ' ' + str(N) + ' ' + directory + 'PW/' + ' ' + str(k) + ' ' + directory + 'seeds.txt')
    finish2Harvester = time.time() - start2Harvester

    start2Spread = time.time()
    # calculate spread
    os.system('make -C ./getSpread/')
    os.system('./getSpread/runCascade ' + dataset + ' ' + str(N) + ' ' + directory + 'seeds.txt ' + str(I) + ' ' + directory + 'spread.txt')
    finish2Spread = time.time() - start2Spread

    finish2exec = time.time() - start2exec

    print '* To get worlds: %s sec' %finish2worlds
    print '* To run Harvester: %s sec' %finish2Harvester
    print '* To get spread: %s sec' %finish2Spread
    print '* Total execution time: %s sec' %finish2exec

    console = []