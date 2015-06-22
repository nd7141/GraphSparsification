'''

'''
from __future__ import division
__author__ = 'sivanov'
import pygal
import matplotlib.pyplot as plt
import matplotlib
# import seaborn as sns
import numpy as np
import pandas as pd
from scipy import stats, optimize

def visualiseResults(x_lst, y_lst, title="Dataset", filename="tempResults.png",):
    matplotlib.rcParams.update({'font.size': 24})
    fig = plt.figure(figsize=(18, 10))

    # ax.set_xscale("log")
    # ax.set_yscale("log")

    legends = ['Harvester', 'IMM']
    colors = ['b', 'r', 'g', 'm', 'k', u'#abfeaa', u'#cccabc', u'#1111ee', 'y', 'c', u'#fe2fb3']
    marks = ["o", "s", "^", "v", 'x', "<", ">", '8', "<", ">", '8']
    colors = colors[::1]
    marks = marks[::1]
    x_lst.reverse()
    y_lst.reverse()
    legends.reverse()
    colors.reverse()
    marks.reverse()

    plots = []
    # print colors
    for i in range(len(x_lst)):
        plt.plot(x_lst[i], y_lst[i], color=colors[i], linewidth=3)
        p, = plt.plot(x_lst[i], y_lst[i], color = colors[i], marker = marks[i], markersize=10)
        plots.append(p)

    # plt.xlim([9, 200])
    # plt.ylim([770, 820])

    plt.legend(plots, legends, loc=4, prop={'size': 24})
    plt.grid()
    plt.xlabel('Seed set k')
    plt.ylabel('Spread')
    plt.title('%s' %(title), fontsize = 24)
    fig.savefig(filename, dpi=fig.dpi)
    plt.show()

if __name__ == "__main__":

    x_lst = []
    y_lst = []
    with open('Harvester_data/spread.txt') as f:
        x = []
        y = []
        for line in f:
            d = map(float, line.split())
            x.append(d[0])
            y.append(d[1])
        x_lst.append(x)
        y_lst.append(y)
    with open('IMM_data/spread.txt') as f:
        x = []
        y = []
        for line in f:
            d = map(float, line.split())
            x.append(d[0])
            y.append(d[1])
        x_lst.append(x)
        y_lst.append(y)

    visualiseResults(x_lst, y_lst, title='Harvester vs IMM (Hep)', filename="figures/hep.png")

    console = []