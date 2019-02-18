"""
A common toolbox containing useful functions for data analysis and plotting
"""
import numpy as np
import matplotlib.pyplot as plt

colors = ['k', 'r', 'g', 'b', 'm', 'y', 'c']*5

def flip_x(ax=plt.gca()):
    """Flip the axis of `Axes` object ax"""
    x1, x2 = ax.get_xlim()
    ax.set_xlim([x2, x1])

def get_range(x, lim1, lim2):
    """Return the index range of x between lim1 and lim2"""
    index1 = np.where(x <= lim2)[0]
    index2 = np.where(x >= lim1)[0]
    index = np.intersect1d(index1, index2)
    return index

def copy_range(x, lim1, lim2):
    """Return a copy of x between lim1 and lim2"""
    x_copy = x[np.where(x >= lim1)[0]]
    x_copy = x_copy[np.where(x_copy <= lim2)[0]]
    return x_copy
