"""
A common toolbox containing useful functions for data analysis and plotting
"""

def flip_x(ax):
    x1, x2 = ax.get_xlim()
    ax.set_xlim([x2, x1])
