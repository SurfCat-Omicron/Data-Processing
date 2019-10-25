"""General plotting settings for matplotlib"""

import matplotlib
import matplotlib.pyplot as plt

RCPARAM = {
    # Default text/line sizes
    'lines.linewidth'   :   2,
    'lines.markersize'  :   10,
    'font.size'         :   18,

    'xtick.labelsize'   :   24,
    'ytick.labelsize'   :   24,
    'axes.labelsize'    :   28,
    'axes.titlesize'    :   32,

    # Default figure size (important for 'plt.savefig'
    'figure.figsize': (6.93*2, 4.15*2),
    'figure.constrained_layout.use': True,

    # Save figure options
    'savefig.dpi': 250,
    }
plt.style.use(RCPARAM)
matplotlib.rc('text', usetex=True)

