"""General plotting settings for matplotlib"""

import matplotlib
import matplotlib.pyplot as plt

RCPARAM = {
    'lines.linewidth'   :   2,
    'lines.markersize'  :   10,
    'font.size'         :   30,
    'xtick.labelsize'   :   28,
    'ytick.labelsize'   :   28,
    }
plt.style.use(RCPARAM)
matplotlib.rc('text', usetex=True)
