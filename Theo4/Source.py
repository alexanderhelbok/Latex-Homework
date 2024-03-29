import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from scipy.stats import norm
from scipy.optimize import curve_fit


def set_square_lim(ax):
    # get x and y limits
    x_lim = ax.get_xlim()
    y_lim = ax.get_ylim()
    # find the biggest and set all to this value to create a square plot
    max_lim = max(abs(x_lim[0]), abs(y_lim[0]), abs(x_lim[1]), abs(y_lim[1]))
    ax.set_xlim(-max_lim, max_lim)
    ax.set_ylim(-max_lim, max_lim)

    return max_lim

# define plot parameters
plt.rcParams['xtick.direction'] = 'in'
plt.rcParams['ytick.direction'] = 'in'
# enable minor ticks
plt.rcParams['xtick.minor.visible'] = True
plt.rcParams['ytick.minor.visible'] = True
# enable ticks on top and right
plt.rcParams['xtick.top'] = True
plt.rcParams['ytick.right'] = True
# increase tick length
plt.rcParams['xtick.major.size'] = 7
plt.rcParams['xtick.minor.size'] = 3.5
plt.rcParams['ytick.major.size'] = 7
plt.rcParams['ytick.minor.size'] = 3.5
# increase tick label pad
plt.rcParams['xtick.major.pad'] = 7
plt.rcParams['xtick.minor.pad'] = 7
plt.rcParams['ytick.major.pad'] = 7
plt.rcParams['ytick.minor.pad'] = 7
# increase border width
plt.rcParams['axes.linewidth'] = 1.25
# increase legend axespad
plt.rcParams['legend.borderaxespad'] = 1
# enable latex font in math mode
# plt.rc('text', usetex=True)  # enable use of LaTeX in matplotlib
# plt.rc('font', family="sans-serif", serif="cm", size=14)  # font settings
# plt.rc('text.latex', preamble=r'\usepackage{sansmath} \usepackage[' + "cm" + r']{sfmath} \sansmath \sffamily')