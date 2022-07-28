import math
import getdist
import getdist.plots as gplot
import matplotlib as mpl
import matplotlib.pyplot as plt
#import os
#import re
import numpy as np
mpl.use('Agg')
mpl.rc('font', **{'family': 'serif', 'serif': ['Computer Modern']})
#mpl.rc('text', usetex=True)

MTHdir = '/project/d/dcurtin/jpbarron/ADM/montepython/Lensing_BAO_Planck/ADM/chains/ADM_FixedParams_mp10GeV_me5MeV_alpha0.01_/'

mcdir = [MTHdir]

samples = getdist.loadMCSamples('/project/d/dcurtin/jpbarron/ADM/montepython/Lensing_BAO_Planck/ADM/chains/ADM_FixedParams_mp10GeV_me5MeV_alpha0.01_/2022-05-08_400000_', settings={'ignore_rows': 0.2})
g = gplot.get_subplot_plotter()


g.triangle_plot(samples,params=['r_all_twin','Delta_N_twin'])

g.export('/project/d/dcurtin/jpbarron/ADM/montepython/Lensing_BAO_Planck/ADM/chains/ADM_FixedParams_mp10GeV_me5MeV_alpha0.01_/plots/r_vs_deltaNeff.png')