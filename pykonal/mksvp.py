#!/home/sgo/anaconda3/bin/python3
import matplotlib
matplotlib.use('Agg')
import matplotlib.gridspec as gs
import matplotlib.pyplot as plt
import matplotlib.cm as cm
from matplotlib.colors import Normalize
import shutil,math
import pandas as pd
import numpy as np
import pkg_resources
import pykonal
import os,time,sys

from matplotlib import markers
from matplotlib.path import Path
from mpl_toolkits.axes_grid1.inset_locator import zoomed_inset_axes, mark_inset

#Set the velocity gradient in 1/s
sv = pd.read_csv('svpA.csv')
dpX = np.arange(10000, dtype = 'float64')
svX = np.arange(10000, dtype = 'float64')
dvX = np.arange(10000, dtype = 'float64')
print(len(sv.speed))
pd.options.display.float_format = '{:.5f}'.format
spli=20
for i in range(len(sv.speed)-1):
  for t in range(spli):
#    dpX[i*spli+t]=100.*float(i)+100.*float(t)/float(spli)
#    svX[i*spli+t]=(sv.speed[i]*float(spli-t)+sv.speed[i+1]*float(t))/float(spli)
    dpX[i*spli+t]=100.*i+100.*t/spli
    X=(sv.speed[i]*(spli-t)+sv.speed[i+1]*t)/spli
    svX[i*spli+t]=X
    dvX[i*spli+t]=1.0
df = pd.DataFrame({ 'depth' : dpX,'speed' : svX,'degree_of_varience_in_svinv' : dvX}, dtype="float64")
df = df.set_index('depth')
df.to_csv('out.csv')
