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
svX= [[] for i in range(10000)]
print(len(sv.speed))
#for i in range(len(sv.speed)-1):
#  for t in range(10):
#    svX[i*10+t]=(sv.speed[i]*(10-t)+sv.speed[i+1]*t)*0.1
