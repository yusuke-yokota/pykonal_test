#!/home/sgo/anaconda3/bin/python3
import matplotlib
matplotlib.use('Agg')
import matplotlib.gridspec as gs
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import shutil
import pandas as pd
import numpy as np
import pykonal
import time
import sys
#####
r=0.5 #v2/v1
#s=0.2
#t=0.5
s=21
st=9
t=st/s
s=s/30.
st=st/30.
#####
b=2

a=30
v1=a/np.sqrt(1+r**2)
#v1=a
v2=v1*r

x1=(s*(1+t)*v1-(2-s-s*t)*b*v2)/(s*s*t)
x2=(s*t*v1-(2-s*t)*b*v2)/(s*s*(t-1))
print("x1: {:4.1f}".format(x1),"x2: {:4.1f}".format(x2))
print(st*30.,s*30.)
print(v1,v2)
