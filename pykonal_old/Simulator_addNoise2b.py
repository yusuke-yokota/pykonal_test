#!/home/sgo/anaconda3/bin/python3
import matplotlib
matplotlib.use('Agg')
import matplotlib.gridspec as gs
import matplotlib.pyplot as plt
import matplotlib.cm as cm
from matplotlib.colors import Normalize
import shutil,math
import configparser
import pandas as pd
import numpy as np
import pykonal
import os,time,sys
from scipy.interpolate import Rbf
from scipy.special import eval_chebyt
import colorednoise as cn
##
sv = pd.read_csv('svpC.csv')
icfg = configparser.ConfigParser()
icfg.read('Settings.ini', 'UTF-8')
ld = pd.read_csv('linesample.csv')
kf='~/sgobs/garpos/simulator/obsdata/IMAG/IMAG.1911.kaiyo_k4-obs.csv'
dkf = pd.read_csv(kf,comment='#')
sf='~/sgobs/garpos/simulator/initcfg/IMAG/IMAG.1911.kaiyo_k4-noise.ini'
skf = pd.read_csv(sf,comment='#')
##
#################################
##Input parameters start #####
# GNSS noise ######
sample= icfg.get("HyperParameters", "Sample")
ghptb = icfg.get("HyperParameters", "Ghptb")
gvptb = icfg.get("HyperParameters", "Gvptb")
knot  = icfg.get("HyperParameters", "Knot")
lknot  = float(knot)*1.852
lhptb = icfg.get("HyperParameters", "Leptb")
lnptb = icfg.get("HyperParameters", "Lnptb")
luptb = icfg.get("HyperParameters", "Luptb")
sample,ghptb,gvptb,lhptb,lnptb,luptb = int(sample),float(ghptb),float(gvptb),float(lhptb),float(lnptb),float(luptb)
# node setting ####
nlon  = icfg.get("HyperParameters", "NodeLon")
ndep  = icfg.get("HyperParameters", "NodeDep")
nlat  = icfg.get("HyperParameters", "NodeLat")
dlon  = icfg.get("HyperParameters", "DnodLon")
ddep  = icfg.get("HyperParameters", "DnodDep")
dlat  = icfg.get("HyperParameters", "DnodLat")
slon  = icfg.get("HyperParameters", "SrceLon")
sdep  = icfg.get("HyperParameters", "SrceDep")
slat  = icfg.get("HyperParameters", "SrceLat")
asiz  = icfg.get("HyperParameters", "ArraSiz")
nlon,ndep,nlat,dlon,ddep,dlat,slon,sdep,slat,asiz = int(nlon),int(ndep),int(nlat),float(dlon),float(ddep),float(dlat),int(slon),int(sdep),int(slat),int(asiz)
# gradient setting ####
m_s_dp = icfg.get("HyperParameters", "Nod_cut").split()
m_s_km = icfg.get("HyperParameters", "Vel_cut").split()
o_chet = icfg.get("HyperParameters", "O_Chebyt")
c_chet = icfg.get("HyperParameters", "Chebyt").split()
d_chet = icfg.get("HyperParameters", "D_Chebyt")
m_s_dp = np.array(list(map(int, m_s_dp)))
m_s_km = np.array(list(map(float, m_s_km)))
c_chet = np.array(list(map(float, c_chet)))
o_chet,d_chet = int(o_chet),float(d_chet)
#################################
##Input parameters  end  #####
#################################
##Initialize the solver. #####
##############################
beta10 = 1.0
beta15 = 1.5
beta20 = 2.0
#################################
ltime  = int(dkf.ST.iloc[-1])+1 # int(3600.*kyori/lknot)
ltime  = 86400
for epi in range(sample):
  gpserr = cn.powerlaw_psd_gaussian(beta15, 3*ltime)
  gpe = (gpserr[       :  ltime]) # - np.average(gpserr[       :  ltime]))
  gpn = (gpserr[  ltime:2*ltime]) # - np.average(gpserr[  ltime:2*ltime]))
  gpu = (gpserr[2*ltime:3*ltime]) # - np.average(gpserr[2*ltime:3*ltime]))
  df= pd.DataFrame({'e':gpe,'n':gpn,'u':gpu})
  df = df.round({'e': 10, 'n': 10, 'u': 10})
  of='%s.noise' % str(epi).zfill(5)
  hd = "# noisefile\n"
  outf = open(of,"w")
  outf.write(hd)
  outf.close()
  df.to_csv(of,mode='a')
