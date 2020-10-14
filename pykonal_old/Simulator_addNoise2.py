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
slon1,sdep1,slat1=slon,     sdep,slat+asiz
slon2,sdep2,slat2=slon+asiz,sdep,slat
slon3,sdep3,slat3=slon,     sdep,slat-asiz
slon4,sdep4,slat4=slon-asiz,sdep,slat
##############################
beta10 = 1.0
beta15 = 1.5
beta20 = 2.0
#################################
ltime  = int(dkf.ST.iloc[-1])+1 # int(3600.*kyori/lknot)
dkf.ant_e0 = dkf.ant_e0/1.e3
dkf.ant_n0 = dkf.ant_n0/1.e3
dkf.ant_u0 = dkf.ant_u0/1.e3
for epi in range(sample):
  gpserr = cn.powerlaw_psd_gaussian(beta15, 3*ltime)
  gpe = ghptb * (gpserr[       :  ltime]) # - np.average(gpserr[       :  ltime]))
  gpn = ghptb * (gpserr[  ltime:2*ltime]) # - np.average(gpserr[  ltime:2*ltime]))
  gpu = gvptb * (gpserr[2*ltime:3*ltime]) # - np.average(gpserr[2*ltime:3*ltime]))
#################################
  for i in range(len(dkf)):
    dkf.ant_e1[i] = dkf.ant_e0[i]+gpe[int(dkf.ST[i])]
    dkf.ant_n1[i] = dkf.ant_n0[i]+gpn[int(dkf.ST[i])]
    dkf.ant_u1[i] = dkf.ant_u0[i]+gpu[int(dkf.ST[i])]
#    nin = int(df.ant_n0[i]/dlat)
#    nie = int(df.ant_e0[i]/dlon)
    if dkf.MT[i] == 'M01':
      dkf.ResiTT[i] = dkf.TT[i] + 2.*(np.sqrt((dkf.ant_e1[i]-skf.s01[i])**2 + (dkf.ant_n1[i]-skf.s02[i])**2 + (dkf.ant_u1[i]-skf.s03[i])**2) - np.sqrt((dkf.ant_e0[i]-skf.s01[i])**2 + (dkf.ant_n0[i]-skf.s02[i])**2 + (dkf.ant_u0[i]-skf.s03[i])**2)) / skf.dd[i]
    elif dkf.MT[i] == 'M02':
      dkf.ResiTT[i] = dkf.TT[i] + 2.*(np.sqrt((dkf.ant_e1[i]-skf.s11[i])**2 + (dkf.ant_n1[i]-skf.s12[i])**2 + (dkf.ant_u1[i]-skf.s13[i])**2) - np.sqrt((dkf.ant_e0[i]-skf.s11[i])**2 + (dkf.ant_n0[i]-skf.s12[i])**2 + (dkf.ant_u0[i]-skf.s13[i])**2)) / skf.dd[i]
    elif dkf.MT[i] == 'M03':
      dkf.ResiTT[i] = dkf.TT[i] + 2.*(np.sqrt((dkf.ant_e1[i]-skf.s21[i])**2 + (dkf.ant_n1[i]-skf.s22[i])**2 + (dkf.ant_u1[i]-skf.s23[i])**2) - np.sqrt((dkf.ant_e0[i]-skf.s21[i])**2 + (dkf.ant_n0[i]-skf.s22[i])**2 + (dkf.ant_u0[i]-skf.s23[i])**2)) / skf.dd[i]
    elif dkf.MT[i] == 'M04':
      dkf.ResiTT[i] = dkf.TT[i] + 2.*(np.sqrt((dkf.ant_e1[i]-skf.s31[i])**2 + (dkf.ant_n1[i]-skf.s32[i])**2 + (dkf.ant_u1[i]-skf.s33[i])**2) - np.sqrt((dkf.ant_e0[i]-skf.s31[i])**2 + (dkf.ant_n0[i]-skf.s32[i])**2 + (dkf.ant_u0[i]-skf.s33[i])**2)) / skf.dd[i]
#################################
#################################
  df = pd.DataFrame({'SET': dkf.SET, 'LN': dkf.LN,'MT': dkf.MT,'TT': dkf.ResiTT,'ResiTT': np.nan,'TakeOff': np.nan,'ST': dkf.ST,'ant_e0': dkf.ant_e0*1.e3,'ant_n0': dkf.ant_n0*1.e3,'ant_u0': dkf.ant_u0*1.e3,'head0': dkf.head0,'pitch0': dkf.pitch0,'roll0': dkf.roll0,'dSV0': dkf.dSV0,'RT': dkf.RT,'ant_e1': dkf.ant_e1*1.e3,'ant_n1': dkf.ant_n1*1.e3,'ant_u1': dkf.ant_u1*1.e3,'head1': dkf.head1,'pitch1': dkf.pitch1,'roll1': dkf.roll1,'dSV1': dkf.dSV1,'flag': dkf.flag})
  df = df.round({'TT': 6, 'ST': 6, 'ant_e0': 5, 'ant_n0': 5, 'ant_u0': 5, 'RT': 6, 'ant_e1': 5, 'ant_n1': 5, 'ant_u1': 5})
  ix = 2001 + epi % 12 + 100 * (epi // 12)
  of='~/sgobs/garpos/simulator/obsdata/IMAG/IMAG.%s.kaiyo_k4-obs.csv' % ix
  cfgfile='./initcfg/IMAG/IMAG.%s.kaiyo_k4-initcfg.ini' % ix
  hd = "# cfgfile = %s\n" % cfgfile
  df.to_csv(of)
  os.system('sed -i -e "1i ' + hd + '" ' + of )
  basfile='~/sgobs/garpos/simulator/initcfg/IMAG/IMAG.1911.kaiyo_k4-initcfg.ini'
  cfsfile='~/sgobs/garpos/simulator/initcfg/IMAG/IMAG.%s.kaiyo_k4-initcfg.ini' % ix
  os.system('cp ' + basfile + ' ' + cfsfile)
  os.system("sed -i -e 's/1911.kaiyo_k4-obs/" + str(ix) + ".kaiyo_k4-obs/g' " + cfsfile )
  os.system("sed -i -e 's/Campaign    = 1911/Campaign    = " + str(ix) + "/g' " + cfsfile )
  os.system("sed -i -e 's/M01_dPos    =      0.0000    500.0000  -1000.0000/M01_dPos    =  " + str(round((slon1-slon)*dlon*1.e3,4)).rjust(10) + "  " + str(round((slat1-slat)*dlat*1.e3,4)).rjust(10) + "  " + str(round(-sdep1*ddep*1.e3,4)).rjust(10) + "/g' " + cfsfile )
  os.system("sed -i -e 's/M02_dPos    =    500.0000      0.0000  -1000.0000/M02_dPos    =  " + str(round((slon2-slon)*dlon*1.e3,4)).rjust(10) + "  " + str(round((slat2-slat)*dlat*1.e3,4)).rjust(10) + "  " + str(round(-sdep2*ddep*1.e3,4)).rjust(10) + "/g' " + cfsfile )
  os.system("sed -i -e 's/M03_dPos    =      0.0000   -500.0000  -1000.0000/M03_dPos    =  " + str(round((slon3-slon)*dlon*1.e3,4)).rjust(10) + "  " + str(round((slat3-slat)*dlat*1.e3,4)).rjust(10) + "  " + str(round(-sdep3*ddep*1.e3,4)).rjust(10) + "/g' " + cfsfile )
  os.system("sed -i -e 's/M04_dPos    =   -500.0000      0.0000  -1000.0000/M04_dPos    =  " + str(round((slon4-slon)*dlon*1.e3,4)).rjust(10) + "  " + str(round((slat4-slat)*dlat*1.e3,4)).rjust(10) + "  " + str(round(-sdep4*ddep*1.e3,4)).rjust(10) + "/g' " + cfsfile )
  os.system("sed -i -e 's/Center_ENU  =      0.0000      0.0000  -1000.0000/Center_ENU  =      0.0000      0.0000  " + str(round(-sdep4*ddep*1.e3,4)).rjust(10) + "/g' " + cfsfile )
  os.system("sed -i -e 's/N_shot      =   720/N_shot      =" + str(len(dkf)).rjust(6) + "/g' " + cfsfile )
  print(ix)
