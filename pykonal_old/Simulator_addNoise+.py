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
tsv=0.
ai=0
for i in range(1,ndep):
  tsv = tsv+float(sv.speed[i])
tov = np.arange(ndep)
tsv=(tsv+(sv.speed[ndep]+sv.speed[0])*0.5)/float(ndep)
tov = np.full(ndep, tsv, dtype="float64")

solverb1 = pykonal.EikonalSolver(coord_sys="cartesian")
solverc1 = pykonal.EikonalSolver(coord_sys="cartesian")
solverb2 = pykonal.EikonalSolver(coord_sys="cartesian")
solverc2 = pykonal.EikonalSolver(coord_sys="cartesian")
solverb3 = pykonal.EikonalSolver(coord_sys="cartesian")
solverc3 = pykonal.EikonalSolver(coord_sys="cartesian")
solverb4 = pykonal.EikonalSolver(coord_sys="cartesian")
solverc4 = pykonal.EikonalSolver(coord_sys="cartesian")

########### base model
# This time we want a 3D computational grid, so set the number of grid nodes
# in the z direction to 8 as well.
solverb1.velocity.min_coords = 0, 0, 0
solverb1.velocity.node_intervals = dlon, ddep, dlat
solverb1.velocity.npts = nlon, ndep, nlat
solverb2.velocity.min_coords = 0, 0, 0
solverb2.velocity.node_intervals = dlon, ddep, dlat
solverb2.velocity.npts = nlon, ndep, nlat
solverb3.velocity.min_coords = 0, 0, 0
solverb3.velocity.node_intervals = dlon, ddep, dlat
solverb3.velocity.npts = nlon, ndep, nlat
solverb4.velocity.min_coords = 0, 0, 0
solverb4.velocity.node_intervals = dlon, ddep, dlat
solverb4.velocity.npts = nlon, ndep, nlat
solverc1.velocity.min_coords = 0, 0, 0
solverc1.velocity.node_intervals = dlon, ddep, dlat
solverc1.velocity.npts = nlon, ndep, nlat
solverc2.velocity.min_coords = 0, 0, 0
solverc2.velocity.node_intervals = dlon, ddep, dlat
solverc2.velocity.npts = nlon, ndep, nlat
solverc3.velocity.min_coords = 0, 0, 0
solverc3.velocity.node_intervals = dlon, ddep, dlat
solverc3.velocity.npts = nlon, ndep, nlat
solverc4.velocity.min_coords = 0, 0, 0
solverc4.velocity.node_intervals = dlon, ddep, dlat
solverc4.velocity.npts = nlon, ndep, nlat
####################

svdo = [[] for i in range(ndep)]
svd = [[] for i in range(ndep)]
for i in range(ndep):
  svdo[i] = np.array([(tov[i])*0.0010000] * nlat)
  svd[i] = np.array([(sv.speed[i])*0.0010000] * nlat)
svpo = [svdo[0]]
svp = [svd[0]]
for i in range(1,ndep):
  svpo = np.append(svpo, [svdo[i]], axis=0)
  svp = np.append(svp, [svd[i]], axis=0)
o0 = np.array([svpo] * nlon)
a0 = np.array([svp] * nlon)
########## com model
svp  = [[] for s in range(nlon)]
chet = np.arange(o_chet+1, dtype=float)
for s in range(nlon):
  achet=0
  for c in range(o_chet+1):
    chet[c] = eval_chebyt(c, (s-slon)/slon)
    achet = c_chet[c]*chet[c] + achet
  achet = achet/d_chet
  svc = [[] for i in range(ndep)]
  for i in range(0,m_s_dp[0]):
    svc[i] = np.array([m_s_km[0]*0.001*0.1*0.5*achet] * nlat)
#    svc[i] = np.array([m_s_km[0]*0.001*0.1*(1./float(2*slon))*float(s-slon)] * nlat)
  for i in range(m_s_dp[0],m_s_dp[1]):
    svc[i] = np.array([m_s_km[1]*0.001*0.1*0.5*achet] * nlat)
#    svc[i] = np.array([m_s_km[1]*0.001*0.1*(1./float(2*slon))*float(s-slon)] * nlat)
  for i in range(m_s_dp[1],ndep):
    svc[i] = np.array([m_s_km[2]*0.001*0.1*0.5*achet] * nlat)
#    svc[i] = np.array([m_s_km[2]*0.001*0.1*(1./float(2*slon))*float(s-slon)] * nlat)
  svp[s] = [svc[0]]
  for i in range(1,ndep):
    svp[s] = np.append(svp[s], [svc[i]], axis=0)
a1 = [svp[0]]
for s in range(1,nlon):
  a1 = np.append(a1, [svp[s]], axis=0)
a1=a1+a0
s01=solverc1.velocity.nodes[slon1,sdep1,slat1,0]
s02=solverc1.velocity.nodes[slon1,sdep1,slat1,2]
s03=solverc1.velocity.nodes[slon1,sdep1,slat1,1]
s11=solverc2.velocity.nodes[slon2,sdep2,slat2,0]
s12=solverc2.velocity.nodes[slon2,sdep2,slat2,2]
s13=solverc2.velocity.nodes[slon2,sdep2,slat2,1]
s21=solverc3.velocity.nodes[slon3,sdep3,slat3,0]
s22=solverc3.velocity.nodes[slon3,sdep3,slat3,2]
s23=solverc3.velocity.nodes[slon3,sdep3,slat3,1]
s31=solverc4.velocity.nodes[slon4,sdep4,slat4,0]
s32=solverc4.velocity.nodes[slon4,sdep4,slat4,2]
s33=solverc4.velocity.nodes[slon4,sdep4,slat4,1]
####################
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
    nin = int(dkf.ant_n0[i]/dlat)
    nie = int(dkf.ant_e0[i]/dlon)
    if dkf.MT[i] == 'M01':
      dkf.ResiTT[i] = dkf.TT[i] + 2.*(np.sqrt((dkf.ant_e1[i]-s01)**2 + (dkf.ant_n1[i]-s02)**2 + (dkf.ant_u1[i]-s03)**2) - np.sqrt((dkf.ant_e0[i]-s01)**2 + (dkf.ant_n0[i]-s02)**2 + (dkf.ant_u0[i]-s03)**2)) / a1[nin,0,nie]
    elif dkf.MT[i] == 'M02':
      dkf.ResiTT[i] = dkf.TT[i] + 2.*(np.sqrt((dkf.ant_e1[i]-s11)**2 + (dkf.ant_n1[i]-s12)**2 + (dkf.ant_u1[i]-s13)**2) - np.sqrt((dkf.ant_e0[i]-s11)**2 + (dkf.ant_n0[i]-s12)**2 + (dkf.ant_u0[i]-s13)**2)) / a1[nin,0,nie]
    elif dkf.MT[i] == 'M03':
      dkf.ResiTT[i] = dkf.TT[i] + 2.*(np.sqrt((dkf.ant_e1[i]-s21)**2 + (dkf.ant_n1[i]-s22)**2 + (dkf.ant_u1[i]-s23)**2) - np.sqrt((dkf.ant_e0[i]-s21)**2 + (dkf.ant_n0[i]-s22)**2 + (dkf.ant_u0[i]-s23)**2)) / a1[nin,0,nie]
    elif dkf.MT[i] == 'M04':
      dkf.ResiTT[i] = dkf.TT[i] + 2.*(np.sqrt((dkf.ant_e1[i]-s31)**2 + (dkf.ant_n1[i]-s32)**2 + (dkf.ant_u1[i]-s33)**2) - np.sqrt((dkf.ant_e0[i]-s31)**2 + (dkf.ant_n0[i]-s32)**2 + (dkf.ant_u0[i]-s33)**2)) / a1[nin,0,nie]
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
