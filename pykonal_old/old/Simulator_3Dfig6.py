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
import pykonal
import os,time,sys
from scipy.interpolate import Rbf
import colorednoise as cn

t1 = time.time() 
##
sv = pd.read_csv('svpC.csv')
##
# Initialize the solver.
#m_s_dp=[   3,   6,   30]
#m_s_km=[ -72.0, 504.0,  0.0] # cm/s/km

#m_s_dp=[  40,  80, 120]
m_s_dp=[  10,  20,  30]
m_s_km=[  0.0, 0.0,  0.0] # cm/s/km
#m_s_km=[ 108.0, 0.0,  0.0] # cm/s/km
#m_s_dp=[    9,   21,   30]
#m_s_km=[  0.0, 134.2,  0.0] # cm/s/km
#m_s_dp=[   36,   84,  120]
#m_s_km=[  0.0, 134.2,  0.0] # cm/s/km
#m_s_dp=[   15,   24,   30]
#m_s_km=[ -50.3, 307.5,  0.0] # cm/s/km
#m_s_dp=[   15,   24,   30]
#m_s_km=[ 118.5, -137.5,  0.0] # cm/s/km

#m_s_dp=[   10,   20,   30]
#m_s_km=[ 108.0,  0.0,  0.0] # cm/s/km
#m_s_dp=[    9,   21,   30]
#m_s_km=[  0.0, 150.0,  0.0] # cm/s/km
#m_s_dp=[   15,   24,   30]
#m_s_km=[-112.5, 687.5,  0.0] # cm/s/km
#m_s_dp=[   15,   24,   30]
#m_s_km=[ 307.5,-812.5,  0.0] # cm/s/km
#m_s_dpB=[   23,   24,   30]
#m_s_kmB=[ 86.9,-702.0,  0.0] # cm/s/km
#nlat,ndep,nlon=100,51,100#8000,8000,4000
#dlat,ddep,dlon=0.02,0.02,0.02
#slat,sdep,slon=50,50,50
nlat,ndep,nlon=200,101,200#8000,8000,4000
dlat,ddep,dlon=0.01,0.01,0.01
slat,sdep,slon=100,100,100
#nlat,ndep,nlon=400,201,400#8000,8000,4000
#dlat,ddep,dlon=0.005,0.005,0.005
#slat,sdep,slon=200,200,200
solverb = pykonal.EikonalSolver(coord_sys="cartesian")
solverc = pykonal.EikonalSolver(coord_sys="cartesian")

########### base model
solverb.velocity.min_coords = 0, 0, 0
solverb.velocity.node_intervals = dlat, ddep, dlon
# This time we want a 3D computational grid, so set the number of grid nodes
# in the z direction to 8 as well.
solverb.velocity.npts = nlat, ndep, nlon
#solver.velocity.values = np.ones(solver.velocity.npts)
#solver.velocity.values = np.full(solver.velocity.npts,1.480554)
#solverb.velocity.values = np.array([[[1.480554] * nlon] * ndep] * nlat)
svd = [[] for i in range(ndep)]
for i in range(ndep):
  svd[i] = np.array([(sv.speed[i]+10.)*0.001] * nlon)
svp = [svd[0]]
for i in range(1,ndep):
  svp = np.append(svp, [svd[i]], axis=0)
a0 = np.array([svp] * nlat)
solverb.velocity.values = a0
src_idx = slat, sdep, slon
solverb.traveltime.values[src_idx] = 0
solverb.unknown[src_idx] = False
solverb.trial.push(*src_idx)
# Solve the system.
solverb.solve()

########## com model
solverc.velocity.min_coords = 0, 0, 0
solverc.velocity.node_intervals = dlat, ddep, dlon
solverc.velocity.npts = nlat, ndep, nlon
svp  = [[] for s in range(nlat)]
for s in range(nlat):
  svc = [[] for i in range(ndep)]
  for i in range(0,m_s_dp[0]):
    svc[i] = np.array([m_s_km[0]*0.001*0.1*(1./float(nlat))*float(s-slat)] * nlon)
#    svc[i] = np.array([0] * nlon)
  for i in range(m_s_dp[0],m_s_dp[1]):
    svc[i] = np.array([m_s_km[1]*0.001*0.1*(1./float(nlat))*float(s-slat)] * nlon)
#    svc[i] = np.array([0] * nlon)
  for i in range(m_s_dp[1],ndep):
    svc[i] = np.array([m_s_km[2]*0.001*0.1*(1./float(nlat))*float(s-slat)] * nlon)
#    svc[i] = np.array([0] * nlon)
  svp[s] = [svc[0]]
  for i in range(1,ndep):
    svp[s] = np.append(svp[s], [svc[i]], axis=0)
a1 = [svp[0]]
for s in range(1,nlat):
  a1 = np.append(a1, [svp[s]], axis=0)
solverc.velocity.values = a0 + a1
# Initialize the source.
src_idx = slat, sdep, slon
solverc.traveltime.values[src_idx] = 0
solverc.unknown[src_idx] = False
solverc.trial.push(*src_idx)
# Solve the system.
solverc.solve()

# And finally, get the traveltime values.
#print(solver.traveltime.values)
### Plot the results
plt.close('all')
fig = plt.figure(figsize=(5, 7.5))
ax = fig.add_subplot(1, 1, 1, frameon=False)
plt.tick_params(labelcolor='none', top=False, bottom=False, left=False, right=False)
ax.set_xlabel("Eastward offset [km]")
ax.set_ylabel("Northward offset [km]                             Downward offset [km]")
gridspec = gs.GridSpec(nrows=4, ncols=1, height_ratios=[0.04,0.32,0.06,1.05])
cax00 = fig.add_subplot(gridspec[0, 0])
cax01 = fig.add_subplot(gridspec[2, 0])
ax10 = fig.add_subplot(gridspec[1, 0])
ax11 = fig.add_subplot(gridspec[3, 0])
panel = ord("a")

for ax in (ax10, ax11):
    ax.text(-0.05, 0.8, f"({chr(panel)})", ha="right", va="top", transform=ax.transAxes)
    panel += 1
ttm = solverb.traveltime.values # - solverb.traveltime.values
ttc = solverc.traveltime.values # - solverb.traveltime.values
tt0 = solverc.traveltime.values - solverb.traveltime.values
#tt1 = solverd.traveltime.values - solverc.traveltime.values
qmesh = ax10.pcolormesh(
    solverc.velocity.nodes[:,:,0,0],
    solverc.velocity.nodes[:,:,0,1],
    a1[:,:,int(nlon/2)]*10.,
    cmap=plt.get_cmap("seismic"),
    norm=Normalize(vmin=-0.1, vmax=0.1)
)
ax10.set_ylim(0,1)
ax10.invert_yaxis()
cbar = fig.colorbar(qmesh, cax=cax00, orientation="horizontal", ticks=[-0.1,-0.05,0,0.05,0.1])
cbar.set_label("d-Velocity [m/s]")
cbar.ax.xaxis.tick_top()
cbar.ax.xaxis.set_label_position("top")
cbar.set_clim(-0.1,0.1)
qmesh = ax11.pcolormesh(
    solverb.traveltime.nodes[:,0,:,0],
    solverb.traveltime.nodes[:,0,:,2],
#    tt1[:,0,:],
    tt0[:,0,:],
    cmap=plt.get_cmap("coolwarm"),
    norm=Normalize(vmin=-0.005,vmax=0.005)
)
ax11.contour(
    solverb.traveltime.nodes[:,0,:,0],
    solverb.traveltime.nodes[:,0,:,2],
#    tt1[:,0,:],
    tt0[:,0,:],
    colors="k",
    levels=np.arange(tt0.min(), tt0.max(), 0.0003),
    linewidths=0.2,
    linestyles="--"
)
#print(tt0[:,0,50])
np.savetxt('centerline.csv',tt0[:,0,slon],delimiter=',')
#np.savetxt('line1500.csv',tt0[:,0,65],delimiter=',')
ax11.plot([1, 1], [0, 2], 'w-', lw=2)
ax11.plot([0, 2], [1, 1], 'w-', lw=2)
ax11.plot([1-1/math.sqrt(2), 1+1/math.sqrt(2)], [1-1/math.sqrt(2), 1+1/math.sqrt(2)], 'w-', lw=2)
ax11.plot([1-1/math.sqrt(2), 1+1/math.sqrt(2)], [1+1/math.sqrt(2), 1-1/math.sqrt(2)], 'w-', lw=2)
ax11.plot([1, 2, 1, 0, 1], [0, 1, 2, 1, 0], 'w-', lw=2)
ax11.plot([1, 1.5, 1, 0.5, 1], [0.5, 1, 1.5, 1, 0.5], 'w-', lw=2)
ax11.plot(1, 1.0, marker='*',markersize=4)
ax11.set_xlim([0,2])
ax11.set_ylim([0,2])
cbar = fig.colorbar(qmesh, cax=cax01, orientation="horizontal", ticks=[-0.005,0.0,0.005])
#cbar = fig.colorbar(qmesh, cax=cax01, orientation="horizontal", ticks=[-0.01,-0.005,0,0.005,0.01])
cbar.set_label("")
cbar.set_clim(-0.005,0.005)
#cbar.ax.xaxis.tick_top()
#cbar.ax.xaxis.set_label_position("top")
plt.savefig('figure.png')
shutil.copyfile("figure.png", "/mnt/owncloud_webdav/webdav/figureD8.png")
#################################
#################################
#xn=solverb.traveltime.nodes[50:150,0:4,50:150,0].reshape(-1,)
#yn=solverb.traveltime.nodes[50:150,0:4,50:150,1].reshape(-1,)
#zn=solverb.traveltime.nodes[50:150,0:4,50:150,2].reshape(-1,)
xn=solverb.traveltime.nodes[:,0,slon,0].reshape(-1,)
yn=solverb.traveltime.nodes[:,0,slon,1].reshape(-1,)
zn=solverb.traveltime.nodes[:,0,slon,2].reshape(-1,)
tn=ttc[:,0,slon].reshape(-1,)
t0=np.sqrt((xn-1)**2+(yn-1)**2+(zn-1)**2)/1.51
df = pd.DataFrame({'a':xn,'b':yn,'c':zn,'d':tn,'e':t0})
df.to_csv('test_li_00.csv')
shutil.copyfile("test_li_00.csv", "/mnt/owncloud_webdav/webdav/test_li_11.csv")
stop
rbfi = Rbf(xn,yn,zn,tn)
#################################
lhptb  = 0.0 #5000 #(km-order)
lnptb  = 0.0 #0250 #(km-order)
lzptb  = 0.0 #0250 #(km-order)
ghptb  = 0.0 #0001 #(km-order)
gvptb  = 0.0 #0003 #(km-order)
#################################
lknot  = 7.0*1.852
beta10 = 1.0
beta15 = 1.5
beta20 = 2.0
#################################
ld     = pd.read_csv('linesample.csv')
for n in range(len(ld['se'])):
  lstart = np.array([ld['se'][n],ld['sn'][n],ld['su'][n]])
  lend   = np.array([ld['ee'][n],ld['en'][n],ld['eu'][n]])
  line   = lend-lstart
  kyori  = np.linalg.norm(line)
  ltime  = int(3600.*kyori/lknot)+30
  gpe = cn.powerlaw_psd_gaussian(beta15, ltime)
  gpe = gpe - np.average(gpe)
  gpn = cn.powerlaw_psd_gaussian(beta15, ltime)
  gpn = gpn - np.average(gpn)
  gpu = cn.powerlaw_psd_gaussian(beta15, ltime)
  gpu = gpu - np.average(gpu)
  l_e = cn.powerlaw_psd_gaussian(beta20, ltime)
  l_e = l_e - np.average(l_e)
  l_n = cn.powerlaw_psd_gaussian(beta20, ltime)
  l_n = l_n - np.average(l_n)
  l_z = cn.powerlaw_psd_gaussian(beta10, ltime)
  l_z = l_z - np.average(l_z)
#np.savetxt('test_dit.csv',l_h,delimiter=',')
#################################
  ang = math.atan2(ld['ee'][n]-ld['se'][n],ld['en'][n]-ld['sn'][n])
  sewr =                   (np.cos(ang)*lhptb+np.sin(ang)*lnptb)*l_e[::10]+np.linspace(lstart[0],lend[0],len(l_e[::10]))
  sewn = ghptb*gpe[::10] + sewr
  snsr =                   (np.sin(ang)*lhptb+np.cos(ang)*lnptb)*l_n[::10]+np.linspace(lstart[1],lend[1],len(l_e[::10]))
  snsn = ghptb*gpn[::10] + snsr
  sudr =                                                   lzptb*l_z[::10]+np.linspace(lstart[2],lend[2],len(l_e[::10]))
  sudn = gvptb*gpu[::10] + sudr
  sdi = rbfi(sewr, sudr, snsr)
#################################
  rewr = sewr + sdi*lknot*(ld['ee'][n]-ld['se'][n])/(3600.*kyori) + (np.cos(ang)*lhptb+np.sin(ang)*lnptb)*l_e[3::10]
  rewn = ghptb*gpe[3::10]+ rewr
  rnsr = snsr + sdi*lknot*(ld['en'][n]-ld['sn'][n])/(3600.*kyori) + (np.sin(ang)*lhptb+np.cos(ang)*lnptb)*l_n[3::10]
  rnsn = ghptb*gpn[3::10]+ rnsr
  rudr = lzptb*l_z[3::10]+np.linspace(lstart[2],lend[2],len(l_e[::10]))
  rudn = gvptb*gpu[3::10]+ rudr
  rdi  = np.sqrt((sewr-1)**2+(snsr-1)**2+(sudr-1)**2)/1.5
#  rdi = rbfi(rewr, rudr, rnsr)
#################################
  df = pd.DataFrame({'sEW': sewr,'sNS': snsr,'sUD': sudr,'sEWn': sewn,'sNSn': snsn,'sUDn': sudn,'sTT': sdi,'rEW': rewr,'rNS': rnsr,'rUD': rudr,'rEWn': rewn,'rNSn': rnsn,'rUDn': rudn,'rTT': rdi})
#  df = pd.DataFrame({'sEW': sewr,'sNS': snsr,'sUD': sudr,'sEWn': sewn,'sNSn': snsn,'sUDn': sudn,'sTT': sdi,'rEW': rewr,'rNS': rnsr,'rUD': rudr,'rEWn': rewn,'rNSn': rnsn,'rUDn': rudn,'rTT': rdi})
  df.to_csv('test_li_%s.csv'%(str(n).zfill(2)))
shutil.copyfile("test_li_00.csv", "/mnt/owncloud_webdav/webdav/test_li_00.csv")
