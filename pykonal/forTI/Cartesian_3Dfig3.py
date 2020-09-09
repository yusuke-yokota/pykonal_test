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

m_s_dp=[   10,   20,   30]
m_s_km=[ 108.0, 0.0,  0.0] # cm/s/km
#m_s_dp=[    9,   21,   30]
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
nlat,ndep,nlon=100,50,100#8000,8000,4000
dlat,ddep,dlon=0.1,0.1,0.1
slat,sdep,slon=50,30,50
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
  svd[i] = np.array([sv.speed[i]*0.001] * nlon)
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
    svc[i] = np.array([m_s_km[0]*0.001*0.01*0.1*float(s-50)] * nlon)
#    svc[i] = np.array([0] * nlon)
  for i in range(m_s_dp[0],m_s_dp[1]):
    svc[i] = np.array([m_s_km[1]*0.001*0.01*0.1*float(s-50)] * nlon)
#    svc[i] = np.array([0] * nlon)
  for i in range(m_s_dp[1],ndep):
    svc[i] = np.array([m_s_km[2]*0.001*0.01*0.1*float(s-50)] * nlon)
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
#tt0 = solverb.traveltime.values # - solverb.traveltime.values
tt0 = solverc.traveltime.values - solverb.traveltime.values
#tt1 = solverd.traveltime.values - solverc.traveltime.values
qmesh = ax10.pcolormesh(
    solverc.velocity.nodes[:,:,0,0],
    solverc.velocity.nodes[:,:,0,1],
    a1[:,:,int(nlon/2)]*1000.,
    cmap=plt.get_cmap("seismic"),
    norm=Normalize(vmin=-10, vmax=10)
)
ax10.set_ylim(0,3)
ax10.invert_yaxis()
cbar = fig.colorbar(qmesh, cax=cax00, orientation="horizontal", ticks=[-10,-5,0,5,10])
cbar.set_label("d-Velocity [m/s]")
cbar.ax.xaxis.tick_top()
cbar.ax.xaxis.set_label_position("top")
cbar.set_clim(-10,10)
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
#np.savetxt('centerline.csv',tt0[:,0,50],delimiter=',')
#np.savetxt('line1500.csv',tt0[:,0,65],delimiter=',')
ax11.plot([5, 5], [2, 8], 'w-', lw=2)
ax11.plot([2, 8], [5, 5], 'w-', lw=2)
ax11.plot([5-3/math.sqrt(2), 5+3/math.sqrt(2)], [5-3/math.sqrt(2), 5+3/math.sqrt(2)], 'w-', lw=2)
ax11.plot([5-3/math.sqrt(2), 5+3/math.sqrt(2)], [5+3/math.sqrt(2), 5-3/math.sqrt(2)], 'w-', lw=2)
ax11.plot([5, 8, 5, 2, 5], [2, 5, 8, 5, 2], 'w-', lw=2)
ax11.plot([5, 6.5, 5, 3.5, 5], [3.5, 5, 6.5, 5, 3.5], 'w-', lw=2)
ax11.plot(5, 5.0, marker='*',markersize=4)
cbar = fig.colorbar(qmesh, cax=cax01, orientation="horizontal", ticks=[-0.005,0.0,0.005])
cbar.set_label("")
cbar.set_clim(-0.005,0.005)
#cbar.ax.xaxis.tick_top()
#cbar.ax.xaxis.set_label_position("top")
plt.savefig('figure.png')
#shutil.copyfile("figure.png", "/mnt/owncloud_webdav/webdav/figure.png")
