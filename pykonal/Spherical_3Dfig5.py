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

sv = pd.read_csv('svpC.csv')
sva = pd.read_csv('svpA.csv')


Erad=6371.
m_s_dp=[   10,   13,   16]
m_s_km=[ 40.0,  0.0,  0.0] # cm/s/km
#m_s_dp=[   10,   20,   30]
#m_s_km=[ 1.08,  0.0,  0.0] # cm/s/km
#m_s_dp=[    3,    6,   30]
#m_s_km=[-0.60, 4.20,  0.0] # cm/s/km
nrh,nth,nph=50,100,100
drh,dth,dph=0.1, (10./Erad)/float(nth-1), (10./Erad)/float(nth-1)
srh,sth,sph=30,50,50

vv = pykonal.fields.ScalarField3D(coord_sys='spherical')
dv = pykonal.fields.ScalarField3D(coord_sys='spherical')
vv.min_coords = dv.min_coords = [Erad-5., np.pi/2, np.pi/2]
vv.node_intervals = dv.node_intervals = [drh,dth,dph]
vv.npts = dv.npts = [nrh,nth,nph]
v2 = pykonal.fields.ScalarField3D(coord_sys='spherical')
d2 = pykonal.fields.ScalarField3D(coord_sys='spherical')
v2.min_coords = d2.min_coords = [Erad-5., np.pi/2, np.pi/2]
v2.node_intervals = d2.node_intervals = [drh,dth,dph]
v2.npts = d2.npts = [nrh,nth,nph]

a0 = [[] for j in range(nrh)]
#svp = [[] for j in range(nrh)]
for j in range(nrh):
  svd = [[] for i in range(nth)]
  for i in range(nth):
    svd[i] = np.array([sv.speed[j]*0.001] * nph)
  a0[j] = [svd[0]]
  for i in range(1,nth):
    a0[j] = np.append(a0[j], [svd[i]], axis=0)
###
a1 = [[] for j in range(nrh)]
#svp = [[] for j in range(nrh)]
for j in range(0,m_s_dp[0]):
  svd = [[] for i in range(nth)]
  for i in range(nth):
    for p in range(nph):
      a      = [sv.speed[j]*0.001+m_s_km[0]*0.001*float(p-50)]
      svd[i] = np.append(svd[i], a)
  a1[j] = [svd[0]]
  for i in range(1,nth):
    a1[j] = np.append(a1[j], [svd[i]], axis=0)
for j in range(m_s_dp[0],m_s_dp[1]):
  svd = [[] for i in range(nth)]
  for i in range(nth):
    for p in range(nph):
      a      = [sv.speed[j]*0.001+m_s_km[1]*0.001*float(p-50)]
      svd[i] = np.append(svd[i], a)
  a1[j] = [svd[0]]
  for i in range(1,nth):
    a1[j] = np.append(a1[j], [svd[i]], axis=0)
for j in range(m_s_dp[1],nrh):
  svd = [[] for i in range(nth)]
  for i in range(nth):
    for p in range(nph):
      a      = [sv.speed[j]*0.001+m_s_km[2]*0.001*float(p-50)]
      svd[i] = np.append(svd[i], a)
  a1[j] = [svd[0]]
  for i in range(1,nth):
    a1[j] = np.append(a1[j], [svd[i]], axis=0)
vv.values = a0
dv.values = a0
v2.values = a1
d2.values = a1

# Resample the velocity model (64, 1, 128) --> (1024, 1, 2048)
rho_min, theta_min, phi_min = vv.min_coords
rho_max, theta_max, phi_max = vv.max_coords
nrho, ntheta, nphi = nrh, nth, nph

drho = (rho_max - rho_min) / (nrho - 1)
rho = np.linspace(rho_min, rho_max, nrho)

dtheta = (theta_max - theta_min) / (ntheta - 1)
theta = np.linspace(theta_min, theta_max, ntheta)

dphi = (phi_max - phi_min) / (nphi - 1)
phi = np.linspace(phi_min, phi_max, nphi)

rtp = np.moveaxis(
    np.stack(np.meshgrid(rho, theta, phi, indexing="ij")),
    0, 
    -1
)
###
vv_new = vv.resample(rtp.reshape(-1, 3)).reshape(rtp.shape[:-1])
dv_new = dv.resample(rtp.reshape(-1, 3)).reshape(rtp.shape[:-1])
vv = pykonal.fields.ScalarField3D(coord_sys="spherical")
dv = pykonal.fields.ScalarField3D(coord_sys="spherical")
vv.min_coords = dv.min_coords = rho_min, theta_min, phi_min
vv.node_intervals = dv.node_intervals = drho, dtheta, dphi
vv.npts = dv.npts = nrho, ntheta, nphi
vv.values = vv_new
dv.values = dv_new
velocity = vv
###
v2_new = v2.resample(rtp.reshape(-1, 3)).reshape(rtp.shape[:-1])
d2_new = d2.resample(rtp.reshape(-1, 3)).reshape(rtp.shape[:-1])
v2 = pykonal.fields.ScalarField3D(coord_sys="spherical")
d2 = pykonal.fields.ScalarField3D(coord_sys="spherical")
v2.min_coords = d2.min_coords = rho_min, theta_min, phi_min
v2.node_intervals = d2.node_intervals = drho, dtheta, dphi
v2.npts = d2.npts = nrho, ntheta, nphi
v2.values = v2_new
d2.values = d2_new
velocity2 = v2
###

SRC_IDX = np.array([srh, sth, sph])

traveltime_fields = dict()
traveltime2_fields = dict()
for decimation_factor in range(7, -1, -1):
    decimation_factor = 2**decimation_factor
    
    vv = velocity.values[::decimation_factor, :, ::decimation_factor]
    solver = pykonal.EikonalSolver(coord_sys="spherical")
    solver.vv.min_coords = velocity.min_coords
    solver.vv.node_intervals = velocity.node_intervals * decimation_factor
    solver.vv.npts = vv.shape
    solver.vv.values = vv
    src_idx = tuple(SRC_IDX // decimation_factor - [1, 0, 1])
    solver.traveltime.values[src_idx] = 0
    solver.unknown[src_idx] = False
    solver.trial.push(*src_idx)
    solver.solve()
    traveltime_fields[decimation_factor] = solver.traveltime
    
    v2 = velocity2.values[::decimation_factor, :, ::decimation_factor]
    solver2 = pykonal.EikonalSolver(coord_sys="spherical")
    solver2.vv.min_coords = velocity2.min_coords
    solver2.vv.node_intervals = velocity2.node_intervals * decimation_factor
    solver2.vv.npts = v2.shape
    solver2.vv.values = v2
    src_idx = tuple(SRC_IDX // decimation_factor - [1, 0, 1])
    solver2.traveltime.values[src_idx] = 0
    solver2.unknown[src_idx] = False
    solver2.trial.push(*src_idx)
    solver2.solve()
    traveltime2_fields[decimation_factor] = solver2.traveltime

plt.close('all')
fig = plt.figure(figsize=(5, 7.5))

fig.text(0, 0.5, "$y$ [km]", ha="left", va="center", rotation=90)
fig.text(0.5, 0.05, "$x$ [km]", ha="center", va="bottom")

gridspec = gs.GridSpec(nrows=4, ncols=1, height_ratios=[0.04,0.32,0.06,1.05])
#gridspec = gs.GridSpec(nrows=4, ncols=1, height_ratios=[0.03,0.82,0.03,1.64])
cax00 = fig.add_subplot(gridspec[0, 0])
cax01 = fig.add_subplot(gridspec[2, 0])
ax10 = fig.add_subplot(gridspec[1, 0])
ax11 = fig.add_subplot(gridspec[3, 0])

panel = ord("a")
for ax in (ax10, ax11):
    ax.text(-0.05, 1.1, f"({chr(panel)})", ha="right", va="top", transform=ax.transAxes)
    panel += 1

tt0 = traveltime_fields[1]
tt1 = traveltime2_fields[1]
nodes = tt0.nodes
xx = Erad-nodes[...,0] * np.sin(nodes[...,1]) * np.cos(nodes[...,2])
yy = Erad-nodes[...,0] * np.sin(nodes[...,1]) * np.sin(nodes[...,2])-2.
#print(nodes[...,0]*(1.-np.sin(nodes[...,1])))
#print(xx-nodes[...,0]*(1.-np.sin(nodes[...,1])))
zz = nodes[...,0]*(1.-np.sin(nodes[...,1]))
#print(tt0.values[:,0])
ax10.set_ylim(0,3)
qmesh = ax10.pcolormesh(
    xx[:,0],
    yy[:,0],
    v2_new[:,0,:]-vv_new[:,0,:],
    cmap=plt.get_cmap("seismic"),
    shading="gouraud",
    norm=Normalize(vmin=-0.1, vmax=0.1)
)
cbar = fig.colorbar(qmesh, cax=cax00, orientation="horizontal", ticks=[-0.1,-0.05,0,0.05,0.1])
cbar.set_label("d-Velocity [m/s]")
cbar.ax.xaxis.tick_top()
cbar.ax.xaxis.set_label_position("top")
cbar.set_clim(-0.1,0.1)

x1 = nodes[...,2]*Erad-nodes[0,0,0,2]*Erad
y1 = nodes[...,1]*Erad-nodes[0,0,0,2]*Erad
qmesh = ax11.pcolormesh(
    x1[0,:],
    y1[0,:],
    tt1.values[0,:,:]-tt0.values[0,:,:],
    cmap=plt.get_cmap("coolwarm"),
    shading="gouraud",
    norm=Normalize(vmin=-0.02, vmax=0.02)
)
ax11.contour(
    x1[0,:],
    y1[0,:],
    tt1.values[0,:,:]-tt0.values[0,:,:],
    colors="k",
    levels=np.arange(-0.1, 0.1, 0.001),
    linewidths=0.2,
    linestyles="--"
)
cbar = fig.colorbar(qmesh, cax=cax01, orientation="horizontal", ticks=[-0.02,-0.01,0.01,0.02])
cbar.set_label("Traveltime [s]")
cbar.ax.xaxis.tick_top()
cbar.ax.xaxis.set_label_position("top")

ax10.set_xticklabels([])
ax11.plot([5, 5], [2, 8], 'w-', lw=2)
ax11.plot([2, 8], [5, 5], 'w-', lw=2)
ax11.plot([5-3/math.sqrt(2), 5+3/math.sqrt(2)], [5-3/math.sqrt(2), 5+3/math.sqrt(2)], 'w-', lw=2)
ax11.plot([5-3/math.sqrt(2), 5+3/math.sqrt(2)], [5+3/math.sqrt(2), 5-3/math.sqrt(2)], 'w-', lw=2)
ax11.plot([5, 8, 5, 2, 5], [2, 5, 8, 5, 2], 'w-', lw=2)
ax11.plot([5, 6.5, 5, 3.5, 5], [3.5, 5, 6.5, 5, 3.5], 'w-', lw=2)
ax11.plot(5, 5, marker='*',markersize=4)
cbar.set_clim(-0.02,0.02)
plt.savefig('figure.png')
shutil.copyfile("figure.png", "/mnt/owncloud_webdav/webdav/figureS.png")
