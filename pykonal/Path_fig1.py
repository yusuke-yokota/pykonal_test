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

def align_marker(marker, halign='center', valign='middle',):
    # Define the base marker
    bm = markers.MarkerStyle(marker)

    # Get the marker path and apply the marker transform to get the
    # actual marker vertices (they should all be in a unit-square
    # centered at (0, 0))
    m_arr = bm.get_path().transformed(bm.get_transform()).vertices

    # Shift the marker vertices for the specified alignment.
    m_arr[:, 0] += halign / 2
    m_arr[:, 1] += valign / 2

    return Path(m_arr, bm.get_path().codes)

#Set the velocity gradient in 1/s
#sv = pd.read_csv('svpC.csv')
sv = pd.read_csv('svpA.csv')

# Initialize the solver.
m_s_dp=[  400,  800, 1200]
m_s_km=[  10.,  0.0,  0.0] # cm/s/km
#m_s_dp=[    3,    6,   30]
#m_s_km=[-0.60, 4.20,  0.0] # cm/s/km
#m_s_dp=[   15,   24,   30]
#m_s_km=[ 0.87,-0.25,  0.0] # cm/s/km
velocity_gradient = 1/1.5

nlat,ndep,nlon=4000,2000,1 #4000#8000,8000,4000
dlat,ddep,dlon=0.0025,0.0025,0.1
slat,sdep,slon=2000,1200,1 #2000
####
svX= [[] for i in range(ndep)]
for i in range(len(sv.speed)-1):
  for t in range(int(ndep/(len(sv.speed)-1))):
    svX[i*int(ndep/(len(sv.speed)-1))+t]=(sv.speed[i]*(int(ndep/(len(sv.speed)-1))-t)+sv.speed[i+1]*t)/(ndep/(len(sv.speed)-1))
####
svd = [[] for i in range(ndep)]
for i in range(ndep):
  svd[i] = np.array([svX[i]*0.001] * nlon)
svp = [svd[0]]
for i in range(1,ndep):
  svp = np.append(svp, [svd[i]], axis=0)
a0 = np.array([svp] * nlat)
####
svp  = [[] for s in range(nlat)]
for s in range(nlat):
  svc = [[] for i in range(ndep)]
  for i in range(0,m_s_dp[0]):
    svc[i] = np.array([m_s_km[0]*0.01*0.001] * nlon) #float(s-1000)] * nlon)
#    svc[i] = np.array([0] * nlon)
  for i in range(m_s_dp[0],m_s_dp[1]):
    svc[i] = np.array([m_s_km[1]*0.01*0.001] * nlon) #float(s-1000)] * nlon)
#    svc[i] = np.array([0] * nlon)
  for i in range(m_s_dp[1],m_s_dp[2]):
    svc[i] = np.array([m_s_km[2]*0.01*0.001] * nlon) #float(s-1000)] * nlon)
#    svc[i] = np.array([0] * nlon)
  for i in range(m_s_dp[2],ndep):
    svc[i] = np.array([m_s_km[2]*0.01*0.001] * nlon) #float(s-1000)] * nlon)
#    svc[i] = np.array([0] * nlon)
  svp[s] = [svc[0]]
  for i in range(1,ndep):
    svp[s] = np.append(svp[s], [svc[i]], axis=0)
a1 = [svp[0]]
for s in range(1,nlat):
  a1 = np.append(a1, [svp[s]], axis=0)
####

velocity = pykonal.fields.ScalarField3D(coord_sys="cartesian")
velocity.min_coords = 0, 0, 0
velocity.npts = nlat,ndep,1 #4096, 2048, 1
velocity.node_intervals = [10, 5, 1] / velocity.npts
velocity.values = a0  + a1
#velocity.values = np.full(velocity.npts, fill_value=1.5)

#for iy in range(velocity.npts[1]):
#    velocity.values[:,iy] += velocity_gradient * velocity.nodes[0,iy,0,1]


# Set the take-off angle in radians.
takeoff_angle = np.radians(70)

# Set the index of the node corresponding to the source location.
src_idx = np.array([slat, sdep, 0])

# Get the source coordinates.
src_coords = velocity.nodes[tuple(src_idx)]

# Compute the radius of the circle defining the raypath.
#print(velocity.value(src_coords),velocity_gradient,np.sin(takeoff_angle))
radius = velocity.value(src_coords) / (velocity_gradient * np.sin(takeoff_angle))

# Compute the coordinates of the center of the circle.
center_coords = src_coords + radius*np.array([np.cos(takeoff_angle), np.sin(takeoff_angle), 0])

# Define function mapping horizontal coordinate to vertical
def vertical_coords(horizontal_coord, src_coords=src_coords, radius=radius, takeoff_angle=takeoff_angle):
    sqrt = np.sqrt(radius**2 - (horizontal_coord - src_coords[0] - radius*np.cos(takeoff_angle))**2)
    first_two_terms = src_coords[1]  -  radius * np.sin(takeoff_angle)
    return (first_two_terms + sqrt, first_two_terms - sqrt)

# Define function mapping vertical coordinate to horizontal
def horizontal_coords(vertical_coord, src_coords=src_coords, radius=radius, takeoff_angle=takeoff_angle):
    sqrt = np.sqrt(radius**2 - (src_coords[1] - vertical_coord - radius*np.sin(takeoff_angle))**2)
    first_two_terms = src_coords[0]  +  radius * np.cos(takeoff_angle)
    return (first_two_terms + sqrt, first_two_terms - sqrt)

rec_coords = [horizontal_coords(0)[0], 0]

rays = dict()

for decimation_factor in range(4, 1, -1):
    decimation_factor = 2**decimation_factor
    
    vv = velocity.values[::decimation_factor, ::decimation_factor]

    solver = pykonal.EikonalSolver(coord_sys="cartesian")

    solver.velocity.min_coords = 0, 0, 0
    solver.velocity.node_intervals = velocity.node_intervals * decimation_factor
    solver.velocity.npts = vv.shape
    solver.velocity.values = vv

    idx = tuple((src_idx / decimation_factor).astype(np.int) - [1, 1, 0])
    solver.traveltime.values[idx] = 0
    solver.unknown[idx] = False
    solver.trial.push(*idx)

    solver.solve()
    rays[decimation_factor] = solver.trace_ray(np.array([rec_coords[0], 0, 0]))

zoom = 22
dx = 0.5
dy = 0.5
anchor_y = 1.1


plt.close("all")
fig = plt.figure(figsize=(4, 2))

# Set up the main Axes.
gridspec = gs.GridSpec(nrows=2, ncols=1, height_ratios=[0.04,0.46])
cax00 = fig.add_subplot(gridspec[0, 0])
ax0 = fig.add_subplot(gridspec[1, 0])
#ax0 = fig.add_subplot(2, 1, 2)
ax0.set_xlabel("Horizontal offset [km]")
ax0.set_ylabel("Depth [km]")
ax0.set_xlim(2, 8)
ax0.set_ylim(0, 3)
ax0.invert_yaxis()

x = np.arange(2, 8, dlat) #x軸の描画範囲の生成。0から10まで0.05刻み。
y = np.arange(0, 3, ddep) #y軸の描画範囲の生成。0から10まで0.05刻み。
X, Y = np.meshgrid(x, y)
Z = np.sin(X) + np.cos(Y)   # 表示する計算式の指定。等高線はZに対して作られる。
for j in range(len(x)):
    for i in range(0,m_s_dp[0]):
        Z[i][j] = m_s_km[0]*0.01
    for i in range(m_s_dp[0],m_s_dp[1]):
        Z[i][j] = 0
    for i in range(m_s_dp[1],m_s_dp[2]):
        Z[i][j] = 0
    for i in range(m_s_dp[2],sdep):
        Z[i][j] = 0
    
qmesh = ax0.pcolormesh(
    X,Y,Z,
    cmap=plt.get_cmap("seismic"),
    norm=Normalize(vmin=-0.1, vmax=0.1)
)
cbar = fig.colorbar(qmesh, cax=cax00, orientation="horizontal", ticks=[-0.1,-0.05,0,0.05,0.1])
cbar.set_label("d-Velocity [m/s]")
cbar.set_clim(-0.1,0.1)

# Plot the source location.
ax0.scatter(
    src_coords[0], src_coords[1], 
    marker="*",
    s=250,
    facecolor="w",
    edgecolor="k"
)
#print(rays[4][:,0])
#print(rays[4][:,1])
print(rays[4][:,:])
np.savetxt('out.csv',rays[4],delimiter=',')

# Plot the receiver location.
ax0.scatter(
    rec_coords[0], rec_coords[1], 
#        marker=align_marker("v", valign="bottom"),
    s=250,
    facecolor="w",
    edgecolor="k"
)
    
label = True
# Plot the synthetic raypaths.
#for decimation_factor in rays:
#    ax0.plot(
#        rays[decimation_factor][:,0], 
#        rays[decimation_factor][:,1],
#        linewidth=1,
#        label=f"d={decimation_factor}" if label is True else None
#    )
ax0.plot(
    rays[4][:,0], 
    rays[4][:,1],
    linewidth=1,
    label=f"d={decimation_factor}" if label is True else None
)

label = ord("a")
ax0.text(
    0, 1.05, f"{chr(label)})",
    ha="center",
    va="bottom",
    transform=ax0.transAxes
)
#fig.legend(loc="center right")
#fig.tight_layout()

plt.savefig('figure.png')
shutil.copyfile("figure.png", "/mnt/owncloud_webdav/webdav/figureP.png")
