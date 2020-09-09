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
##
sv = pd.read_csv('svpC.csv')
##
# Initialize the solver.######
##Input parameters start #####
m_s_dp=[  10,  20,  30]
m_s_km=[  0.0, 0.0,  0.0] # cm/s/km
nlon,ndep,nlat=201,101,201#8000,8000,4000
dlon,ddep,dlat=0.010,0.010,0.010
slon,sdep,slat=100,100,100
#nlon,ndep,nlat=400,201,400#8000,8000,4000
#dlon,ddep,dlat=0.005,0.005,0.005
#slon,sdep,slat=200,200,200
##Input parameters  end  #####
##############################
slon1,sdep1,slat1=slon,    sdep,int(slat*1.5)
slon2,sdep2,slat2=int(slon*1.5),sdep,slat
slon3,sdep3,slat3=slon,    sdep,int(slat*0.5)
slon4,sdep4,slat4=int(slon*0.5),sdep,slat
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
for s in range(nlon):
  svc = [[] for i in range(ndep)]
  for i in range(0,m_s_dp[0]):
    svc[i] = np.array([m_s_km[0]*0.001*0.1*(1./float(2*slon))*float(s-slon)] * nlat)
  for i in range(m_s_dp[0],m_s_dp[1]):
    svc[i] = np.array([m_s_km[1]*0.001*0.1*(1./float(2*slon))*float(s-slon)] * nlat)
  for i in range(m_s_dp[1],ndep):
    svc[i] = np.array([m_s_km[2]*0.001*0.1*(1./float(2*slon))*float(s-slon)] * nlat)
  svp[s] = [svc[0]]
  for i in range(1,ndep):
    svp[s] = np.append(svp[s], [svc[i]], axis=0)
a1 = [svp[0]]
for s in range(1,nlon):
  a1 = np.append(a1, [svp[s]], axis=0)
####################
#solverb1.velocity.values = a0
solverb1.velocity.values = o0
src_idx = slon1, sdep1, slat1
solverb1.traveltime.values[src_idx] = 0
solverb1.unknown[src_idx] = False
solverb1.trial.push(*src_idx)
solverb1.solve()
solverc1.velocity.values = a0 + a1
solverc1.traveltime.values[src_idx] = 0
solverc1.unknown[src_idx] = False
solverc1.trial.push(*src_idx)
solverc1.solve()
#solverb2.velocity.values = a0
solverb2.velocity.values = o0
src_idx = slon2, sdep2, slat2
solverb2.traveltime.values[src_idx] = 0
solverb2.unknown[src_idx] = False
solverb2.trial.push(*src_idx)
solverb2.solve()
solverc2.velocity.values = a0 + a1
solverc2.traveltime.values[src_idx] = 0
solverc2.unknown[src_idx] = False
solverc2.trial.push(*src_idx)
solverc2.solve()
#solverb3.velocity.values = a0
solverb3.velocity.values = o0
src_idx = slon3, sdep3, slat3
solverb3.traveltime.values[src_idx] = 0
solverb3.unknown[src_idx] = False
solverb3.trial.push(*src_idx)
solverb3.solve()
solverc3.velocity.values = a0 + a1
solverc3.traveltime.values[src_idx] = 0
solverc3.unknown[src_idx] = False
solverc3.trial.push(*src_idx)
solverc3.solve()
#solverb4.velocity.values = a0
solverb4.velocity.values = o0
src_idx = slon4, sdep4, slat4
solverb4.traveltime.values[src_idx] = 0
solverb4.unknown[src_idx] = False
solverb4.trial.push(*src_idx)
solverb4.solve()
solverc4.velocity.values = a0 + a1
solverc4.traveltime.values[src_idx] = 0
solverc4.unknown[src_idx] = False
solverc4.trial.push(*src_idx)
solverc4.solve()
#####################
#####################
ttb = solverb1.traveltime.values
ttc = solverc1.traveltime.values
ttm = 1000.*np.sqrt((solverc1.velocity.nodes[:,:,:,0]-solverc1.velocity.nodes[slon1,sdep1,slat1,0])**2+(solverc1.velocity.nodes[:,:,:,1]-solverc1.velocity.nodes[slon1,sdep1,slat1,1])**2+(solverc1.velocity.nodes[:,:,:,2]-solverc1.velocity.nodes[slon1,sdep1,slat1,2])**2)/tsv
ttd1 = ttc - ttb + ttm
#####################
ttb = solverb2.traveltime.values
ttc = solverc2.traveltime.values
ttm = 1000.*np.sqrt((solverc2.velocity.nodes[:,:,:,0]-solverc2.velocity.nodes[slon2,sdep2,slat2,0])**2+(solverc2.velocity.nodes[:,:,:,1]-solverc2.velocity.nodes[slon2,sdep2,slat2,1])**2+(solverc2.velocity.nodes[:,:,:,2]-solverc2.velocity.nodes[slon2,sdep2,slat2,2])**2)/tsv
ttd2 = ttc - ttb + ttm
#####################
ttb = solverb3.traveltime.values
ttc = solverc3.traveltime.values
ttm = 1000.*np.sqrt((solverc3.velocity.nodes[:,:,:,0]-solverc3.velocity.nodes[slon3,sdep3,slat3,0])**2+(solverc3.velocity.nodes[:,:,:,1]-solverc3.velocity.nodes[slon3,sdep3,slat3,1])**2+(solverc3.velocity.nodes[:,:,:,2]-solverc3.velocity.nodes[slon3,sdep3,slat3,2])**2)/tsv
ttd3 = ttc - ttb + ttm
#####################
ttb = solverb4.traveltime.values
ttc = solverc4.traveltime.values
ttm = 1000.*np.sqrt((solverc4.velocity.nodes[:,:,:,0]-solverc4.velocity.nodes[slon4,sdep4,slat4,0])**2+(solverc4.velocity.nodes[:,:,:,1]-solverc4.velocity.nodes[slon4,sdep4,slat4,1])**2+(solverc4.velocity.nodes[:,:,:,2]-solverc4.velocity.nodes[slon4,sdep4,slat4,2])**2)/tsv
ttd4 = ttc - ttb + ttm
#####################

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
qmesh = ax10.pcolormesh(
    solverc1.velocity.nodes[:,:,0,0],
    solverc1.velocity.nodes[:,:,0,1],
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
    solverc1.traveltime.nodes[:,0,:,0],
    solverc1.traveltime.nodes[:,0,:,2],
    ttd1[:,0,:],
    cmap=plt.get_cmap("coolwarm"),
    norm=Normalize(vmin=-0.0002,vmax=0.0002)
)
ax11.contour(
    solverc1.traveltime.nodes[:,0,:,0],
    solverc1.traveltime.nodes[:,0,:,2],
    ttd1[:,0,:],
    colors="k",
    levels=np.arange(ttd1.min(), ttd1.max(), 0.1),
    linewidths=0.2,
    linestyles="--"
)
#np.savetxt('centerline.csv',tt0[:,0,slon],delimiter=',')
#print(ttd[:,0,slon])
ax11.plot([1, 1], [0, 2], 'w-', lw=2)
ax11.plot([0, 2], [1, 1], 'w-', lw=2)
ax11.plot([1-1/math.sqrt(2), 1+1/math.sqrt(2)], [1-1/math.sqrt(2), 1+1/math.sqrt(2)], 'w-', lw=2)
ax11.plot([1-1/math.sqrt(2), 1+1/math.sqrt(2)], [1+1/math.sqrt(2), 1-1/math.sqrt(2)], 'w-', lw=2)
ax11.plot([1, 2, 1, 0, 1], [0, 1, 2, 1, 0], 'w-', lw=2)
ax11.plot([1, 1.5, 1, 0.5, 1], [0.5, 1, 1.5, 1, 0.5], 'w-', lw=2)
ax11.plot(1, 1.0, marker='*',markersize=4)
ax11.set_xlim([0,2])
ax11.set_ylim([0,2])
cbar = fig.colorbar(qmesh, cax=cax01, orientation="horizontal", ticks=[-1.,0.0,1.])
cbar.set_label("")
cbar.set_clim(-1.,1.)
#cbar.ax.xaxis.tick_top()
#cbar.ax.xaxis.set_label_position("top")
plt.savefig('figure.png')
#################################
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
da     = [[[] for j in range(13)] for i in range(len(ld['se']))]
for n in range(len(ld['se'])):
  lstart = np.array([int(ld['se'][n]*sdep),int(ld['sn'][n]*sdep),0])
  lend   = np.array([int(ld['ee'][n]*sdep),int(ld['en'][n]*sdep),0])
  line   = (lend-lstart)/sdep
  kyori  = np.linalg.norm(line)
#  ltime  = int(3600.*kyori/lknot)+30
#  gpe = cn.powerlaw_psd_gaussian(beta15, ltime)
#  gpe = gpe - np.average(gpe)
#  gpn = cn.powerlaw_psd_gaussian(beta15, ltime)
#  gpn = gpn - np.average(gpn)
#  gpu = cn.powerlaw_psd_gaussian(beta15, ltime)
#  gpu = gpu - np.average(gpu)
#  l_e = cn.powerlaw_psd_gaussian(beta20, ltime)
#  l_e = l_e - np.average(l_e)
#  l_n = cn.powerlaw_psd_gaussian(beta20, ltime)
#  l_n = l_n - np.average(l_n)
#  l_z = cn.powerlaw_psd_gaussian(beta10, ltime)
#  l_z = l_z - np.average(l_z)
#################################
  sewr1 = np.linspace(0, 10, ld['shot'][n])
  snsr1 = np.linspace(0, 10, ld['shot'][n])
  sudr1 = np.linspace(0, 10, ld['shot'][n])
  sdi1  = np.linspace(0, 10, ld['shot'][n])
  sdi2  = np.linspace(0, 10, ld['shot'][n])
  sdi3  = np.linspace(0, 10, ld['shot'][n])
  sdi4  = np.linspace(0, 10, ld['shot'][n])
  for i in range(ld['shot'][n]):
    nin = int(lstart[0]+i*(lend[0]-lstart[0])/ld['shot'][n])
    nie = int(lstart[1]+i*(lend[1]-lstart[1])/ld['shot'][n])
    sewr1[i] = solverc1.velocity.nodes[nin,0,nie,0]-solverc1.velocity.nodes[slon,0,slat,0]
    snsr1[i] = solverc1.velocity.nodes[nin,0,nie,2]-solverc1.velocity.nodes[slon,0,slat,2]
    sudr1[i] = solverc1.velocity.nodes[nin,0,nie,1]-solverc1.velocity.nodes[slon,0,slat,1]
    sdi1[i] = ttd1[nin,0,nie]
    sdi2[i] = ttd2[nin,0,nie]
    sdi3[i] = ttd3[nin,0,nie]
    sdi4[i] = ttd4[nin,0,nie]
  sewr = np.concatenate([sewr1[0::4], sewr1[1::4]], 0)
  sewr = np.concatenate([sewr, sewr1[2::4]], 0)
  sewr = np.concatenate([sewr, sewr1[3::4]], 0)
  snsr = np.concatenate([snsr1[0::4], snsr1[1::4]], 0)
  snsr = np.concatenate([snsr, snsr1[2::4]], 0)
  snsr = np.concatenate([snsr, snsr1[3::4]], 0)
  sudr = np.concatenate([sudr1[0::4], sudr1[1::4]], 0)
  sudr = np.concatenate([sudr, sudr1[2::4]], 0)
  sudr = np.concatenate([sudr, sudr1[3::4]], 0)
  stime = np.linspace(n*700., n*700.+3600.*kyori/lknot, len(sewr))
  stim = np.concatenate([stime[0::4], stime[1::4]], 0)
  stim = np.concatenate([stim, stime[2::4]], 0)
  stim = np.concatenate([stim, stime[3::4]], 0)
  sdi  = np.concatenate([sdi1[0::4], sdi2[1::4]], 0)
  sdi  = np.concatenate([sdi, sdi3[2::4]], 0)
  sdi  = np.concatenate([sdi, sdi4[3::4]], 0)
  sudr = sudr1
  sewn = sewr+0.
  snsn = snsr+0.
  sudn = sudr+0.
  rtim = stim
  rewr = sewr
  rnsr = snsr
  rudr = sudr
  rewn = rewr+0.
  rnsn = rnsr+0.
  rudn = rudr+0.
  rdi  = sdi
  da[n][0] = np.full(len(sewr),"S%s"%(str(ld['set'][n]).zfill(2)))
  da[n][1] = np.full(len(sewr),"L%s"%(str(n+1).zfill(2)))
  da[n][2] = np.full(len(sewr),"M01")
  for i in range(len(sewr1[0::4]),len(sewr1[0::4])+len(sewr1[1::4])):
    da[n][2][i] = "M02"
  for i in range(len(sewr1[0::4])+len(sewr1[1::4]),len(sewr1[0::4])+len(sewr1[1::4])+len(sewr1[2::4])):
    da[n][2][i] = "M03"
  for i in range(len(sewr1[0::4])+len(sewr1[1::4])+len(sewr1[2::4]),len(sewr)):
    da[n][2][i] = "M04"
  da[n][3] = sdi + rdi
  da[n][4] = stim
  da[n][5] = sewr
  da[n][6] = snsr
  da[n][7] = sudr
  da[n][8] = rtim
  da[n][9] = rewr
  da[n][10] = rnsr
  da[n][11] = rudr
  da[n][12] = np.full(len(sewr),"False")
for n in range(1,len(ld['se'])):
  for i in range(13):
    da[0][i] = np.concatenate([da[0][i], da[n][i]], 0)
#################################
df = pd.DataFrame({'SET': da[0][0], 'LN': da[0][1],'MT': da[0][2],'TT': da[0][3],'ResiTT': np.nan,'TakeOff': np.nan,'ST': da[0][4],'ant_e0': da[0][5]*1000.,'ant_n0': da[0][6]*1000.,'ant_u0': da[0][7]*1000.,'head0': 0.0,'pitch0': 0.0,'roll0': 0.0,'dSV0': 0.0,'RT': da[0][8],'ant_e1': da[0][9]*1000.,'ant_n1': da[0][10]*1000.,'ant_u1': da[0][11]*1000.,'head1': 0.0,'pitch1': 0.0,'roll1': 0.0,'dSV1': 0.0,'flag': da[0][12]})
df = df.sort_values('ST')
df = df.reset_index()
df = df.drop('index',axis=1)
df = df.round({'TT': 6, 'ST': 6, 'ant_e0': 5, 'ant_n0': 5, 'ant_u0': 5, 'RT': 6, 'ant_e1': 5, 'ant_n1': 5, 'ant_u1': 5})
of='~/sgobs/garpos/simulator/obsdata/IMAG/IMAG.1911.kaiyo_k4-obs.csv'
cfgfile='./initcfg/IMAG/IMAG.1901.kaiyo_k4-initcfg.ini'
hd = "# cfgfile = %s\n" % cfgfile
df.to_csv(of)
os.system('sed -i -e "1i ' + hd + '" ' + of )
