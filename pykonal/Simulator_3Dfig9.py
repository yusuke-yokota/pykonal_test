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
##
#################################
##Input parameters start #####
# GNSS noise ######
ghptb = icfg.get("HyperParameters", "Ghptb")
gvptb = icfg.get("HyperParameters", "Gvptb")
knot  = icfg.get("HyperParameters", "Knot")
lknot  = float(knot)*1.852
lhptb = icfg.get("HyperParameters", "Leptb")
lnptb = icfg.get("HyperParameters", "Lnptb")
luptb = icfg.get("HyperParameters", "Luptb")
ghptb,gvptb,lhptb,lnptb,luptb = float(ghptb),float(gvptb),float(lhptb),float(lnptb),float(luptb)
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
beta10 = 1.0
beta15 = 1.5
beta20 = 2.0
#################################
ltime  = 43200 # int(3600.*kyori/lknot)
gpe = cn.powerlaw_psd_gaussian(beta15, ltime)
gpe = ghptb * (gpe - np.average(gpe))
gpn = cn.powerlaw_psd_gaussian(beta15, ltime)
gpn = ghptb * (gpn - np.average(gpn))
gpu = cn.powerlaw_psd_gaussian(beta15, ltime)
gpu = gvptb * (gpu - np.average(gpu))
#################################
ld     = pd.read_csv('linesample.csv')
da     = [[[] for j in range(13)] for i in range(len(ld['se']))]
for n in range(len(ld['se'])):
  lstart = np.array([int(ld['se'][n]*sdep),int(ld['sn'][n]*sdep),0])
  lend   = np.array([int(ld['ee'][n]*sdep),int(ld['en'][n]*sdep),0])
  line   = (lend-lstart)/sdep
  kyori  = np.linalg.norm(line)
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
  sewn  = np.linspace(0, 10, ld['shot'][n])
  snsn  = np.linspace(0, 10, ld['shot'][n])
  sudn  = np.linspace(0, 10, ld['shot'][n])
  for i in range(len(sewr)):
    sewn[i] = sewr[i]+gpe[int(stim[i])]
    snsn[i] = snsr[i]+gpn[int(stim[i])]
    sudn[i] = sudr[i]+gpu[int(stim[i])]
  rtim = stim
  rewr = sewr
  rnsr = snsr
  rudr = sudr
  rewn = sewn
  rnsn = snsn
  rudn = sudn
  rdi  = sdi
  din1  = np.linspace(0, 10, ld['shot'][n])
  din2  = np.linspace(0, 10, ld['shot'][n])
  din3  = np.linspace(0, 10, ld['shot'][n])
  din4  = np.linspace(0, 10, ld['shot'][n])
  for i in range(ld['shot'][n]):
    nin = int(lstart[0]+i*(lend[0]-lstart[0])/ld['shot'][n])
    nie = int(lstart[1]+i*(lend[1]-lstart[1])/ld['shot'][n])
    din1[i]  = (np.sqrt((sewn[i]-solverc1.velocity.nodes[slon1,sdep1,slat1,0])**2 + (snsn[i]-solverc1.velocity.nodes[slon1,sdep1,slat1,2])**2 + (sudn[i]-solverc1.velocity.nodes[slon1,sdep1,slat1,1])**2) - np.sqrt((sewr[i]-solverc1.velocity.nodes[slon1,sdep1,slat1,0])**2 + (snsr[i]-solverc1.velocity.nodes[slon1,sdep1,slat1,2])**2 + (sudr[i]-solverc1.velocity.nodes[slon1,sdep1,slat1,1])**2)) / solverc1.velocity.values[nin,0,nie]
    din2[i]  = (np.sqrt((sewn[i]-solverc2.velocity.nodes[slon2,sdep2,slat2,0])**2 + (snsn[i]-solverc2.velocity.nodes[slon2,sdep2,slat2,2])**2 + (sudn[i]-solverc2.velocity.nodes[slon2,sdep2,slat2,1])**2) - np.sqrt((sewr[i]-solverc2.velocity.nodes[slon2,sdep2,slat2,0])**2 + (snsr[i]-solverc2.velocity.nodes[slon2,sdep2,slat2,2])**2 + (sudr[i]-solverc2.velocity.nodes[slon2,sdep2,slat2,1])**2)) / solverc2.velocity.values[nin,0,nie]
    din3[i]  = (np.sqrt((sewn[i]-solverc3.velocity.nodes[slon3,sdep3,slat3,0])**2 + (snsn[i]-solverc3.velocity.nodes[slon3,sdep3,slat3,2])**2 + (sudn[i]-solverc3.velocity.nodes[slon3,sdep3,slat3,1])**2) - np.sqrt((sewr[i]-solverc3.velocity.nodes[slon3,sdep3,slat3,0])**2 + (snsr[i]-solverc3.velocity.nodes[slon3,sdep3,slat3,2])**2 + (sudr[i]-solverc3.velocity.nodes[slon3,sdep3,slat3,1])**2)) / solverc3.velocity.values[nin,0,nie]
    din4[i]  = (np.sqrt((sewn[i]-solverc4.velocity.nodes[slon4,sdep4,slat4,0])**2 + (snsn[i]-solverc4.velocity.nodes[slon4,sdep4,slat4,2])**2 + (sudn[i]-solverc4.velocity.nodes[slon4,sdep4,slat4,1])**2) - np.sqrt((sewr[i]-solverc4.velocity.nodes[slon4,sdep4,slat4,0])**2 + (snsr[i]-solverc4.velocity.nodes[slon4,sdep4,slat4,2])**2 + (sudr[i]-solverc4.velocity.nodes[slon4,sdep4,slat4,1])**2)) / solverc4.velocity.values[nin,0,nie]
  da[n][0] = np.full(len(sewr),"S%s"%(str(ld['set'][n]).zfill(2)))
  da[n][1] = np.full(len(sewr),"L%s"%(str(n+1).zfill(2)))
  da[n][2] = np.full(len(sewr),"M01")
  da[n][3] = sdi + rdi
  for i in range(0,len(sewr1[0::4])):
#    da[n][2][i] = "M01"
    da[n][3][i] = da[n][3][i] + 2.*din1[i]
  for i in range(len(sewr1[0::4]),len(sewr1[0::4])+len(sewr1[1::4])):
    da[n][2][i] = "M02"
    da[n][3][i] = da[n][3][i] + 2.*din2[i]
  for i in range(len(sewr1[0::4])+len(sewr1[1::4]),len(sewr1[0::4])+len(sewr1[1::4])+len(sewr1[2::4])):
    da[n][2][i] = "M03"
    da[n][3][i] = da[n][3][i] + 2.*din3[i]
  for i in range(len(sewr1[0::4])+len(sewr1[1::4])+len(sewr1[2::4]),len(sewr)):
    da[n][2][i] = "M04"
    da[n][3][i] = da[n][3][i] + 2.*din4[i]
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
