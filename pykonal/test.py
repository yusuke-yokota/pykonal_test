#!/home/sgo/anaconda3/bin/python3
from scipy.interpolate import Rbf
import numpy as np
x, y, z, d = np.random.rand(4, 50)
print(x)
print(y)
print(z)
print(d)
rbfi = Rbf(x, y, z, d)  # radial basis function interpolator instance
xi = yi = zi = np.linspace(0, 1, 100)
di = rbfi(xi, yi, zi)   # interpolated values
xi = yi = zi = np.linspace(0, 1, 1000)
dit = rbfi(xi, yi, zi)   # interpolated values
#di.shape
np.savetxt('test_x.csv',x,delimiter=',')
np.savetxt('test_y.csv',y,delimiter=',')
np.savetxt('test_z.csv',z,delimiter=',')
np.savetxt('test_dit.csv',dit,delimiter=',')
np.savetxt('test_di.csv',di,delimiter=',')
#print(di)
