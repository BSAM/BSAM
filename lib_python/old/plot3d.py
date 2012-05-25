from mpl_toolkits.mplot3d import Axes3D
from matplotlib import cm
from matplotlib.ticker import LinearLocator, FixedLocator, FormatStrFormatter
import matplotlib.pyplot as plt
import numpy
from pprint import pprint
from pylab import *

# 3D plots
X, Y = numpy.meshgrid(x, y)
Z = ma.array(uan)
Z[:10,:] = ma.masked
fig = plt.figure()
ax = fig.add_subplot(111,projection='3d')
surf = ax.plot_surface(X, Y, Z, rstride=1, cstride=1, cmap=cm.jet,
                       linewidth=0, antialiased=False)
fig.colorbar(surf, shrink=0.5)
ax.contourf(x,y,rhs)
ax.set_zlim3d(0.1,1)
plt.show()
