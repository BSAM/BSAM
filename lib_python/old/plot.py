from matplotlib import cm
import matplotlib.pyplot as plt
import numpy
from pylab import *

# Open file to get data
f = open('./out/m00001.dat')
f.readline()
f.readline()
f.readline()
tmp = f.readline().split()
ndims = int(tmp[1])
nvars = int(tmp[3])

# Read grid data
dx   = [float(x) for x in f.readline().split()]
xmin = [float(x) for x in f.readline().split()]
xmax = [float(x) for x in f.readline().split()]
n    = [int(x)   for x in f.readline().split()]

# Read variables
f.readline()
u   = numpy.zeros([n[0]+2, n[1]+2])
phi = numpy.zeros([n[0]+2, n[1]+2])
rhs = numpy.zeros([n[0]+2, n[1]+2])
uan = numpy.zeros([n[0]+2, n[1]+2])
for i in range(n[0]+2):
    for j in range(n[1]+2):
        tmp = [float(x) for x in f.readline().split()]
        u[i][j]   = tmp[0]
        phi[i][j] = tmp[1]
        rhs[i][j] = tmp[2]
        uan[i][j] = tmp[3]

# Close file
f.close()

# Create x and y data
x = numpy.linspace(xmin[0], xmax[0], n[0]+2)
y = numpy.linspace(xmin[1], xmax[1], n[1]+2)

plot2D = False
if plot2D:
    fig = plt.figure(figsize=(16,8))
    levels = numpy.arange(-1.1,1.1,0.01)

    ax = fig.add_subplot(121)
    ax.set_title("uan")
    ax.set_aspect(1)
    cs1 = ax.contourf(x, y, uan*phi, levels, antialiased=False, cmap=cm.Greys)
    fig.colorbar(cs1, shrink=0.7)
    
    #ax = fig.add_subplot(122)
    #ax.set_title("phi")
    #ax.set_aspect(1)
    #cs4 = ax.contourf(x, y, phi, levels, antialiased=False, cmap=cm.Greys)
    #fig.colorbar(cs4, shrink=0.7)
    
    # ax = fig.add_subplot(223)
    # ax.set_title("f")
    # ax.set_aspect(1)
    # print rhs[n[0]/2-2:n[0]/2+4, n[1]/2-2:n[1]/2+4]
    # cs3 = ax.contourf(x, y, rhs, levels, antialiased=False, cmap=cm.Greys)
    # fig.colorbar(cs3, shrink=0.7)

    ax = fig.add_subplot(122)
    ax.set_title("u")
    ax.set_aspect(1)
    cs2 = ax.contourf(x, y, u*phi, levels, antialiased=False, cmap=cm.Greys)
    fig.colorbar(cs2, shrink=0.7)
else:
    fig = plt.figure()
    ax = fig.add_subplot(111)
    ax.set_title("profile of u*phi and u_an*phi")
    #ax.set_aspect(1)
    #ax.plot(x, u[:,n[1]/2]*phi[:,n[1]/2],'k--')
    ax.plot(x, u[:,n[1]/2],'k--')
    ax.plot(x, uan[:,n[1]/2],'k')
    ax.legend(('u', 'uan'), 'upper center')
    ax.set_ylim(-0.1,0.30)
    ax.set_xlim(-1.0,1.0)
    #ax.set_aspect(5.71)

plt.show()
