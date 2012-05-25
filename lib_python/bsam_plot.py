from matplotlib import cm
from matplotlib.patches import Rectangle
import numpy

def create_levels(xmin,xmax,n):
    h = (xmax-xmin)/n
    return [xmin + h*i for i in xrange(n+1)]

def plot_limits(ax, limits):
    ax.set_aspect((limits[1][1]-limits[1][0])/(limits[0][1]-limits[0][0]))
    ax.set_ylim(limits[0])
    ax.set_xlim(limits[1])

def plot_interface(ax, patch, ivar=1):
    nx = patch['nx']
    ny = patch['ny']
    x  = numpy.linspace(patch['x'][0], patch['x'][1], nx+2)
    y  = numpy.linspace(patch['y'][0], patch['y'][1], ny+2)
    z  = patch['data'][ivar]
    return ax.contour(x, y, z, [0.5], colors='r', antialiased=False)

def plot_contourf(ax, patch, ivar, levels=False):
    nx = patch['nx']
    ny = patch['ny']
    x  = numpy.linspace(patch['x'][0], patch['x'][1], nx+2)
    y  = numpy.linspace(patch['y'][0], patch['y'][1], ny+2)
    z  = patch['data'][ivar]

    # Create figure
    ax.set_aspect(1)
    if levels:
        ex = 'neither'
        if min(levels) > z.min(): ex = 'min'
        if max(levels) < z.max():
            if ex == 'min':
                ex = 'both'
            else:
                ex = 'max'
        return ax.contourf(x, y, z, levels, extend=ex,
                          antialiased=False, cmap=cm.Greys)
    else:
        return ax.contourf(x, y, z, antialiased=False, cmap=cm.Greys)

def plot_bounding_boxes(ax,data):
    cl = create_levels(0.9,0.1,data['maxlevel'])
    for level in data['levels'][1:]:
        color = str(cl.pop())
        for patch in level:
            x = patch['x']
            y = patch['y']
            ax.add_patch(Rectangle([x[0], y[0]], x[1] - x[0], y[1] - y[0],
                                   ec=color, fc="none"))

def plot_grid(ax,data,maxlevel=3):
    width = 0.1
    for level in data['levels'][:maxlevel]:
        for patch in level:
            x0 = patch['x'][0]
            x1 = patch['x'][1]
            y0 = patch['y'][0]
            y1 = patch['y'][1]
            h  = patch['h']
            xgrid = numpy.mgrid[x0:x1:h]
            ygrid = numpy.mgrid[y0:y1:h]
            ax.vlines(xgrid, y0, y1, lw=width)
            ax.hlines(ygrid, x0, x1, lw=width)

def plot_profile(ax, patch, plot_list, limits=False, direction='x'):
    if direction == 'x':
        n = patch['nx']
        x = numpy.linspace(patch['x'][0], patch['x'][1], n+2)
        a = patch['data'][:,patch['ny']/2,:]
    else:
        n = patch['ny']
        x = numpy.linspace(patch['y'][0], patch['y'][1], n+2)
        a = patch['data'][:,:,patch['nx']/2]

    # Add plots
    for var, style in plot_list:
        if type(var) == tuple:
            ax.plot(x, a[var[0]]*a[var[1]], style)
        else:
            ax.plot(x, a[var], style)

    # Set limits
    if limits:
        ax.set_ylim(limits[0])
        ax.set_xlim(limits[1])
        #ax.set_aspect((limits[1][1]-limits[1][0])/(limits[0][1]-limits[0][0]))

def plot_grid_1d(ax,data,maxlevel=10,coor='x'):
    if coor == 'x':
        y0, y1 = ax.get_ylim()
    else:
        y0, y1 = ax.get_xlim()
    cl = create_levels(0.9,0.0,data['maxlevel'])
    for level in data['levels'][:maxlevel]:
        color = str(cl.pop())
        for patch in level:
            x0 = patch[coor][0]
            x1 = patch[coor][1]
            h  = patch['h']
            grid = numpy.mgrid[x0:x1:h]
            ax.vlines(grid, y0, y1, colors=color, lw=0.2)
