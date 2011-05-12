"""
routines for plotting orbits
"""

import numpy
import pylab
from mpl_toolkits.mplot3d import Axes3D

def plot_orbits_3D(filename='orbits_0.0.npz',
                   N=1000,
                   final_index=-1):
    """
    3D plot of the orbits in the file
    plot the final position of the first N orbits in given file
    """
    X = numpy.load(filename)

    orbits = X['orbits']
    tmax   = X['tmax']

    t = numpy.linspace(0,tmax,orbits.shape[1]-1)
    
    orbits = orbits[:N,:,:]

    x0 = orbits[:,0,0]
    y0 = orbits[:,0,1]
    z0 = orbits[:,0,2]

    xf = orbits[:,final_index,0]
    yf = orbits[:,final_index,1]
    zf = orbits[:,final_index,2]

    fig1 = pylab.figure()
    fig2 = pylab.figure()
    
    ax1 = Axes3D(fig1)
    ax2 = Axes3D(fig2)

    #ax1.scatter3D(x0,y0,z0)
    #ax2.scatter3D(xf,yf,zf)

    ax1.plot3D(x0,y0,z0,'.k',ms=2)
    ax2.plot3D(xf,yf,zf,'.k',ms=2)

    for ax in (ax1,ax2,):
        ax.view_init(14,-73)

    for ax in (ax1,ax2,):
        xlim = ax.get_xlim3d()
        ylim = ax.get_ylim3d()
        zlim = ax.get_xlim3d()

        axmax = max(abs(numpy.concatenate((xlim,ylim,zlim))))
        ax.set_xlim3d(-axmax,axmax)
        ax.set_ylim3d(-axmax,axmax)
        ax.set_zlim3d(-axmax,axmax)

        ax.set_xlabel('x (kpc)')
        ax.set_ylabel('y (kpc)')
        ax.set_zlabel('z (kpc)')

    print "disk center of mass: (R,z) = (%.2g kpc, %.2g kpc)" \
        % (numpy.sqrt(numpy.mean(xf)**2+numpy.mean(yf)**2),numpy.mean(zf))

    return fig1,fig2

def plot_from_angle(filename,
                    view=(14,-73),
                    index=-1):
    """
    2D scatter plot of the projected view of orbits
    view gives the (altitude,azimuth) measured from x=+infinity
    """
    #view=(0,0) is just the data projected onto the y-z plane.
    # from this, we rotate about the z-axis by view[1] degrees
    # then rotate about the new y-axis by view[0] degrees.
    
    #first get the x,y,z data
    X = numpy.load(filename)
    orbits = X['orbits']

    xyz = orbits[:,index,:3].T
    print xyz.shape

    xyz[:,-1] = [4,0,0]

    view = numpy.asarray(view)
    view[1]*=-1 #rotate axes, not coordinates
    
    c = numpy.cos(view*numpy.pi/180.)
    s = numpy.sin(view*numpy.pi/180.)

    Rz = numpy.array( [ [c[1], -s[1], 0],
                        [s[1], c[1] , 0],
                        [0   , 0    , 1] ] )

    Ry = numpy.array( [ [c[0] , 0, s[0]],
                        [0    , 1, 0   ],
                        [-s[0], 0, c[0]] ] )

    x,y,z = numpy.dot(Ry, numpy.dot(Rz,xyz))

    c = numpy.zeros(xyz.shape[1])
    c[-1]=1

    pylab.scatter(y,z,c=c)
    #pylab.plot(y,z,'.k')
    

def plot_infall(filename):
    """
    plot the infall with time
    """
    X = numpy.load(filename)
    orbits = X['orbits']
    tmax   = X['tmax']

    t = numpy.linspace(0,tmax,orbits.shape[1])
    
    pylab.plot(t,orbits[0,:,6],label='x')
    pylab.plot(t,orbits[0,:,7],label='y')
    pylab.plot(t,orbits[0,:,8],label='z')
    
    pylab.legend(loc=0)

def plot_final_shape(filename,final_index=-1):
    """
    plot the final shape of the disk (works only for edge-on disks
    """
    X = numpy.load(filename)
    orbits = X['orbits']
    xf = orbits[:,final_index,0]
    yf = orbits[:,final_index,1]
    zf = orbits[:,final_index,2]

    rf = numpy.sqrt(xf**2+yf**2)
    pylab.plot(rf,zf,'.k')
    pylab.xlabel('r')
    pylab.ylabel('z')

def plot_init_final(filename):
    """
    plot the initial and final positions together
    """
    X = numpy.load(filename)

    orbits = X['orbits']

    x0 = orbits[:,0,0]
    y0 = orbits[:,0,1]
    z0 = orbits[:,0,2]

    xf = orbits[:,-1,0]
    yf = orbits[:,-1,1]
    zf = orbits[:,-1,2]

    x = numpy.concatenate((x0,xf))
    y = numpy.concatenate((y0,yf))
    z = numpy.concatenate((z0,zf))
    
    fig = pylab.figure()
    ax = Axes3D(fig)
    ax.plot3D(x,y,z,'.k',ms=2)
    #ax.scatter3D(x,y,z)
    ax.set_zlim3d(-1,3)
    ax.view_init(14,-73)

    return fig
