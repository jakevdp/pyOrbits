import pylab
import numpy
from pyOrbits.potentials.spherical import *

def plot_from_potential(halos = [ NFW(1E7, 4.0),
                                  cSIS(1E7, 4.0),
                                  cSIS(1E7, 1.0)  ]):
    Npts = 50
    r = 10**numpy.linspace(-2,numpy.log10(10),Npts)
    
    for halo in halos:
        pylab.plot(r,halo.vc(r),label=str(halo) )
    pylab.legend(loc=2)
    pylab.xlabel(r'$\mathdefault{r\ (kpc)}$')
    pylab.ylabel(r'$\mathdefault{v_c\ (km/s)}$')
    pylab.ylim(0,80)

def view_2D(xyz,view):
    """
    xyz is shape (3,N)
    view is (alt, az) from positive x-axis
    returns the x,y components of the 3D collection of points
    as viewed from (alt,az)
    """
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
    
    return x,y,z

def plot_from_angle(filename,
                    view=(14,-73),
                    index=-1,
                    indicator=True):
    """
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

    #indicator point
    if indicator:
        xyz[:,-1] = [4,0,0] 

    x,y,z = view_2D(xyz,view)

    c = numpy.zeros(xyz.shape[1])
    s = 9*numpy.ones(len(c))
    if indicator:
        c[-1]=1
        s[-1]*=4

    pylab.scatter(y,z,c=c,s=s)

def scatter_rotation(filename,view=(14,-73),
                     index=-1):
    """
    display a color-coded line-of-sight velocity plot
    """
    X = numpy.load(filename)
    orbits = X['orbits']

    pos = orbits[:,index,:3].T
    vel = orbits[:,index,3:6].T

    x,y,z = view_2D(pos,view)
    vx,vy,vz = view_2D(vel,view)

    pylab.scatter(y,z,c=vx,lw=0)

def binned_plot(x,y,bins=40):
    if not hasattr(bins,'__len__'):
        bins = numpy.linspace(min(x),max(x)+1E-10,bins+1)
        
    x = numpy.asarray(x)
    y = numpy.asarray(y)

    ym = numpy.zeros(len(bins)-1)
    ystd  = numpy.zeros(len(bins)-1)

    for i in range(len(bins)-1):
        j = numpy.where( (x>=bins[i])&(x<bins[i+1]) )
        yj = y[j]
        ym[i] = numpy.mean(yj)
        ystd[i] = numpy.std(yj)
    
    pylab.errorbar(0.5*(bins[:-1]+bins[1:]),
                   ym,ystd,fmt='-k')

def plot_rotation_curve(filename,
                        view,
                        index=-1):
    X = numpy.load(filename)
    orbits = X['orbits']

    pos = orbits[:,index,:3].T
    vel = orbits[:,index,3:6].T

    x,y,z = view_2D(pos,view)
    vx,vy,vz = view_2D(vel,view)
    
    #We need to find the y-z orientation with the largest variance.
    # use PCA!!
    M = numpy.concatenate([[y],[z]])

    M_m = M.mean(1)
    M -= M_m[:,None]
    U,s,V = numpy.linalg.svd(M)
    
    p1,p2 = numpy.dot(U.T,M)

    if abs(U[0,1])>abs(U[0,0]):
        zl = numpy.linspace(min(z),max(z),10)
        yl = M_m[0] + (zl-M_m[1])*U[0,0]/U[0,1]
    else:
        yl = numpy.linspace(min(y),max(y),10)
        zl = M_m[1] + (yl-M_m[0])*U[0,1]/U[0,0]
    
    pylab.figure(figsize=(8,8))
    pylab.axes( (0.1,0.07,0.8,0.23) )
    binned_plot(p1,vx)
    pylab.xlabel('distance along principal axis')
    pylab.ylabel(r'$\mathdefault{v_{los}}$')

    pylab.axes( (0.1,0.35,0.8,0.6) )
    pylab.scatter(y,z,c=vx,lw=0)
    pylab.plot(yl,zl,'-k',lw=3)
    pylab.colorbar().set_label(r'$\mathdefault{v_{los} (km/s)}$')

    #set axes limits so x and y scales are equal
    xlim = pylab.xlim()
    ylim = pylab.ylim()
    
    x0 = numpy.mean(xlim)
    y0 = numpy.mean(ylim)
    xr = xlim[1]-xlim[0]
    yr = ylim[1]-ylim[0]
    
    d = 0.5*max(xr,yr)
    
    pylab.xlim(x0-d,x0+d)
    pylab.ylim(y0-d,y0+d)


if __name__ == '__main__':
    #plot_from_potential()
    #pylab.savefig('fig/rotation_curve.eps')
    #pylab.show()

    #filename = 'orbits/NFW_c10_rho9E6.npz'
    filename = 'orbits/orbit.npz'
    az = -45#-85
    for index in (0,-1):
        for alt in (0,):
            plot_rotation_curve(filename,(alt,az),index)
            #pylab.xlim(-4.5,4.5)
            #pylab.ylim(-4.5,4.5)
        
        
    pylab.show()
