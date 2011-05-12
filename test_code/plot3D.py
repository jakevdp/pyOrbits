from mpl_toolkits.mplot3d import Axes3D
import pylab
import numpy


def surface():
    fig = pylab.figure()
    ax = Axes3D(fig)
    
    X = numpy.linspace(0,10,50)
    Y = numpy.linspace(0,10,50)
    X,Y = numpy.meshgrid(X,Y)
    Z = numpy.sin(X)*numpy.sin(Y)*numpy.exp(-0.1*((X-5)**2+(Y-5)**2))
    
    ax.plot_wireframe(X, Y, Z, rstride=1, cstride=1, cmap=pylab.cm.jet)

def line():
    t = numpy.linspace(0,10,1000)

    x = numpy.sin(3*t)
    y = numpy.cos(3*t)
    z = numpy.cos(0.1*t*numpy.pi)

    fig = pylab.figure()
    ax = Axes3D(fig)

    ax.plot3D(x,y,z,zorder=10)

    X = numpy.linspace(-1,1,50)
    Y = numpy.linspace(-1,1,50)
    X,Y = numpy.meshgrid(X,Y)
    Z = -1+0.01*numpy.exp(-0.1*X**2-0.2*Y**2)

    ax.contourf(X,Y,Z,50,zorder=1)

def scatter():
    N = 1000
    theta = 2*numpy.pi*numpy.random.random(N)
    r = numpy.random.random(size=N)
    z = 0.2*(numpy.random.random(N)-0.5)/(r**2+1)

    x = r*numpy.cos(theta)
    y = r*numpy.sin(theta)

    fig = pylab.figure()
    ax = Axes3D(fig)

    ax.plot3D(x,y,z,'.k',ms=2)
    #ax.scatter3D(x,y,z)
    ax.set_xlim3d(-1,1)
    ax.set_ylim3d(-1,1)
    ax.set_zlim3d(-2,2)

    

if __name__ == '__main__':
    scatter()
    pylab.show()
