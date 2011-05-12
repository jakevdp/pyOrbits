"""
A set of potentials to use for integrating orbits
"""
import pylab
from pyOrbits.potentials.spherical import *

if __name__ == '__main__':
    
    halos = [ SIS(60),
              cSIS(10**8,4),
              NFW(10**8,4),
              NFW_a(10**8,4,0.5),
              NFW_a(10**8,4,0) ]
    
    x = 10**numpy.linspace(-1,3,1000)
    
    for attr in ('density','vc','mass'):
        pylab.figure()
        for halo in halos:
            y = getattr(halo,attr)(x)
            pylab.loglog(x,y,label=str(halo))
        pylab.legend(loc=0)
        pylab.xlabel('r')
        pylab.ylabel(attr)
    pylab.show()
