"""
code to integrate the infall of one halo on another.
This is useful in determining the appropriate initial 
conditions and integration time for an orbit run.
"""

import numpy
import pylab
from scipy import integrate
from pyOrbits.potentials import *

dG = 1.#/3.
G = G*(1+dG)

def fprime(u,t=None,M=1.6E11):
    """
    u is [d,d']
    d is in kpc
    d' is in km/s
    """
    return numpy.array( [u[1]*1E3*yr/kpc,
                         -numpy.sign(u[0])\
                             *(G*M/u[0]**2)*Msun*yr/kpc**2*1E-3] )

def compute_by_energy(r,v,M):
    #print 1./(r*kpc)**2
    #print (v*1E3)**2/(2*G*M*Msun)
    return 1./( 1./(r*kpc) - (v*1E3)**2/(2*G*M*Msun) )/kpc

def compute_mass(r,v):
    return 0.5 * (r*kpc) * (v*1E3)**2 / G / Msun

def main(u0 = [240,0],
         M=1.6E11,
         tmax = 2.8E9,
         Nsteps = 1000):
    t = numpy.linspace(0,tmax,Nsteps)

    x_t = integrate.odeint(fprime,u0,t,args=(M,))
    
    print x_t[-1]
    
    print compute_by_energy(u0[0],u0[1],M)
    
    pylab.plot(t,x_t[:,0])
    pylab.show()

if __name__ == '__main__':
    main()
