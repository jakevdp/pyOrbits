"""
tools for integrating orbits
"""
import numpy
from scipy import integrate

def integrate_orbits(fprime,
                     pos,vel,
                     t,
                     outfile=None,
                     extra_args=(),
                     quiet=False,
                     param_dict={}):
    """
    |integrate orbits.
    | pos:    a size [D,N] array of initial positions
    | vel:    a size [D,N] array of corresponding initial velocities
    |   D = number of positional phase space dimensions
    |   N = number of stars
    | t:      a length-M array of time steps.  pos and vel are the
    |         values at time t[0]
    | fprime: a function which accepts (y,t,*extra_args)
    |          y is the size 2*D array of positions and velocities at time t
    |          returns: a size 2*D array corresponding to dy/dt (t)
    | 
    | returns: orbits, a size [N,M,2*D] array of the computed orbits
    """
    (D,N) = pos.shape
    assert vel.shape==(D,N)
    M = len(t)
    orbits = numpy.empty((N,M,2*D))

    for i in range(N):
        if not quiet and i and (i%100==0): print "%i/%i" % (i,N)
        
        y0 = numpy.concatenate( (pos[:,i],vel[:,i]) )
        
        orbits[i,:,:] = integrate.odeint(fprime,y0,t)

    if outfile is not None:
        numpy.savez(outfile,orbits=orbits,**param_dict)

    return orbits
    
