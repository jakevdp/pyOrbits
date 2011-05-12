"""
classes which combine potentials to create fprime functions
which can be put into integrators.
"""

import numpy
from .potentials.spherical import *

class fprime_halos:
    def __init__(self,*args):
        """
        arguments are pairs of spherical_potential objects, along with
        cartesian distances x = (x,y,z) in kpc
        """
        for arg in args:
            assert hasattr(arg,'__len__')
            assert len(arg)==2
            assert hasattr(arg[0],'accel')
            assert hasattr(arg[1],'__len__')
            assert len(arg[1])==3
            
        self.halos = [arg[0] for arg in args]
        self.locs = numpy.asarray( [arg[1] for arg in args] )
        self.N = len(args)
        
    def __call__(self,X,t):
        """
        X is the phase-space vector [x,y,z,x',y',z']
        
        distances are in kpc,
        velocities are in km/s
        time step is assumed to be in yr

        return the time derivative for ode integration
        """
        xvals = X[:3]-self.locs
        rvals = numpy.sqrt( (xvals**2).sum(1) )
        
        dVdt = sum([ self.halos[i].accel(rvals[i])*xvals[i]/rvals[i] \
                         for i in range(self.N) ])
        return numpy.concatenate([X[3:] * 1E3 * yr/kpc,
                                  dVdt])

class fprime_halo_pos:
    def __init__(self,halo1,halo2,dG):
        """
        halo2 is falling toward halo1
        gravitational difference is dG
        """            
        self.halo1 = halo1
        self.halo2 = halo2
        self.dG = dG
        
    def __call__(self,X,t):
        """
        X is the phase-space vector [x,y,z,x',y',z',X,Y,Z,X',Y',Z']

        xyz is the distance from halo2 to the test point
        XYZ is the distance from halo1 to halo2
        X'Y'Z' is the motion of halo2 in reference to halo1
        
        distances are in kpc,
        velocities are in km/s
        time step is assumed to be in yr

        return the time derivative for ode integration
        """
        #distance from halo1 to halo2
        x = X[6:9]
        r = numpy.sqrt(numpy.dot(x,x))

        #distance from halo2 to test particle
        x2 = X[:3]
        r2 = numpy.sqrt(numpy.dot(x2,x2))

        #distance from halo1 to test particle
        x1 = X[:3]+X[6:9]
        r1 = numpy.sqrt(numpy.dot(x1,x1))
        
        #acceleration of halo2, relative to halo1
        dVdt = (1.+self.dG)*self.halo1.accel(r)*x/r

        #acceleration of test particle, relative to halo2
        dvdt = self.halo2.accel(r2)*x2/r2 \
            + self.halo1.accel(r1)*x1/r1 \
            - dVdt
        return numpy.concatenate([X[3:6] * 1E3 * yr/kpc,
                                  dvdt,
                                  X[9:12] * 1E3 * yr/kpc,
                                  dVdt])
        
