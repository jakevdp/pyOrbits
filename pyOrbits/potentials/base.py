"""
base classes for computing potentials
"""

import numpy

#------------------------------------------------------------
#  Constants used for the potentials
#
Msun = 1.98892E30 #kg
kpc  = 3.08568E19 #m
G    = 6.673E-11  #m^3/kg/s^2
pi   = numpy.pi
yr   = 3.1556926E7 #sec

#------------------------------------------------------------
#  Helper Routines
#
def calc_rho_crit(H = 72):
    """
    given the Hubble constant in km/s/Mpc,
    calculate the critical density in Msun/kpc^3
    """
    rho_crit = (3.*(H/kpc)**2)/(8.*pi*G) #kg/m^3
    return rho_crit * kpc**3/Msun

def display_exp(x):
    s = '%.2g'%x
    L = s.split('e')
    if len(L)==1:
        return L[0]
    elif len(L)==2:
        L[1] = L[1][0] + L[1][1:].lstrip('0')
        if L[1].startswith('+'):
            L[1] = L[1][1:]
        if L[0]=='1' or L[0]=='1.0':
            return r'10^{%s}' % L[1]
        else:
            return r'%s\times 10^{%s}' % (L[0],L[1])



#------------------------------------------------------------
#  Base Classes
#
class spherical_potential(object):
    """
    |sphereical_potential:
    |
    | An abstract base class for a spherically symmetric potential.
    | Derived classes should define three methods:
    |   density(self,r)   : return the denstiy at a radial distance r
    |   mass(self,r)      : return the mass within a sphere of radius r
    |   potential(self,r) : return the potential at a radius r
    |
    """
    def __init__(self):
        pass

    def density(self,r):
        """
        | calculate the density at a radius r
        |  accepts: r (kpc)
        |  returns: rho (Msun/kpc^3)
        """
        raise ValueError, "density undefined for class %s" \
            % self.__class__.__name__

    def mass(self,r):
        """
        | calculate the mass interior to a distance r
        |  accepts: r (kpc)
        |  returns: M (Msun)
        """
        raise ValueError, "mass undefined for class %s" \
            % self.__class__.__name__

    def potential(self,r):
        """
        | calculate the potential due to the halo at distance r
        |  accepts: r (kpc)
        |  returns: Phi (m^2/s^2)
        """
        raise ValueError, "potential undefined for %s" \
            % self.__class__.__name__

    def accel(self,r):
        """
        | calculate the radial acceleration at a distance r
        |  accepts: r (kpc)
        |  returns: a (km/s/yr)
        """
        return -G * (self.mass(r)*Msun) / (r*kpc)**2 * 1E-3*yr

    def vc2(self,r):
        """
        | calculate the square of the circular velocity at a distance r
        |  accepts: r (kpc)
        |  returns: vc^2 ( (km/s)^2 )
        """
        return G * (self.mass(r)*Msun) / (r*kpc) * 1E-6

    def vc(self,r):
        """
        | calculate the circular velocity at a distance r
        |  accepts: r (kpc)
        |  returns: vc (km/s)
        """
        return numpy.sqrt(self.vc2(r))

    def __str__(self):
        return self.__class__.__name__

class axisym_potential(object):
    """
    |axisym_potential:
    |
    | An abstract base class for an axi-symmetric potential.
    | Derived classes should define three methods:
    |   density(self,R,z)   : return the denstiy at position R,z
    |   accel(self,R,z)     : return the acceleration at position R,z
    |   potential(self,R,z) : return the potential at position R,z
    |
    """
    def __init__(self):
        pass

    def density(self,R,z):
        """
        | return the density at a radius R and height z
        |  accepts: R (kpc), z (kpc)
        |  returns: rho (Msun/kpc^3)
        """
        raise ValueError, "density undefined for %s" \
            % self.__class__.__name__

    def potential(self,R,z):
        """
        | calculate the potential at a radius R and height z
        |  accepts: R (kpc), z (kpc)
        |  returns: Phi (m^2/s^2)
        """
        raise ValueError, "potential undefined for %s" \
            % self.__class__.__name__

    def accel(self,R,z):
        """
        | calculate (R,z) acceleration at a radius R and height z
        |  accepts: R (kpc), z (kpc)
        |  returns: a (km/s/yr)
        """
        raise ValueError, "accel undefined for %s" \
            % self.__class__.__name__

    def effective_potential(self,R,z,Lz):
        """
        | calculate effective potential at a radius R, height z, and
        |   angular momentum Lz
        |  accepts: R (kpc), z (kpc), Lz (kpc km/s)
        |  returns: a (km/s/yr)
        """
        return self.potential(R,z) + 0.5 * (Lz*1E3/R)**2

        
