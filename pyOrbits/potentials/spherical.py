import numpy
from scipy.special import hyp2f1 #for NFW_a
from .base import *

class point_mass(spherical_potential):
    def __init__(self,M):
        """
        | point_mass
        |   M is the mass (Msun)
        """
        self.M = M
        self.M_kg = M*Msun

    def density(self,r):
        return 0*r

    def mass(self,r):
        return self.M*r/r

    def potential(self,r):
        return -G * self.M_kg / (r*kpc)

    def __str__(self):
        return r"$\mathdefault{Point\ Mass\ (M=%s M_\odot)}$" \
            % display_exp(self.M)

class SIS(spherical_potential):
    def __init__(self,
                 sigma_v):
        """
        | SIS
        |  sigma_v is the velocity dispersion in km/s
        | the halo has the form rho(r) = sigma_v^2 / (2 pi G r^2)
        """
        
        #self.rho0r0_2 = 4.*pi/3. * rho0 * r0**2
        #self.sigma_v = numpy.sqrt(2*pi*G*rho0r0_2 * Msun/kpc)
            
        self.sigma_v = sigma_v
        self.rho0r0_2 = 0.5 * (sigma_v*1E3)**2 /pi/G * kpc/Msun #Msun/kpc
        
        self.mass_const = 4*pi*self.rho0r0_2 #Msun/kpc
        self.phi_const = G*self.mass_const * Msun/kpc # m^2/s^2

    def density(self,r):
        return self.rho0r0_2 / r**2

    def mass(self,r):
        return self.mass_const * r

    def potential(self,r):
        return self.phi_const * numpy.log(r) #relative to phi at 1kpc

    def __str__(self):
        return r"$\mathdefault{SIS\ (\sigma_v=%.2g km/s)}$" % self.sigma_v

def SIS_(rho0,r0):
    rho0r0_2 = 4.*pi/3. * rho0 * r0**2
    sigma_v =  numpy.sqrt(2*pi*G*rho0r0_2 * Msun/kpc)
    S = SIS(sigma_v)
    S.rho0 = rho0
    S.r0 = r0
    return S

class cSIS(spherical_potential):
    def __init__(self,rho0,r0):
        """
        | cSIS
        |  rho0 is the central density (Msun/kpc^3)
        |  r0 is the core radius (kpc)
        | 
        | this is a halo with the form rho(r) = rho0 r0^2/(r0^2 + r^2)
        """
        self.rho0 = rho0
        self.r0 = r0
        self.M_const = 4*pi*rho0*r0**3 #units: Msun 
        self.Phi_const = 4*pi*G*(rho0*Msun/kpc)*r0*r0 #units: m^2/s^2

    def density(self,r):
        return self.rho0 / (1. + (r/self.r0)**2)

    def mass(self,r):
        x = r/self.r0
        return self.M_const * (x-numpy.arctan(x))

    def potential(self,r):
        x = r/self.r0
        return self.Phi_const * (0.5*numpy.log(x*x+1)+numpy.arctan(x)/x)
    
    def __str__(self):
        return r'$\mathdefault{cSIS(\rho_0=%sM_\odot/kpc^3,r_0=%.2gkpc)}$' \
            % (display_exp(self.rho0),self.r0)
        

class NFW(spherical_potential):
    def __init__(self,rho0,r0):
        """
        | NFW
        |  rho0 is the central density (Msun/kpc^3)
        |  r0 is the core radius (kpc)
        | 
        | this is a halo with the form rho(r) = rho0 /(r/r0)(1 + r/r0)^2
        """
        self.rho0 = rho0
        self.r0 = r0
        self.M_const = 4*pi*rho0*r0**3 #units: Msun 
        self.Phi_const = 4*pi*G*(rho0*Msun/kpc)*r0*r0 #units: m^2/s^2

    def density(self,r):
        x = r/self.r0
        return self.rho0 / (x*(1+x)**2)

    def mass(self,r):
        x = r/self.r0
        return self.M_const * (numpy.log(x+1)-1./(1. + 1./x))
        
    def potential(self,r):
        x = r/self.r0
        return self.Phi_const * (1. - numpy.log(x+1)/x)
    
    def __str__(self):
        return r'$\mathdefault{NFW(\rho_0=%sM_\odot/kpc^3,r_0=%.2gkpc)}$' \
            % (display_exp(self.rho0),self.r0)

class NFW_a(spherical_potential):
    def __init__(self,rho0,r0,alpha):
        self.rho0 = rho0
        self.r0 = r0
        self.alpha = alpha

        if self.alpha==1:
            raise ValueError, "NFW_a cannot use alpha=1.  Use NFW instead"

        self.M_const = 4*pi*rho0*r0**3 #units: Msun 
    
    def density(self,r):
        return self.rho0 / (r/self.r0)**self.alpha \
            / (1+r/self.r0)**(3-self.alpha)

    def mass(self,r):
        x = r/self.r0
        a = 1.*self.alpha
        return self.M_const * x**(1.-a) / (a-1.) * \
            ((x+1.)**(a-2.)*(2.*a*x+a-3.*x-2.)/(a-2.) - hyp2f1(1.-a,1.-a,2.-a,-x))
        
    
    def __str__(self):
        return r'$\mathdefault{NFWa(\rho_0=%sM_\odot/kpc^3,r_0=%.2gkpc,\alpha=%.1f)}$' % (display_exp(self.rho0),self.r0,self.alpha)
