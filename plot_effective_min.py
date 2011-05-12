from scipy import optimize
import numpy
import pylab
from pyOrbits.potentials.spherical import *


def find_min(halo,
             Lz, #kpc km/s
             F_z, #m/s^2
             guess=[1.,0.] #guess [R,z] in kpc
             ):
    """
    Find the minimum of the effective potential given an angular
    momentum L_z and z-direction force F_z
    """
    C1 = F_z*kpc
    C2 = 0.5*Lz**2*1E6

    def min_func(y):
        # y=[R,z], both given in kpc
        r = numpy.sqrt(y[0]**2+y[1]**2)
        return halo.potential(r) - C1*y[1] + C2/y[0]**2

    return optimize.fmin(min_func,guess)

def compute_omega_z(halo,
                    R,z):
    """
    compute the vertical SHO frequency for a minimum at (R,z)
    """
    r = numpy.sqrt(R**2+z**2)
    tmp =  G * halo.mass(r) / r**3 * (1.-3*(z/r)**2) \
        + 4*pi*G*halo.density(r) * (z/r)**2

    return numpy.sqrt( 0.5 * tmp*Msun/kpc**3 )
               
def plot_effective_minimum( halos,
                            Npts = 50,
                            M1 = 1.6E11, #Msun
                            dG = 1.0):
    r = 10**numpy.linspace(-2,numpy.log10(4),Npts)
    
    for halo in halos:
        pylab.figure()
        for d in range(120,220,20):
            Fz = dG * G * (M1*Msun) / (d*kpc)**2
            
            Lz = r*halo.vc(r) #kpc km/s
            yopt = numpy.zeros( (Npts,2) )
            for i in range(Npts):
                yopt[i] = abs( find_min(halo,Lz[i],Fz,
                                        guess=(r[i],0.1) ) )

            #non-convergent solutions will run to infinity.
            #set these to NaN
            yopt[numpy.where(yopt>100)]=numpy.nan
            
            omega_z = compute_omega_z(halo,yopt[:,0],yopt[:,1])
        
            pylab.subplot(211)
            pylab.plot(yopt[:,0],yopt[:,1],label='d=%ikpc'%d)
            
            pylab.subplot(212)
            pylab.plot(yopt[:,0],2*pi/omega_z/yr/1E6,label='d=%ikpc'%d)


        pylab.subplot(211)
        pylab.title(halo)
        pylab.ylabel('z (kpc)')
        pylab.ylim(0,pylab.ylim()[1])
        
        pylab.subplot(212)
        pylab.xlabel('R (kpc)')
        pylab.ylabel('T (Myr)')
        pylab.ylim(0,pylab.ylim()[1])
        pylab.legend(loc=4)


if __name__ == '__main__':
    plot_effective_minimum( halos = [ NFW(1E7, 4.0),
                                      cSIS(1E7, 4.0),
                                      cSIS(1E7, 1.0) ] )
    for i in range(3):
        pylab.figure(i+1)
        pylab.savefig('fig/Phi_eff_min_%i.eps' % (i+1))
    
    pylab.show()
