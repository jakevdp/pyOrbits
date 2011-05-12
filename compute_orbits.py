"""
Compute orbits for a disk within a dark matter halo
"""
import numpy
import pylab
from time import time
from scipy.integrate import odeint


from pyOrbits.fprime import *
from pyOrbits.potentials.spherical import *
from pyOrbits.disks import get_initial_params
from pyOrbits.integrate import integrate_orbits
from pyOrbits import plotting


def calc_params(halo,
                c,
                rho_crit = 140.):
    try:
        rv = c*halo.r0
    except:
        print "Warning: halo has no r0: setting rv to 40"
        rv = 40.
    Mv = halo.mass(rv)
    Nvir = Mv/(4*pi/3.*rho_crit*rv*rv*rv)
    vcirc = halo.vc(rv)

    return dict(rho0 = halo.rho0,
                r0 = halo.r0,
                Mv = Mv,
                rv = rv,
                Nvir = Nvir,
                vcirc = vcirc)

    
              
def compute_orbits(halo1,
                   halo2,
                   d0   = 140.,     #kpc
                   v0   = 0.,       #km/s : initial velocity of 
                   tmax   = 1E9,    #yr
                   Nsteps = 10,     #number of time steps
                   dG   = 0.3,      #(unitless)
                   #                #halo toward larger one
                   theta = pi/4,    #angle of disk to infall (radians)
                   rmin  = 0.0,
                   rdisk = 4.,      #radius of disk (kpc)
                   Nstars = 10,     #number of stars
                   **kwargs):
    fprime = fprime_halo_pos(halo1,halo2,dG)

    #pos and vel are shape (3,Nstars)
    pos,vel = get_initial_params(halo2.vc,
                                 Nstars = Nstars,
                                 theta = theta,
                                 rdisk = rdisk,
                                 rmin = rmin,
                                 **kwargs)

    t = numpy.linspace(0,tmax,Nsteps+1)

    pos_all = numpy.zeros((6,Nstars))
    vel_all = numpy.zeros((6,Nstars))

    #switched pos and vel: fix this later
    pos_all[:3] = pos
    pos_all[3:] = vel
    vel_all[2] = d0
    vel_all[5] = -v0

    return integrate_orbits(fprime, pos_all, vel_all, t)


def compute_orbits_static(halo1,
                          halo2,
                          d0   = 140.,     #kpc
                          tmax   = 1E9,    #yr
                          Nsteps = 10,     #number of time steps
                          dG   = 0.3,      #(unitless)
                          #                #halo toward larger one
                          theta = pi/4,    #angle of disk to infall (radians)
                          rmin  = 0.0,
                          rdisk = 4.,      #radius of disk (kpc)
                          Nstars = 10,     #number of stars
                          **kwargs):
    halo1.M *= -dG #make halo repel with the correct force

    fprime = fprime_halos((halo1,(0,0,-d0)),(halo2,(0,0,0)))

    #pos and vel are shape (3,Nstars)
    pos,vel = get_initial_params(halo2.vc,
                                 Nstars = Nstars,
                                 theta = theta,
                                 rdisk = rdisk,
                                 rmin = rmin,
                                 **kwargs)

    t = numpy.linspace(0,tmax,Nsteps+1)

    return integrate_orbits(fprime, pos, vel, t)




def main1():
    Nstars = 100
    scatter_vel=(1,1,0.5)
    scatter_pos=(0,0,0.05)
    theta = pi/2
    
    #fiducial parameters
    #halo1
    M=1.6E11

    #halo2
    rho0 = 9E6#1E6
    r0 = 4.
    c = 10.
    #halo2 = cSIS(rho0,r0)
    halo2 = NFW(rho0,r0)
    #halo2 = SIS(20)

    #infall
    d0=140
    tmax=1E9
    Nsteps=100
    
    #gravity modification
    dG=1.0

    #modify any parameters?

    halo1 = point_mass(M)
    

    #------------------------------------------------------------
    params = calc_params(halo2,c,rho_crit = 140.)

    halo_params = ["rho_0  = %.2g Msun/kpc^3" % params['rho0'],
                   "r_0    = %.2g kpc" % params['r0'],
                   "M_v    = %.2g Msun" % params['Mv'],
                   "r_v    = %.2g kpc" % params['rv'],
                   "N_vir  = %.0f"     % params['Nvir'],
                   "v_circ = %.2g km/s at virial radius" % params['vcirc']]

    print "halo parameters:"
    print '   | '+'\n   | '.join(halo_params)
            
    numpy.random.seed(45)
    print "computing orbits for theta=%2.gpi, dG=%.2g" % (theta/pi,dG)
    t0 = time()
    orbits = compute_orbits(halo1,
                            halo2,
                            d0 = d0,
                            tmax = tmax,
                            Nsteps = Nsteps,
                            dG = dG,
                            theta = theta,
                            Nstars = Nstars,
                            scatter_vel = scatter_vel,
                            scatter_pos = scatter_pos)
    tf = time()-t0
    print "computation time:"
    print " - %.2g sec for %i orbits (%.3g orbits/sec)" \
        % (tf,
           orbits.shape[0],
           orbits.shape[0]/tf)
    
    d0 = orbits[0,0,8]
    df = orbits[0,-1,8]
    vf = orbits[0,-1,11]
    
    print "infall parameters:"
    print "   | M1 = %.2g Msun" % M
    print "   | d0 = %.0f kpc" % d0
    print "   | df = %.0f kpc" % df
    print "   | vf = %.0f km/s" % vf
    
    outfile = 'orbits/orbit.npz'

    numpy.savez(outfile,
                orbits = orbits,
                dG = dG,
                d0 = d0,
                df = df,
                M1 = M,
                theta = theta,
                tmax = tmax,
                **params)
    
    fig1,fig2 = plotting.plot_orbits_3D(outfile,Nstars)
    fig1.get_axes()[0].set_zlim3d(-2,2)
    fig2.get_axes()[0].set_zlim3d(-1,3)

    pylab.figure()
    plotting.plot_infall(outfile)
        
    pylab.show()



def main2():
    Nstars = 4000
    scatter_vel=(1,1,0.5)
    scatter_pos=(0,0,0.05)
    theta = pi
    
    #fiducial parameters
    #halo1
    M=1.6E11

    #halo2
    rho0 = 9E6#1E6
    r0 = 4.
    c = 10.
    #halo2 = cSIS(rho0,r0)
    halo2 = NFW(rho0,r0)
    #halo2 = SIS(20)

    #infall
    d0=140
    tmax=1E9
    Nsteps=100
    
    #gravity modification
    dG=1.0

    #modify any parameters?

    halo1 = point_mass(M)
    

    #------------------------------------------------------------
    params = calc_params(halo2,c,rho_crit = 140.)

    halo_params = ["rho_0  = %.2g Msun/kpc^3" % params['rho0'],
                   "r_0    = %.2g kpc" % params['r0'],
                   "M_v    = %.2g Msun" % params['Mv'],
                   "r_v    = %.2g kpc" % params['rv'],
                   "N_vir  = %.0f"     % params['Nvir'],
                   "v_circ = %.2g km/s at virial radius" % params['vcirc']]

    print "halo parameters:"
    print '   | '+'\n   | '.join(halo_params)
            
    numpy.random.seed(45)
    print "computing orbits for theta=%2.gpi, dG=%.2g" % (theta/pi,dG)
    t0 = time()
    orbits = compute_orbits_static(halo1,
                                   halo2,
                                   d0 = d0,
                                   tmax = tmax,
                                   Nsteps = Nsteps,
                                   dG = dG,
                                   theta = theta,
                                   Nstars = Nstars,
                                   scatter_vel = scatter_vel,
                                   scatter_pos = scatter_pos)
    tf = time()-t0
    print "computation time:"
    print " - %.2g sec for %i orbits (%.3g orbits/sec)" \
        % (tf,
           orbits.shape[0],
           orbits.shape[0]/tf)
    
    print "infall parameters:"
    print "   | M1 = %.2g Msun" % M
    
    outfile = 'orbits/orbit.npz'

    numpy.savez(outfile,
                orbits = orbits,
                dG = dG,
                M1 = M,
                theta = theta,
                tmax = tmax,
                **params)
    
    fig1,fig2 = plotting.plot_orbits_3D(outfile,Nstars)
    fig1.get_axes()[0].set_zlim3d(-2,2)
    fig2.get_axes()[0].set_zlim3d(-1,3)
        
    pylab.show()

if __name__ == '__main__':
    main1()
