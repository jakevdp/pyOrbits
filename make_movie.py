from pyOrbits.animate import *

if __name__ == '__main__':
    make_movie(filename='orbits/orbit.npz',
               outfile='animations/NFW_spiral.mpg')
    #make_movie(filename='orbits/cSIS_r3_rho1E7_G1.npz',
    #           outfile = 'animations/cSIS_r3_rho1E7_G1.mpg')
    #make_movie(filename='orbits/NFW_c10_rho9E6.npz',
    #           outfile='animations/NFW_c10_rho9E6.mpg')
