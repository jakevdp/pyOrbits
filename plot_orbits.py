import numpy
import pylab
from mpl_toolkits.mplot3d import Axes3D

from pyOrbits.plotting import *

if __name__ == '__main__':
    filename = 'orbits/orbit.npz'
    #filename = 'orbits/NFW_G1.0.npz'
    #plot_init_final(filename)
    #pylab.show()
    #exit()
    
    #filename = 'orbits/NFW_c10_rho9E6.npz'
    #filename = 'orbits/cSIS_r3_rho1E7_G1.npz'
    
    fig1,fig2 = plot_orbits_3D(filename,4000)
    #fig2.get_axes()[0].set_zlim3d(-1,3)
    fig1.get_axes()[0].set_zlim3d(-3,3)

    #pylab.figure()
    #plot_final_shape(filename)

    pylab.show()
