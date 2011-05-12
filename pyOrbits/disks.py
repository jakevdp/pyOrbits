import numpy

def get_initial_params(vc_func,
                       Nstars,
                       theta,
                       rdisk,
                       rmin=0.0,
                       scatter_pos = [0,0,0],
                       scatter_vel = [0,0,0]):
    """
    |Given a halo potential, return the initial conditions positions and
    | velocities in the tilted cylindrical coordinates
    |
    |vc_func: function that takes a radius in kpc and returns
    | the circular velocity in km/s
    |
    |Nstars:  number of stars to return
    |
    |theta (radians): inclination angle of the disk axis relative to z-axis
    |
    |rdisk (kpc): outer radius of disk
    |rmin  (kpc): inner radius of disk. In some potentials, integrating 
    |             close-in orbits causes problems
    |
    |scatter_pos gives the magnitude of positional scatter (kpc)
    |scatter_vel gives the magnitude of velocity scatter (km/s)
    |  both scatter_pos and scatter_vel are in the disk frame
    |   cylindrical coordinates: [r,phi,z]
    """
    #------------------------------------------------------------
    #start by populating the locations and velocities in a 
    # polar coordinate system, centered on the disk
    pos_rpz = numpy.zeros( (3,Nstars) )
    vel_rpz = numpy.zeros( (3,Nstars) )

    pos_rpz[0] = rmin + (rdisk-rmin)*numpy.random.random(Nstars)
    pos_rpz[1] = 2*numpy.pi*numpy.random.random(Nstars)

    vel_rpz[1] = vc_func(pos_rpz[0])

    #------------------------------------------------------------
    #add scatter to positions and velocities if specified
    for i in range(3):
        if scatter_pos[i]>0:
            pos_rpz[i] += numpy.random.normal(scale=scatter_pos[i],
                                              size = Nstars)
        if scatter_vel[i]>0:
            vel_rpz[i] += numpy.random.normal(scale=scatter_vel[i],
                                              size = Nstars)

    #------------------------------------------------------------
    #convert polar coordinates to cartesian coordinates
    pos_xyz = numpy.zeros( (3,Nstars) )
    vel_xyz = numpy.zeros( (3,Nstars) )

    #--------------------
    #positions:
    # x = R cos(phi)
    # y = R sin(phi)
    # z = z
    cos_phi = numpy.cos( pos_rpz[1] )
    sin_phi = numpy.sin( pos_rpz[1] )

    pos_xyz[0] = pos_rpz[0]*cos_phi
    pos_xyz[1] = pos_rpz[0]*sin_phi
    pos_xyz[2] = pos_rpz[2]

    #--------------------
    #velocities:
    #vx = vR cos(phi) - vp sin(phi)
    #vy = vp cos(phi) + vR sin(phi)
    #vz = vz    
    vel_xyz[0] = vel_rpz[0]*cos_phi - vel_rpz[1]*sin_phi
    vel_xyz[1] = vel_rpz[1]*cos_phi + vel_rpz[0]*sin_phi
    vel_xyz[2] = vel_rpz[2]

    #------------------------------------------------------------
    #rotate these cartesian coordinates to a new coordinate system,
    #  which is rotated about the y-axis by an angle theta
    R = numpy.array( [ [ numpy.cos(theta),   0,   numpy.sin(theta) ],
                       [ 0,                  1,   0                ],
                       [ -numpy.sin(theta),  0,   numpy.cos(theta) ] ] )

    pos2_xyz = numpy.dot(R,pos_xyz)
    vel2_xyz = numpy.dot(R,vel_xyz)
    
    return pos2_xyz,vel2_xyz
