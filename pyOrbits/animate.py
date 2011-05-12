import os
import numpy
import pylab
from mpl_toolkits.mplot3d import Axes3D

#ideas: make this more flexible (ie reference a function that takes
#orbits and an index, so we can animate any function)

def make_movie(filename='orbits/orbits_dG0.3_theta0.00.npz',
               outfile='animations/animation.mpg',
               view = (14,-50),
               istep=1,
               vcodec='mpeg4',
               create_frames = True):
    """
    try vcodec=mpeg4, avi, wmv2
    """
    X = numpy.load(filename)

    orbits = X['orbits']
    tmax   = X['tmax']

    N = orbits.shape[1]

    t = numpy.linspace(0,tmax,N)

    fig = pylab.figure()

    xlim = (0,numpy.max(numpy.sqrt(orbits[:,:,0]**2 + orbits[:,:,1]**2)))
    ylim = (numpy.min(orbits[:,:,2]),numpy.max(orbits[:,:,2]))

    if create_frames:
        for i in range( 0,N,istep ):
            fig.clf()
            
            ax = Axes3D(fig)
            ax.scatter3D(orbits[:,i,0],
                         orbits[:,i,1],
                         orbits[:,i,2])
            ax.view_init(view[0],view[1])
            
            ax.set_xlim3d(-4,4)
            ax.set_ylim3d(-4,4)
            ax.set_zlim3d(-2,4)
            
            ax.set_xlabel('x (kpc)')
            ax.set_ylabel('y (kpc)')
            ax.set_zlabel('z (kpc)')
            
            """
            pylab.subplot(111)
            xf = orbits[:,i,0]
            yf = orbits[:,i,1]
            zf = orbits[:,i,2]
            
            pylab.plot(numpy.sqrt(xf**2+yf**2),zf,'.k')
            pylab.xlabel('r')
            pylab.ylabel('z')
            pylab.xlim(xlim)
            pylab.ylim(ylim)
            """
            
            fig.text(0.8,0.9,
                     't=%i Myr' % (t[i]/1E6),
                     fontsize=14)
            
            frame = 'frames/_frame%03d.png' % (i+1)
            frame2 = 'frames/_frame%03d.png' % (N+i+1)
            print "saving",frame
            fig.savefig(frame)
            os.system('cp %s %s' % (frame,frame2))

    print 'Making movie animation.mpg - this make take a while'
    os.system("mencoder 'mf://frames/_frame*.png' -mf type=png:fps=10 -ovc lavc -lavcopts vcodec=%s -oac copy -o %s" % (vcodec,outfile))
    
