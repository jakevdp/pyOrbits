import numpy
import pylab

"""
integrators take a function F which, given a vector x,
returns the components of acceleration
"""

def setup(x0,v0,tmax,Nsteps):
    assert len(x0)==len(v0)

    t = numpy.linspace(0,tmax,Nsteps+1)
    h = t[1]-t[0]

    x_t = numpy.zeros((len(t),len(x0)))
    v_t = x_t.copy()

    x_t[0] = x0
    v_t[0] = v0

    return t,h,x_t,v_t
    

def Euler(F,x0,v0,tmax,Nsteps,args=()):
    t,h,x_t,v_t = setup(x0,v0,tmax,Nsteps)

    for i in range(1,Nsteps+1):
        x_t[i] = x_t[i-1] + h*v_t[i-1]
        v_t[i] = v_t[i-1] + h*F(x_t[i-1],*args)

    return x_t,v_t

def DriftKick(F,x0,v0,tmax,Nsteps,args=()):
    t,h,x_t,v_t = setup(x0,v0,tmax,Nsteps)

    for i in range(1,Nsteps+1):
        x_t[i] = x_t[i-1] + h*v_t[i-1]
        v_t[i] = v_t[i-1] + h*F(x_t[i],*args)

    return x_t,v_t

def KickDrift(F,x0,v0,tmax,Nsteps,args=()):
    t,h,x_t,v_t = setup(x0,v0,tmax,Nsteps)

    for i in range(1,Nsteps+1):
        v_t[i] = v_t[i-1] + h*F(x_t[i-1],*args)
        x_t[i] = x_t[i-1] + h*v_t[i]

    return x_t,v_t

def LeapFrog(F,x0,v0,tmax,Nsteps,args=()):
    t,h,x_t,v_t = setup(x0,v0,tmax,Nsteps)

    for i in range(1,Nsteps+1):
        xprime = x_t[i-1] + 0.5*h*v_t[i-1]
        v_t[i] = v_t[i-1] + h*F(xprime,*args)
        x_t[i] = xprime + 0.5*h*v_t[i]
    
    return x_t,v_t

def RungeKutta(F,x0,v0,tmax,Nsteps,args=()):
    t,h,x_t,v_t = setup(x0,v0,tmax,Nsteps)
    
    ndim = len(x0)
    xv_F = lambda x: numpy.concatenate((x[ndim:],F(x[:ndim],*args)))
    
    xv = numpy.concatenate((x0,v0))

    for i in range(1,Nsteps+1):
        k1 = h*xv_F(xv)
        k2 = h*xv_F(xv+k1/2.)
        k3 = h*xv_F(xv+k2/2.)
        k4 = h*xv_F(xv+k3)
        xv += k1/6. + k2/3. + k3/3. + k4/6.
        x_t[i] = xv[:ndim]
        v_t[i] = xv[ndim:]
    
    return x_t,v_t


def VelocityVerlet(F,x0,v0,tmax,Nsteps,args=()):
    t,h,x_t,v_t = setup(x0,v0,tmax,Nsteps)

    F_last = F(x_t[0],*args)
    F_next = None
    
    for i in range(1,Nsteps+1):
        x_t[i] = x_t[i-1] + h*v_t[i-1] + 0.5*F_last*h**2
        F_next = F(x_t[i],*args)
        v_t[i] = v_t[i-1] + 0.5*(F_last+F_next)*h
        F_last = F_next

    return x_t,v_t




def test_integrators(tmax=100,
                     Nsteps=1000,
                     x0=[0.0,1.0],
                     v0=[1.0,0.0]):
    mag = lambda x: numpy.sqrt(numpy.dot(x,x))
    fprime = lambda x: -x/mag(x)**3
    energy = lambda x,v: 0.5*(v**2).sum(1) - 1./numpy.sqrt((x**2).sum(1))

    pylab.figure(1)
    
    for (method,N) in ((Euler,Nsteps),
                       (DriftKick,Nsteps),
                       (LeapFrog,Nsteps),
                       (VelocityVerlet,Nsteps),
                       (RungeKutta,Nsteps/4)):
        t = numpy.linspace(0,tmax,N+1)
        x_t,v_t = method(fprime,x0,v0,tmax,N)
        E = energy(x_t,v_t)

        delta_E = abs((E-E[0])/E[0])
        
        pylab.loglog(t[1:],delta_E[1:],
                     label=method.__name__)
    
    pylab.ylim(1E-8,1E4)
    pylab.legend(loc=2)
    pylab.xlabel(r'$\mathdefault{t}$')
    pylab.ylabel(r'$\mathdefault{|\Delta E/E|}$')
    
    pylab.show()
                    
if __name__ == '__main__':
    test_integrators()
    
