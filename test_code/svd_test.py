import numpy
import pylab

#create random distribution
X = numpy.random.normal(size=(2,4000))
X[0]-=10

#create mixing matrix
M = numpy.random.random([2,2])-0.5

#mix it
MX = numpy.dot(M,X)

#perform SVD to find principal axes
MX_m = MX.mean(1)
MX_c = MX-MX_m[:,None]
U,S,V = numpy.linalg.svd(MX_c,full_matrices=0)

#determine value along principal axis
UMX = numpy.dot(U,MX_c)

#scatterplot, color-coded by principal axis value
pylab.scatter(MX[0],MX[1],c=UMX[0],s=9,lw=0)
pylab.colorbar()

#plot principal axis
x = numpy.linspace(min(MX[0]),max(MX[0]),10)
y = MX_m[1] + (x-MX_m[0])*U[0,1]/U[0,0]
pylab.plot(x,y,'-k')

#set axes limits so x and y scales are equal
xlim = pylab.xlim()
ylim = pylab.ylim()

x0 = numpy.mean(xlim)
y0 = numpy.mean(ylim)
xr = xlim[1]-xlim[0]
yr = ylim[1]-ylim[0]

d = 0.5*max(xr,yr)

pylab.xlim(x0-d,x0+d)
pylab.ylim(y0-d,y0+d)

pylab.show()
