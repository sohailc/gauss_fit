"""
A 2D gaussian function is sampled randomely at 10000 points and noise is
added to the resulting dataset. The parameters of the gaussian function 
is also chosen at random. A gaussian fit is performed by a library 
which I wrote in C, which is accessed through ctypes. 

NB: I know a gauss fit could also have been performed by 
scipy.optimize.cirve_fit but I chose to use a library which I wrote 
myself to show my skills in C.  

The result of this fit is shown by drawing a voronoi cell around each of 
the 10000 sample points, with the color indicating the height of the of 
the function. A solid green line shows the result of the fit.
"""


from scipy.spatial import Voronoi
from matplotlib.collections import PolyCollection
import matplotlib.pyplot as plt
import numpy
from gauss2DFit import main as gauss2DFit

def visualizeData(ax, intensities, X, Y):

    X = X.reshape(-1,1)
    Y = Y.reshape(-1,1)
    X -= numpy.mean(X)
    Y -= numpy.mean(Y)

    layout = numpy.append(X, Y, axis=1)

    cm = numpy.min(intensities)
    cx = numpy.max(intensities)
    colorValues = [float(c-cm)/(cx-cm) for c in intensities]
    
    minx = numpy.min(layout[:,0])
    maxx = numpy.max(layout[:,0])
    miny = numpy.min(layout[:,1])
    maxy = numpy.max(layout[:,1])

    vor = Voronoi(layout)

    polygons = []
    idx = []

    for count, regionIdx in enumerate(vor.point_region):

        region = vor.regions[regionIdx]
            
        if not -1 in region:
            polygon = [vor.vertices[i] for i in region]

            px, py = zip(*polygon)
            
            if (px > maxx).any() or (px < minx).any() or (py > maxy).any() or (py < miny).any():
                continue

            polygons.append(polygon)
            idx.append(count)

    vals = numpy.array(colorValues)[idx]
    coll = PolyCollection(polygons, array=vals, cmap=plt.cm.jet, edgecolors='none')
    ax.add_collection(coll)
    ax.set_aspect(1)
    mapable=ax.get_children()[2] 
    clbr = plt.colorbar(mapable, ax=ax)
    clbr.set_label("intensity")
    #ax.plot(X, Y, 'xr', label="random samples")


def gaussianFunction(x0, sigx, y0, sigy, flr, A):
	
	def f(X,Y):
		Xrel = X-x0
		Yrel = Y-y0
		Z = A*numpy.exp(-Xrel**2/(2*sigx**2)-Yrel**2/(2*sigy**2)) + flr
		return Z
	
	return f
	
def drawFit(ax, fitParams):
	
	x0 = fitParams["mux"]
	y0 = fitParams["muy"]
	r1 = fitParams["sigx"]
	r2 = fitParams["sigy"]
	
	t = numpy.linspace(0,1,100)
	x = r1*numpy.cos(2*numpy.pi*t) + x0
	y = r2*numpy.sin(2*numpy.pi*t) + y0
	
	ax.plot(x,y,'g', label="gaussian fit")

def test(f, xmin, xmax, ymin, ymax, N, noise):
	
	X = numpy.random.uniform(xmin, xmax, N)
	Y = numpy.random.uniform(xmin, xmax, N)
	Z = f(X,Y) + numpy.random.normal(0.0, noise, N)
	
	fitParams = gauss2DFit(X,Y,Z)
	if not fitParams:
		return
	
	fig, ax = plt.subplots(1,1)
	
	visualizeData(ax, Z, X, Y)
	drawFit(ax, fitParams)
	
	ax.set_xlim([xmin,xmax])
	ax.set_ylim([ymin,ymax])
	ax.legend()
	plt.show()

def main():
	
	xmin = -5
	xmax = 5
	ymin = -5
	ymax = 5
	N = 100**2
	
	x0 = numpy.random.uniform(xmin,xmax)
	y0 = numpy.random.uniform(ymin, ymax)
	sigx = numpy.random.uniform(0.3,3.0)
	sigy = numpy.random.uniform(0.3,3.0)
	
	A = 1.0
	noise = A/50
	flr = 0
	
	f = gaussianFunction(x0, sigx, y0, sigy, flr, A)
	
	test(f,xmin,xmax,ymin,ymax, N, noise)

if __name__ == "__main__":
	main()
