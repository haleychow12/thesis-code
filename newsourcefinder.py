# use (/) for floating point division and (//) for integer division
from __future__ import division

import numpy as np
import scipy as scipy
import random as random
import math
import matplotlib.pyplot as plt
import sys

"""the magnetic field strength of an avalanche transceiver is 
required to be between 0.5 uA/m and 2.16 uA/m at a distance of 10m

the magnetic field vector (in uA/m) at a given location.

                1     3 r (m * r)        m
        H(x) = ---- ( -----------  -  ------- ) A/m (Wikipedia: Magnetic Dipole)
               4*pi     ||r||^5       ||r||^3
"""

#Point class, each point is represented by an x value, and y value
#The point class also contains an array that has an error value 
#at some sample theta
class Point:
    def __init__(self, xvalue, yvalue):
        self.x = xvalue
        self.y = yvalue
        self.error = (-1,-1)

    def setError(self, e, degree):
        #self.error.append((e, degree))
        self.error = (e, degree)

#rotates x and y around 0,0 by theta radians
def rotate(x, y, T):
    return (x*math.cos(T) - y*math.sin(T), x*math.sin(T) + y*math.cos(T))

#rotates x,y about rotate_point_x, rotate_point_y by theta radians
def rotate_point(x, y, rotate_point_x, rotate_point_y, T):
    a = rotate(x-rotate_point_x, y-rotate_point_y, T) 
    return (a[0] + rotate_point_x, a[1] + rotate_point_y)

#returns the value of the magnetic moment
def calc_magnetic_moment():
	power = 1
	moment_min, moment_max = (0.00628318531, 0.02714336053)
	return moment_min + power * (moment_max - moment_min)

#Returns a tuple [x,y] with the value of the electric field 
#tloc (tuple): location of the transmitter
#myloc (tuple): location where the field will be solved for
#theta : in degrees, the orientation of the field
#m : magnetic moment constant
def field(tloc, myloc, theta, m):
	theta = math.radians(theta)
	x,y = rotate_point(myloc[0], myloc[1], tloc[0], tloc[1], -theta)

	#find the real location
	rloc = [x, y]
	m = [0, m]

	r = rloc - tloc

	r_norm = np.linalg.norm(r) #magnitude of the r vector

	f1 = 1 / (4 * scipy.pi)
	f2 = 3 * r * np.dot(m, r) / (r_norm ** 5) #np.dot = rx*mx + ry*my
	f3 = m / (r_norm ** 3)
	f = f1 * (f2 - f3)

	return f * (10 ** 6)

#Plot the vector field of the transmitter (in the xy-plane).
#tloc (tuple): location of the transmitter
#xlim (tuple): (x_min, x_max) representing the x-range of the plot.
#ylim (tuple): (y_min, y_max) representing the y-range of the plot.
#n (int): number of points between each of (x_min, x_max) and
#    (y_min, y_max), inclusive.
def plot_field(tloc, theta, m, xlim, ylim, n):
	x0, x1 = xlim
	y0, y1 = ylim

	nc = n * 1j
	Y, X = np.ogrid[y0:y1:nc, x0:x1:nc]

	xfield = np.zeros((n, n))
	yfield = np.zeros((n, n))
	mag = np.zeros((n, n))
	for i in range(n):
	    for j in range(n):
	        location = np.array([X[0][i], Y[j][0]])
	        (x, y) = field(tloc, location, theta, m)

	        xfield[j][i] = x
	        yfield[j][i] = y
	        mag[j][i] = np.linalg.norm(np.array([x, y]))

	plt.streamplot(X, Y, xfield, yfield)

	#plt.show()

#Returns the straight line distance between loc1 and loc2
#loc1 (tuple): location of the source
#loc2 (tuple): location of the destination
def distance(loc1, loc2):
    return np.linalg.norm(loc1 - loc2)

    

#Returns a tuple that contains the x,y location of the next step
#tloc (tuple): location of the transmitter
#myloc (tuple): location of the searcher
#theta : in degrees
#stepsize: self explanatory
def step(tloc, myloc, theta, stepsize = 0.1):
	(lx, ly) = myloc
	moment = calc_magnetic_moment()
	(fx, fy) = field(tloc, myloc=myloc, theta=theta, m=moment)
	c = stepsize / np.linalg.norm(np.array([fx, fy]))
	newloc = (lx + c * fx, ly + c * fy)

	return newloc

#updates the bestGuess array to hold transmitter locations in order
#bestGuess (list of tuples): list of points with the smallest error value
#	organized by (error, point, degree)
#error: error calculated when p is the transmitter location
#p: Point where the test transmitter is located
#degree: Theta that the test transmitter was oriented in
def storeGuess(bestGuess, error, p, degree):
    bestGuessLength = 4

    #append if list is too short
    #errorList.append((error, p, degree))
    #p.setError(error, degree)
    count = len(bestGuess)

    if count < bestGuessLength:
        bestGuess.append((error, p, degree))
        bestGuess.sort(key = lambda t: t[0])
    #else check if smaller than the largest item in the list
    elif error < bestGuess[count-1][0]:
        bestGuess.append((error, p, degree))
        bestGuess.sort(key = lambda t: t[0])
        del bestGuess[count]

    #print "Source: %.4f, %.4f, Error: %.4f, Degree: %d" % (p.x,p.y,error, degree);

#xguess, yguess = findAvg(bestGuess)
#Calculates the average of the bestGuess array and prints findings
#bestGuess (list of tuples): list of points with the smallest error value
#	organized by (error, point, degree)
#stepNumber: the number of steps that were taken for this data
def findAvg(bestGuess, stepNumber):
	avgx = 0
	avgy = 0

	errorThreshold = .75

	for e,p,d in bestGuess:
	    print("Best Guesses x: %.4f, y: %.4f, error: %.2f, degree: %d" % (p.x, p.y, e, d))
	    avgx += p.x
	    avgy += p.y

	    if e > errorThreshold:
	    	print "error too large"
	    	return 3, 3

	#check if less than 5 iterations of the loop 
	print ("Taken %d steps" % (stepNumber+1))
	if stepNumber < 4:
	    print "started too close"
	    return 3, 3

	xguess = avgx/len(bestGuess)
	yguess = avgy/len(bestGuess)
	print("Source guess: %.4f, %.4f" % (xguess,yguess))

	return xguess, yguess

#draws a plot with a transmitter at tloc, with searchx and searchy locations
#marked by x, points sized by magnitude of their error and with the bestGuess
#locations marked with a *
def drawplot(tloc, searchx, searchy, pointsList, bestGuess):
	m = calc_magnetic_moment()
	for count in range(0, 36):
	    print count*5
	    #test = Transmitter(1, location=np.array([xguess, yguess, 0]), theta=count*5)
	    theta = count*5
	    plot_field(tloc, theta, m, xlim = (-2, 2), ylim = (-1, 1), n = 100)
	    #draw search path
	    for x,y in zip(searchx, searchy):    
	        plt.annotate("x", (x, y))

	    #draw the points with fontsize modulated by magnitude of the error
	    # for p in pointsList:
	    #     e,d = zip(*p.error)
	    #     #x = min(e)
	    #     x = e[count]
	    #     if (x > 6):
	    #         size = 30
	    #     else: 
	    #         size = (x * 5)
	    #     plt.annotate(".", (p.x, p.y), fontsize = size)
	        #print ("x: %.4f, y: %.4f, error: %.2f" % (p.x, p.y, x))

	    #draws the top 4 guesses and places the source at their average
	    for e,p,d in bestGuess:
	        plt.annotate("*", (p.x, p.y))
	    #plt.annotate("o", (xguess, yguess))

	    plt.xlim([-2, 2])
	    plt.ylim([-1, 1])

	    plt.show()

def quickdrawplot(tloc, theta, bestGuess):
	m = calc_magnetic_moment()
	plot_field(tloc, theta, m, xlim = (-20, 20), ylim = (-10, 10), n = 1000)
	#draws the top 4 guesses and places the source at their average
	for e,p,d in bestGuess:
	    plt.annotate("*", (p.x, p.y))
	plt.annotate("o", (tloc[0], tloc[1]))
	plt.xlim([-20, 20])
	plt.ylim([-10, 10])
	plt.show()

def fillPointsList(xlim, ylim, xPoints, yPoints):
	pointsList = []
	pixelX = np.linspace(xlim[0], xlim[1], xPoints)
	pixelY = np.linspace(ylim[0], ylim[1], yPoints)

	for i in range(0, len(pixelX)):
		pointsList.append([])
		for j in range(0, len(pixelY)):
			pointsList[i].append(Point(pixelX[i], pixelY[j]))

	return pointsList


def getXYVals(xlim, ylim, xPoints, yPoints):
	pixelX = np.linspace(xlim[0], xlim[1], xPoints)
	pixelY = np.linspace(ylim[0], ylim[1], yPoints)

	return pixelX, pixelY

#Returns a tuple [testslope, testr] with the test values computed
#tloc (tuple): location of the test transmitter
#myloc (tuple): location where the field will be solved for
#theta : in degrees, the orientation of the field
#m : magnetic moment constant
def getTestVals(tloc, myloc, theta, m):
	(testx, testy) = field(tloc, myloc=myloc, theta=theta, m=m)
	#find the distance between test and search location k
	testr = distance(myloc, tloc)
	testslope = testy/testx

	return testslope, testr

def findNeighbors(x, y, pointsList):
	neighbors = []
	lenx = len(pointsList)
	leny = len(pointsList[0])
	if (x > 0):
		neighbors.append(pointsList[x-1][y])
		if (y > 0):
			neighbors.append(pointsList[x-1][y-1])
		if (y+1 < leny):
			neighbors.append(pointsList[x-1][y+1])
	if (x+1 < lenx):
		neighbors.append(pointsList[x+1][y])
		if (y > 0):
			neighbors.append(pointsList[x+1][y-1])
		if (y+1 < leny):
			neighbors.append(pointsList[x+1][y+1])

	if (y > 0):
		neighbors.append(pointsList[x][y-1])
	if (y+1 < leny):
		neighbors.append(pointsList[x][y+1])

	return neighbors

def findRandomNeighbor(x, y, pointsList):
	neighbors = []
	lenx = len(pointsList)
	leny = len(pointsList[0])
	if (x > 0):
		neighbors.append(pointsList[x-1][y])
		if (y > 0):
			neighbors.append(pointsList[x-1][y-1])
		if (y+1 < leny):
			neighbors.append(pointsList[x-1][y+1])
	if (x+1 < lenx):
		neighbors.append(pointsList[x+1][y])
		if (y > 0):
			neighbors.append(pointsList[x+1][y-1])
		if (y+1 < leny):
			neighbors.append(pointsList[x+1][y+1])

	if (y > 0):
		neighbors.append(pointsList[x][y-1])
	if (y+1 < leny):
		neighbors.append(pointsList[x][y+1])

	rand = random.random()*len(neighbors)

	#print "oldpoint: (%.4f, %.4f)" % (pointsList[x][y].x,pointsList[x][y].y)
	#print "newpoint: (%.4f, %.4f)" % (neighbors[int(rand)].x,neighbors[int(rand)].y)
	return neighbors[int(rand)]

#normalize the slope error
def snorm(error):
	if (error > 20):
		return 1
	return error/20

#normalize the distance error
def dnorm(error):
	if (error > 10):
		return 1
	return error/10


#Returns a tuple that estimates the x and y location of the transciever using
#search data
#tloc (tuple): original location of the transmitter
#myloc (tuple): location of the searcher
#theta : in degrees
#stepsize : self explanatory
#steps : number of steps to take
#visualize : boolean that determines whether plots are drawn
def bruteforce_search(tloc, myloc, theta, radius = .3, stepsize = 0.1, steps = 5, visualize = False):
	searchx = []
	searchy = []
	r = []
	
	slopeList = []
	bestGuess = []

	pointsList = fillPointsList(xlim = (-2,2), ylim = (-1,1), xPoints=100, yPoints=50)
	interval = 2
	m = calc_magnetic_moment()

	for i in range(steps):
		del bestGuess[:]
		(x,y) = myloc
		searchx.append(x)
		searchy.append(y)
		dist = distance(myloc, tloc)
		r.append(dist)

		if (i > 0):
			#solve for the slope
			slope = (y-searchy[i-1])/(x-searchx[i-1])
			slopeList.append(slope)

			for p in pointsList:
				#check if the point is on the wrong side of the line
				del p.error[:]
				for degree in range(0, 180, interval):
					#test = Transmitter(1, location=np.array([p.x, p.y, 0]), theta=degree)
					testloc = np.array([p.x,p.y])
					error = 0
					for k in range(0, len(slopeList)):
						#location is the searchx and searchy
						searchloc = np.array([searchx[k], searchy[k]])
						(testx, testy) = field(testloc, myloc=searchloc, theta=degree, m=m)

						#find the distance between test and search location k
						testr = distance(searchloc, testloc)
						testslope = testy/testx

						#need to compare testy/testx with slope at last location
						#using RMSE 
						error += math.sqrt((slopeList[k] - testslope)**2 + (r[k] - testr)**2)

					if (i > 1):
					    storeGuess(bestGuess, error, p, degree)
					    #print ("x: %.4f, y: %.4f, error: %.2f" % (p.x, p.y, error))
		if dist < radius:
			break 
		myloc = step(tloc, myloc=myloc, theta=theta, stepsize=stepsize)

	if visualize:
		drawplot(tloc, searchx=searchx, searchy=searchy, pointsList=pointsList, bestGuess=bestGuess)

	return findAvg(bestGuess, i)

#Returns a tuple that estimates the x and y location of the transciever using
#search data
#tloc (tuple): original location of the transmitter
#myloc (tuple): location of the searcher
#theta : in degrees
#stepsize : self explanatory
#steps : number of steps to take
#visualize : boolean that determines whether plots are drawn
def tiered_search(tloc, myloc, theta, radius = .3, stepsize = 0.1, steps = 5, visualize = False):
	searchx = []
	searchy = []
	r = []
	
	slopeList = []
	bestGuess = []

	pointsList = fillPointsList(xlim =(-2,2), ylim =(-1,1), xPoints=80, yPoints=40)
	interval = 10
	m = calc_magnetic_moment()

	for i in range(steps):
		del bestGuess[:]
		(x,y) = myloc
		searchx.append(x)
		searchy.append(y)
		dist = distance(myloc, tloc)
		r.append(dist)

		if (i > 0):
			#solve for the slope
			slope = (y-searchy[i-1])/(x-searchx[i-1])
			slopeList.append(slope)

			for p in pointsList:
				del p.error[:]
				for degree in range(0, 180, interval):
					#test = Transmitter(1, location=np.array([p.x, p.y, 0]), theta=degree)
					testloc = np.array([p.x,p.y])
					error = 0
					for k in range(0, len(slopeList)):
						#location is the searchx and searchy
						searchloc = np.array([searchx[k], searchy[k]])
						(testx, testy) = field(testloc, myloc=searchloc, theta=degree, m=m)

						#find the distance between test and search location k
						testr = distance(searchloc, testloc)
						testslope = testy/testx

						#need to compare testy/testx with slope at last location
						#using RMSE 
						error += math.sqrt((slopeList[k] - testslope)**2 + (r[k] - testr)**2)

					if (i > 1):
					    storeGuess(bestGuess, error, p, degree)
					    #print ("x: %.4f, y: %.4f, error: %.2f" % (p.x, p.y, error))
		if dist < radius:
			break 
		myloc = step(tloc, myloc=myloc, theta=theta, stepsize=stepsize)

	xguess, yguess = findAvg(bestGuess, i)
	if (xguess == 3 and yguess == 3):
		return (3,3)
	print ("Tier 1")
	print ("dist: %f" % np.linalg.norm((xguess, yguess) - tloc))
	
	#2nd tier of searching, not necessary if guess < .01 away
	pointsList = fillPointsList(xlim=(xguess-.2,xguess+.2), ylim=(yguess-.1, yguess+.1), xPoints=30, yPoints=15)
	for p in pointsList:
		for degree in range(0, 180, 2):
			testloc = np.array([p.x,p.y])
			error = 0
			for k in range(0, len(slopeList)):
				searchloc = np.array([searchx[k], searchy[k]])
				(testx, testy) = field(testloc, myloc=searchloc, theta=degree, m=m)
				testr = distance(searchloc, testloc)
				testslope = testy/testx
				error += math.sqrt((slopeList[k] - testslope)**2 + (r[k] - testr)**2)
				storeGuess(bestGuess, error, p, degree)

	if visualize:
		drawplot(tloc, searchx=searchx, searchy=searchy, pointsList=pointsList, bestGuess=bestGuess)

	print("Tier 2")
	return findAvg(bestGuess, i)

#visualize : boolean that determines whether plots are drawn
def crawling_search(tloc, myloc, theta, radius = .3, stepsize = 0.2, steps = 7, visualize = False):
	searchx = []
	searchy = []
	r = []
	
	slopeList = []
	bestGuess = []

	#Accuracy, within .1 on x: (-2,2) y:(-1,1) with 100 xpoints and 100 ypoints
	# ~1.08 on x:(-20, 20) y:(-10, 10) with 500 xpoints and 500ypoints
	pointsList = fillPointsList(xlim = (-20,20), ylim = (-10,10), xPoints=500, yPoints=500)
	xList, yList = getXYVals(xlim = (-20,20), ylim = (-10,10), xPoints=500, yPoints=500)

	interval = 2
	samples = 50
	m = calc_magnetic_moment()

	for i in range(steps):
		(x,y) = myloc
		searchx.append(x)
		searchy.append(y)
		dist = distance(myloc, tloc)
		r.append(dist)
		#plt.annotate("x", (x,y))

		if (i > 0):
			#solve for the slope
			slope = (y-searchy[i-1])/(x-searchx[i-1])
			slopeList.append(slope)

		if (i > 4): #for every step
			del bestGuess[:]
			for degree in range(0, 180, interval):
				#errorList = np.zeros(len(pointsList)) 

				#check (#sample) points
				errormin = sys.maxint
				for count in range(0, samples):
					testxIndex = int(random.random()*len(pointsList))
					testyIndex = int(random.random()*len(pointsList[0]))
					testPoint = pointsList[testxIndex][testyIndex]

					#print ("Test x: %.4f, y: %.4f" % (testPoint.x, testPoint.y))
					testloc = np.array([testPoint.x,testPoint.y])
					error = 0
					#test every point
					for k in range(0, len(slopeList)):
						(testslope, testr) = getTestVals(testloc, myloc=np.array([searchx[k], searchy[k]]), theta=degree, m=m)
						error += math.sqrt((slopeList[k] - testslope)**2 + (r[k] - testr)**2)
					pointsList[testxIndex][testyIndex].setError(error, degree)

					if (error < errormin):
						errormin = error
						crawlstart = testPoint

				newmin = sys.maxint
				while (True): #some condition
					#check bounds of the x and y
					x = xList.tolist().index(crawlstart.x)
					y = yList.tolist().index(crawlstart.y)

					#compute the error at all the neighbors
					neighbors = findNeighbors(x, y, pointsList)
					for p in neighbors:
						testloc = np.array([p.x,p.y])
						e,d = p.error
						if (d == degree):
							error = e
						else:
							error = 0
							#test every point
							for k in range(0, len(slopeList)):
								(testslope, testr) = getTestVals(testloc, myloc=np.array([searchx[k], searchy[k]]), theta=degree, m=m)
								error += math.sqrt((slopeList[k] - testslope)**2 + (r[k] - testr)**2)
								#4*the slope part might help
							p.setError(error, degree)
							storeGuess(bestGuess, error, p, degree)
						#store the values from neighbor with smallest error
						if (error < newmin):
							newmin = error 
							newcrawl = p

					#check against old errormin value
					if (newmin >= errormin):
						break
					crawlstart = newcrawl
					errormin = newmin

	 	myloc = step(tloc, myloc=myloc, theta=theta, stepsize=stepsize)

	if visualize:
		#quickdrawplot(tloc, theta, bestGuess)
		e,p,d = bestGuess[0]
		quickdrawplot(np.array([p.x, p.y]), d, bestGuess)

	return findAvg(bestGuess, i)

def annealing_search(tloc, myloc, theta, radius = .3, stepsize = .2, steps = 8, visualize = False):
	searchx = []
	searchy = []
	r = []
	slopeList = []
	bestGuess = []

	pointsList = fillPointsList(xlim = (-20,20), ylim = (-10,10), xPoints=500, yPoints=500)
	xList, yList = getXYVals(xlim = (-20,20), ylim = (-10,10), xPoints=500, yPoints=500)

	interval = 2
	samples = 50
	const = 16
	m = calc_magnetic_moment()
	serrorList = []
	derrorList = []

	for i in range(steps):
		(x,y) = myloc
		searchx.append(x)
		searchy.append(y)
		dist = distance(myloc, tloc)
		r.append(dist)
		#plt.annotate("x", (x,y))

		if (i > 0):
			#solve for the slope
			slope = (y-searchy[i-1])/(x-searchx[i-1])
			slopeList.append(slope)

		if (i > 6): #for every step
			for degree in range(0, 180, interval):
				#annealing parameters
				#print degree
				T = 9000 
				alpha = .9
				jmax = 	5000
				errormax = .0001
				

				errormin = sys.maxint
				for count in range(0, samples):
					startxIndex = int(random.random()*len(pointsList))
					startyIndex = int(random.random()*len(pointsList[0]))
					testPoint = pointsList[startxIndex][startyIndex]
					testloc = np.array([testPoint.x, testPoint.y])

					#get cost of initial point
					olderror = 0
					slopeerror = 0
					disterror = 0
					for k in range(0, len(slopeList)):
						(testslope, testr) = getTestVals(testloc, myloc=np.array([searchx[k], searchy[k]]), theta=degree, m=m)
						#olderror += math.sqrt((slopeList[k] - testslope)**2 + (r[k] - testr)**2)
						slopeerror += (slopeList[k] - testslope)**2
						disterror += (r[k] - testr)**2
					olderror = math.sqrt(const*snorm(slopeerror) + dnorm(disterror))

					#print ("slope error: %.4f, dist error: %.4f" % (slopeerrorr, disterror))
					# serrorList.append(slopeerror)
					# derrorList.append(disterror)
					pointsList[startxIndex][startyIndex].setError(olderror, degree)


					if (olderror < errormin):
						errormin = olderror
						startPoint = testPoint

				j = 1
				while j <= jmax and olderror > errormax:
					x = xList.tolist().index(startPoint.x)
					y = yList.tolist().index(startPoint.y)

					#compute the error at all the neighbors
					#print str(startPoint.x) + "," + str(startPoint.y)
					nextPoint =  findRandomNeighbor(x,y, pointsList)#choose random neighbor
					e,d = nextPoint.error
					if (d == degree):
						nexterror = e
					else: 
						nextloc = np.array([nextPoint.x, nextPoint.y])
						nexterror = 0
						slopeerror = 0
						disterror = 0
						for k in range(0, len(slopeList)):
							(testslope, testr) = getTestVals(nextloc, myloc=np.array([searchx[k], searchy[k]]), theta=degree, m=m)
							#nexterror += math.sqrt((slopeList[k] - testslope)**2 + (r[k] - testr)**2)
							slopeerror += (slopeList[k] - testslope)**2
							disterror += (r[k] - testr)**2
						#print ("(%.4f, %.4f): slope error: %.4f, dist error: %.4f" % (nextPoint.x, nextPoint.y, snorm(slopeerrorr), snorm(disterror)))
						nexterror = math.sqrt(const*snorm(slopeerror) + dnorm(disterror))
						#serrorList.append(slopeerror)
						#derrorList.append(disterror)

						nextPoint.setError(nexterror,degree)
						storeGuess(bestGuess, nexterror, nextPoint, degree)

					#print str(nextPoint.x) + "," + str(nextPoint.y) + ", error: " + str(nexterror)

					#don't want to store same guess more than once

					delta = nexterror - olderror
					if delta < 0:
						startPoint = nextPoint
						olderror = nexterror
					else:
						#print ("%.4f/%.4f" % (-delta, T))
						p = math.exp(-delta/T)
						#print ("%.4f > %.4f, go to with prob p: %.4f" % (nexterror, olderror, p))
						#print p
						if (random.random() < p):
							startPoint = nextPoint
							olderror = nexterror

					T = T*alpha
					#print T
					j += 1


		myloc = step(tloc, myloc=myloc, theta=theta, stepsize=stepsize)

	if visualize:
		#quickdrawplot(tloc, theta, bestGuess)
		e,p,d = bestGuess[0]
		quickdrawplot(np.array([p.x, p.y]), d, bestGuess)

	# print ("AVGS")
	# print (sum(serrorList)/len(serrorList))
	# print (sum(derrorList)/len(derrorList))
	return findAvg(bestGuess, i)




####################################################################   
#conducts the search with source-finding analysis at each step
#def main():
iterations = 1
# estepsize = 2
# xarray = np.zeros(iterations)
# start = .0625
# for i in range(0,len(xarray)):
#  	xarray[i] = start
#  	start = float(start * estepsize)

# #print xarray
# yarray = np.zeros(iterations)

for j in range(iterations):
	#print xarray[j]
	i = 0
	trials = 1
	avgDistance = np.zeros(trials) 
	while (i < trials):
		#set the beginning parameters

		#Test for java x: -14.7067, y: 1.8595, theta: 66.12
	    theta = 66.12#random.random()*180
	    sourcex = -14.7067#random.random()*40 - 20
	    sourcey = 1.8595#random.random()*20 - 10

	    # Define a transmitter
	    tloc = np.array([sourcex, sourcey])

	    #Determine random starting location for searcher
	    startx = 15
	    starty = -7.5

	    print ("x: %.4f, y: %.4f, theta: %.2f" % (sourcex, sourcey, theta))

	    # Define a searcher
	    myloc = np.array([startx, starty])

	    # Perform a search
	    xguess, yguess = annealing_search(tloc, myloc=myloc, theta=theta)
	    
	    #redo iteration if not enough steps
	    if (xguess == 3 and yguess == 3):
	        print("redoing iteration") #if error is too large, might want to add additional step?
	        continue

	    avgDistance[i] = math.sqrt((sourcex - xguess)**2 + (sourcey - yguess)**2)
	    print("Distance between the guess and source is %.4f" % avgDistance[i])
	    i += 1

	print ("Avg Distance: %.4f" % (np.mean(avgDistance)))
	#yarray[j] = np.mean(avgDistance)

# plt.plot(xarray,yarray)
# plt.show()


	   









