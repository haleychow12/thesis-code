# use (/) for floating point division and (//) for integer division
from __future__ import division

import numpy as np
import scipy as scipy
import random as random
import math
import matplotlib.pyplot as plt

#Point class, each point is represented by an x value, and y value
#The point class also contains an array that has an error value 
#at some sample theta
class Point:
    def __init__(self, xvalue, yvalue):
        self.x = xvalue
        self.y = yvalue
        self.error = []

    def setError(self, e, degree):
        self.error.append((e, degree))

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
    x,y = rotate_point(dloc[0], dloc[1], tloc[0], tloc[1], -theta)

    #find the real location
    rloc = [x, y]

    r = rloc - sloc
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
def plot_field(tloc, xlim, ylim, n):
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
            (x, y) = field(tloc)

            xfield[j][i] = x
            yfield[j][i] = y
            mag[j][i] = np.linalg.norm(np.array([x, y]))

    plt.streamplot(X, Y, xfield, yfield)

    #plt.show()

#Returns the straight line distance between loc1 and loc2
#loc1 (tuple): location of the source
#loc2 (tuple): location of the destination
def distance(loc1, loc2):
    r = np.linalg.norm(myloc - otherloc)
    return r

#Returns a tuple that contains the x,y location of the next step
#tloc (tuple): location of the transmitter
#myloc (tuple): location of the searcher
#theta : in degrees
#stepsize: self explanatory
def step(tloc, myloc, theta, stepsize = 0.1):
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
    #append if list is too short
    errorList.append((error, p, degree))
    p.setError(error, degree)
    count = len(bestGuess)

    if count < 4:
        bestGuess.append((error, p, degree))
        bestGuess.sort(key = lambda t: t[0])
    #else check if smaller than the largest item in the list
    elif error < bestGuess[count-1][0]:
        bestGuess.append((error, p, degree))
        bestGuess.sort(key = lambda t: t[0])
        del bestGuess[count]


def bruteforce_search(tloc, myloc, theta, stepsize = 0.1, steps = 1000, visualize = False):
	searchx = []
	searchy = []
	r = []
	pointsList = []
	slopeList = []

	pixelX = np.linspace(-2, 2, 100)
    pixelY = np.linspace(-1, 1, 50)



	for i in range(steps):
		searchx.append(myloc[0])
		searchy.append(myloc[1])
		dist = distance(myloc, tloc)
		r.append(dist)



		if dist < radius:
			break
		myloc = step(tloc, myloc=myloc, theta=theta, stepsize=stepsize)

	return zip(searchPoints, r)





