# use (/) for floating point division and (//) for integer division
from __future__ import division

import numpy as np
import scipy as scipy
import random as random
import math
import matplotlib.pyplot as plt

def rotate(x, y, T):
    return (x*math.cos(T) - y*math.sin(T), x*math.sin(T) + y*math.cos(T))

def rotate_point(x, y, rotate_point_x, rotate_point_y, T):
    a = rotate(x-rotate_point_x, y-rotate_point_y, T) 
    return (a[0] + rotate_point_x, a[1] + rotate_point_y)

class Transmitter:
    """Models an avalanche transmitter at the origin."""

    def __init__(self, power, location, theta = 0):
        """Instantiate a transmitter, modeled as a magnetic dipole.

        power (float): A value between 0 and 1 (inclusive) representing the
            transmitted field strength. This argument defines the transmitter's
            power in accordance with the ETS 300 718 avalanche transceiver
            standard, which requires the magnetic field strength to be between
            0,5 uA/m and 2.16 uA/m at a distance of 10m from the transmitter.
            A value of 0 corresponds to the minimum allowable power while a
            value of 1 corresponds to the maximum allowable power.

        The transmitter is modeled as a magnetic dipole at the origin with
        moment c * <0, 1, 0>, where c is set according to the specified power.

        The magnetic dipole moment represents the strength and orientation of
        the transmitter. In particular, one can interpret a magnetic dipole as
        a current loop, so that a magnetic moment <0, 1, 0> specifies a loop
        with current flowing out-of and into the page:

                     \  |  /
                      | | |
                    o -> -> x    (moment of <0, 1, 0>)
                      | | |
                     /  |  \

        theta (value in degrees): determines the orientation of the trasmitter
        which determines the orientation of the field based on a degree value

        """

        if (power < 0) or (power > 1):
            raise Exception("Transmitter power must be between 0 and 1.")

        moment_min, moment_max = (0.00628318531, 0.02714336053)
        y = moment_min + power * (moment_max - moment_min)
        self.moment = np.array([0, y, 0])

        self.location = location
        self.theta = theta


    def field(self, location):
        """Calculate the magnetic field vector (in uA/m) at a given location.

                1     3 r (m * r)        m
        H(x) = ---- ( -----------  -  ------- ) A/m (Wikipedia: Magnetic Dipole)
               4*pi     ||r||^5       ||r||^3

        location (NumPy vector): A vector, np.array([x, y, z]), representing a
            location in space (e.g. the location of an avalanche receiver that
            receives this transmitter's signal).

        Returns: A NumPy vector, np.array([x, y, z]), representing the magnitude
            and direction of the transmitter's magnetic field (in uA/m) at the
            given location."""

        #print("field orientation: %d" %self.theta)

        xprime = location[0]
        yprime = location[1]

        theta = math.radians(self.theta)

        T = math.tan(theta)
        S = math.sin(theta)
        C = math.cos(theta)

        x = (xprime + T*yprime)/(T*S+C)
        y = (yprime-T*xprime)/(T*S+C)

        location = [x, y, 0]

        m = self.moment
        r = location - self.location
        r_norm = np.linalg.norm(r) #magnitude of the r vector

        f1 = 1 / (4 * scipy.pi)
        f2 = 3 * r * np.dot(m, r) / (r_norm ** 5) #np.dot = rx*mx + ry*my
        f3 = m / (r_norm ** 3)
        f = f1 * (f2 - f3)

        return f * (10 ** 6)



    def plot_field(self, xlim, ylim, n):
        """Plot the vector field of the transmitter (in the xy-plane).

        xlim (tuple): (x_min, x_max) representing the x-range of the plot.
        ylim (tuple): (y_min, y_max) representing the y-range of the plot.
        n (int): number of points between each of (x_min, x_max) and
            (y_min, y_max), inclusive.

        Displays the vector field illustration in a new window."""

        x0, x1 = xlim
        y0, y1 = ylim

        nc = n * 1j
        Y, X = np.ogrid[y0:y1:nc, x0:x1:nc]

        xfield = np.zeros((n, n))
        yfield = np.zeros((n, n))
        mag = np.zeros((n, n))
        for i in range(n):
            for j in range(n):
                location = np.array([X[0][i], Y[j][0], 0])
                (x, y, z) = self.field(location)

                xfield[j][i] = x
                yfield[j][i] = y
                mag[j][i] = np.linalg.norm(np.array([x, y, z]))

        #TODO: Make line width (or density) proportional to magnitude.
        #plt.quiver(X, Y, xfield, yfield)#, linewidth=0.5)
        plt.streamplot(X, Y, xfield, yfield)

        #plt.show()

class Point:
    def __init__(self, xvalue, yvalue):
        self.x = xvalue
        self.y = yvalue
        self.error = []

    def setError(self, e, degree):
        # for i in range(0, len(self.error)):
        #     d = self.error[i][1]
        #     if d == degree:
        #         self.error[i][0] = error
        # for e,d in self.error:
        #     if d == degree:

        #     else:
        self.error.append((e, degree))


searchx = []
searchy = []
bestGuess = []
errorList = []
class Searcher:
    """Models a searcher with an avalanche receiver."""

    def __init__(self, transmitter, location):
        """Initialize the searcher at a given location.

        transmitter: a specific transmitter, being located by the searcher.
        location: a vector, np.array([x, y, z]), representing the searcher's
            starting location."""

        self.transmitter = transmitter
        self.location = location


    def distance_from_transmitter(self):
        """Calculates the searcher's straight-line distance to the transmitter.

        Returns: The straight-line distance from the searcher's current location
            to the specified transmitter."""

        s = self.location
        t = self.transmitter.location
        r = np.linalg.norm(s - t)

        #np.linalg.norm is the distance formula
        return r


    def step(self, stepsize = 0.2):
        """Take a step in the direction of the field line.

        stepsize (float): the size of the step to take in the direction of the
            field line at the searcher's current location. Default value of 0.2.

        Effect: Modifies the searcher's current location.

        Returns: The searcher's new location."""

        (lx, ly, lz) = self.location
        (fx, fy, fz) = self.transmitter.field(self.location)
        c = stepsize / np.linalg.norm(np.array([fx, fy, fz]))
        self.location = (lx + c * fx, ly + c * fy, lz + c * fz)

        return self.location


    def search(self, radius, stepsize = 0.2, steps = 1000, visualize = False):
        """Perform a search from the searcher's current location towards the
        specified avalanche transmitter.

        radius (float): the distance away from the transmitter at which to
            terminate the search.
        stepsize (float): the size of the step to take in the direction of the
            field line at the searcher's corrent location. Default value is 0.2.
        steps (int): the maximum number of steps that can be taken before
            terminating the search.
        visualize (bool): if true, a visualization is displayed in a new window
            showing the transmitter's field and each of the searcher's
            locations.

        Effect: Modifies the searcher's current location.

        Returns: The searcher's final location at the termination of the search.
        """

        if visualize:
            t.plot_field(xlim = (-2, 2), ylim = (-1, 1), n = 100)

        for i in range(steps):
            if self.distance_from_transmitter() < radius:
                break
            if visualize:
                (x, y, z) = self.location
                searchx.append(x)
                searchy.append(y)
                #plt.annotate("x", (x, y))
            self.step(stepsize = stepsize)

        if visualize:
            plt.draw()
            #plt.show()

        return self.location

    #updates the bestGuess array to hold transmitter locations in order
    def storeGuess(self, error, p, degree):
        #append if list is too short
        errorList.append((error, p, degree))
        p.setError(error, degree)
        count = len(bestGuess)
        if count < 4:
            bestGuess.append((error, p))
            bestGuess.sort(key = lambda t: t[0])
        #else check if smaller than the largest item in the list
        elif error < bestGuess[count-1][0]:
            bestGuess.append((error, p))
            bestGuess.sort(key = lambda t: t[0])
            del bestGuess[count]

        #print("x: %.2f, y: %.2f, e: %.3f" %(p.x, p.y, error))

#######################################################################

    #Searcher class runs a more advanced version of search
    def advancedSearch(self, radius, stepsize = 0.2, steps = 1000, visualize = False):
    #Actual brute force algorithm calculation
    #at each step of the search, calculate the likelihood of a transmitter
    #at each location, constraining the tested locations by perpendicular 
    #lines
        searchx = []
        searchy = []
        r = []

        pixelX = np.linspace(x_lower_bound, x_upper_bound, 50)
        pixelY = np.linspace(y_lower_bound, y_upper_bound, 50)

        pointsList = []
        slopeList = []
        mList = []
        bList = []

        for px in pixelX:
            for py in pixelY:
                pointsList.append(Point(px, py))

        for i in range(steps):
            del bestGuess[:]
            (x, y, z) = self.location
            searchx.append(x)
            searchy.append(y)
            r.append(self.distance_from_transmitter())
        
            if (i > 0):
                #solve for the slope
                slope = (y-searchy[i-1])/(x-searchx[i-1])
                slopeList.append(slope)

                #find the perpendicular line constants
                m = -1.0/slope
                b = searchy[i-1] - m*searchx[i-1]

                #check which side the current point is on
                nexty = y
                eqn = m*x + b

                mList.append(m)
                bList.append(b)
                
                tempList = []
                for p in pointsList:
                    #check if the point is on the wrong side of the line
                    del p.error[:]
                    #if ((nexty > eqn and p.y > m*p.x + b) or 
                    #    (nexty < eqn and p.y < m*p.x + b)):
                    tempList.append(p)
                    for degree in range(0, 90, 5):
                        test = Transmitter(1, location=np.array([p.x, p.y, 0]), theta=degree)
                        error = 0
                        for k in range(1, len(slopeList)):
                            #location is the searchx and searchy
                            location = np.array([searchx[k], searchy[k], 0])
                            (testx, testy, testz) = test.field(location)
                            
                            #find the distance between test and search location k
                            testr = np.linalg.norm(test.location - location)
                            testslope = testy/testx
                            
                            #need to compare testy/testx with slope at last location
                            #using RMSE 
                            error += math.sqrt((slopeList[k] - testslope)**2 + (r[k] - testr)**2)

                        if (i > 1):
                            self.storeGuess(error, p, degree)
                            #print ("x: %.4f, y: %.4f, error: %.2f" % (p.x, p.y, error))

                pointsList = tempList
            if self.distance_from_transmitter() < radius:
                break

            self.step(stepsize = stepsize)
        
        avgx = 0
        avgy = 0        
        for e,p in bestGuess:
            avgx += p.x
            avgy += p.y
        xguess = avgx/len(bestGuess)
        yguess = avgy/len(bestGuess)
        print("Source guess: %.2f, %.2f" % (xguess,yguess))

        #draws the plot
        if visualize:
            for count in range(0, 18):
                print count*5
                t.theta = count*5
                t.plot_field(xlim = (-2, 2), ylim = (-1, 1), n = 100)
                #draw search path
                for x,y in zip(searchx, searchy):    
                    plt.annotate("x", (x, y))

                #add perpendicular lines  
                # for m,b in zip(mList, bList):
                #     xarray = np.arange(x_lower_bound, x_upper_bound, .1)
                #     yarray = np.zeros(len(xarray))
                #     yarray = m*xarray + b 
                #     plt.plot(xarray, yarray)
                #     plt.draw()

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
                for e,p in bestGuess:
                    plt.annotate("*", (p.x, p.y))
                plt.annotate("o", (xguess, yguess))

                plt.xlim([x_lower_bound, x_upper_bound])
                plt.ylim([y_lower_bound, y_upper_bound])

                plt.show()


        return self.location


####################################################################
#Initiates a generic search


# #Define a transmitter
# t = Transmitter(1, location=np.array([0, 0, 0]))

# # Define a searcher
# l = np.array([2, 0, 0])
# s = Searcher(t, l)

# # Perform a search
# s.search(.3, .2, 1000, True)


# ####################################################################
# #Analysis of search by creating transmitters at each location and 
# #determining their error, and removing points that lie on the 
# #wrong side of the perpendicular line

# x_upper_bound = 2
# x_lower_bound = -2
# y_upper_bound = 1
# y_lower_bound = -1


# searchx = searchx[1:7] #1:7 when 14 spots
# searchy = searchy[1:7]

# #for x,y in zip(searchx,searchy):
# #    print x,y

# #create vector array 
# slope = []
# for i in range(1, len(searchy)):
#     s = (searchy[i]-searchy[i-1])/(searchx[i]-searchx[i-1])
#     slope.append(s)
# slope.append(s)
# #print slope

# #Brute force search

# pixelX = np.linspace(x_lower_bound, x_upper_bound, 40)
# pixelY = np.linspace(y_lower_bound, y_upper_bound, 40)

# #minerror = 5000
# errorRank = []
# xvals = []
# yvals = []

# for i in range(0, len(pixelX)):
#     for j in range(0, len(pixelY)):

#         test = Transmitter(1, location=np.array([pixelX[i], pixelY[j], 0]))
#         error = 0
#         for k in range(1, len(searchx)):
#             location = np.array([searchx[k], searchy[k], 0])
#             (testx, testy, testz) = test.field(location)
            
#             #need to compare testy/testx with slope at that location
#             testslope = testy/testx
#             error += abs(slope[k] - testslope)

#             #to rotate the graph, compare to all the angles plus c

#         #print ("(%.2f,%.2f): %f" %(pixelX[i], pixelY[j], error))
#         #minerror = min(error, minerror)
#         xvals.append(pixelX[i])
#         yvals.append(pixelY[j])
#         errorRank.append(error)


# zipped = zip(errorRank, xvals, yvals)

# zipped.sort(key = lambda t: t[0])

# #for e,x,y in zipped:
#     #print ("(%.2f,%.2f): %f" % (x,y,e))

# #Key observations: 
# #   more points to test is better
# #   more of the path means a more accurate guess

# pixelX = np.linspace(x_lower_bound, x_upper_bound, 40)
# pixelY = np.linspace(y_lower_bound, y_upper_bound, 40)

# pointsList = []

# for px in pixelX:
#     for py in pixelY:
#         pointsList.append(Point(px, py))


# #for each point, remove pixels that are not on the same side of the
# # perpendicular line as the previous point

# for i in range(1, len(slope)-1):
#     #get perpendicular line equation, y = mx+b
#     m = -1.0/slope[i]
#     b = searchy[i] - m*searchx[i]

#     #check which side the last point is on
#     lasty = searchy[i-1]
#     eqn = m*searchx[i-1] + b

#     # xarray = np.arange(x_lower_bound, x_upper_bound, .1)
#     # yarray = np.zeros(len(xarray))
#     # yarray = m*xarray + b 
#     # plt.plot(xarray, yarray)
#     # plt.draw()

#     if (lasty > eqn):
#         #get rid of points above, eqn should be y <= mx+b
#         tempList = []
#         for p in pointsList:
#             if (p.y < m*p.x + b):
#                 tempList.append(p)
#         pointsList = tempList


#     elif (lasty < eqn):
#         #get rid of points below, eqn should be y >= mx+b
#         tempList = []
#         for p in pointsList:
#             if (p.y > m*p.x + b):
#                 tempList.append(p)
#         pointsList = tempList


# #for p in pointsList:
#     #plt.annotate(".", (p.x, p.y))



# plt.show()
# #plt.clf()
####################################################################   
#conducts the search with source-finding analysis at each step
def main():
    x_upper_bound = 2
    x_lower_bound = -2
    y_upper_bound = 1
    y_lower_bound = -1

    theta = int(random.random()*90)
    sourcex = random.random()*4 - 2
    sourcey = random.random()*2 - 1

    # Define a transmitter
    t = Transmitter(1, location=np.array([0, 0, 0]), theta=theta)

    #Determine random starting location for searcher
    startx = 2#random.random()*4 - 2 #prolbem with x: .3949 y: -.0701
    starty = 1#random.random() - 1

    print ("x: %.4f, y: %.4f, theta: %d" % (sourcex, sourcey, theta))

    # Define a searcher
    l = np.array([startx, starty, 0])
    s = Searcher(t, l)

    # Perform a search
    s.advancedSearch(.3, .1, 5, True)


    errorList.sort(key=lambda t: t[0]) 
    for e,p,d in errorList:
        print ("x: %.4f, y: %.4f, error: %.2f, degree: %d" % (p.x, p.y, e, d))



    plt.show()











