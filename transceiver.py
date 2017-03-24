# use (/) for floating point division and (//) for integer division
from __future__ import division

import numpy as np
import scipy as scipy
import matplotlib.pyplot as plt


class Transmitter:
    """Models an avalanche transmitter at the origin."""

    def __init__(self, power):
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

        """

        if (power < 0) or (power > 1):
            raise Exception("Transmitter power must be between 0 and 1.")

        moment_min, moment_max = (0.00628318531, 0.02714336053)
        y = moment_min + power * (moment_max - moment_min)
        self.moment = np.array([0, y, 0])

        self.location = np.array([0, 0, 0])


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
        plt.quiver(X, Y, xfield, yfield)#, linewidth=0.5)
        
        #plt.show()


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

        #assuming np.linalg.norm is the distance formula

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

        # for i in range(steps):
        #     if self.distance_from_transmitter() < radius:
        #         break
        #     if visualize:
        #         (x, y, z) = self.location
        #         plt.annotate("x", (x, y))
        #     self.step(stepsize = stepsize)

        if visualize:
            plt.draw()
            plt.show()

        return self.location


# Define a transmitter
t = Transmitter(1)

# Define a searcher
l = np.array([2, 0, 0])
s = Searcher(t, l)

# Perform a search
s.search(.3, .2, 1000, True)
