#!/usr/bin/python3

import numpy as np

centre = np.array([3, 0.05, 0.4])
xAxis = np.array([2.2, -0.2, -0.1])
yAxis = np.cross(xAxis, [0, 0.2, -0.2])

def position(t):
    return centre + np.cos(2*np.pi*t)*xAxis + np.sin(2*np.pi*t)*yAxis
