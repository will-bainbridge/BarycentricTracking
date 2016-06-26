#!/usr/bin/python3

import numpy as np

def mag(x, axis=None):
    return np.sqrt(np.sum(np.square(x), axis=axis))

def position(t, x):
    fT = np.sin(2*np.pi*t)
    fX = np.expand_dims(np.sin(np.pi*x[:,0]/6), axis=-1)
    c, s = np.cos(0.1*2*np.pi*fT), np.sin(0.1*2*np.pi*fT)
    R = np.array([[1, 0, 0], [0, c, -s], [0, s, c]])
    o = np.array([3, 0.35416667, -0.08333333])
    temp = (x - o).dot(R) + o
    return fX*((x - o).dot(R) + o) + (1 - fX)*x
