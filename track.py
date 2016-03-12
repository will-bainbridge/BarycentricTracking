#!/usr/bin/python3

import numpy as np
import sys

np.set_printoptions(suppress=True)

#------------------------------------------------------------------------------#

def mag(x, axis=None):
    return np.sqrt(np.sum(np.square(x), axis=axis))

def flat(x):
    return [ item for sub in x for item in sub ]

#------------------------------------------------------------------------------#

mesh = __import__(sys.argv[1])

points = np.array(mesh.points, dtype=float)
faces = np.array(mesh.faces)
faceOwners = np.array(mesh.faceOwners)
faceNeighbours = np.array(mesh.faceNeighbours)
faceTetPoints = np.array([ np.random.randint(len(face)) for face in faces ])

cells = [ [] for i in range(max(faceOwners) + 1) ]
for f, c in enumerate(faceOwners):
    cells[c].append(f)
for f, c in enumerate(faceNeighbours):
    cells[c].append(f)
cells = np.array(cells)

faceAreas = np.zeros((len(faces), 3))
faceCentres = np.zeros((len(faces), 3))
for f in range(len(faces)):
    tri1 = points[faces[f]]
    tri2 = points[np.append(faces[f][1:], faces[f][0])]
    base = np.mean(tri1, axis=0)
    triAreas = np.cross(tri1 - base, tri2 - base)
    triCentres = tri1/3 + tri2/3 + base/3
    faceAreas[f] = np.sum(triAreas, axis=0)/2
    m = mag(faceAreas[f])
    if m != 0:
        n = faceAreas[f]/m
        faceCentres[f] = np.sum(triCentres*np.expand_dims(triAreas.dot(n), axis=1), axis=0)/2
        faceCentres[f] /= faceAreas[f].dot(n)
    else:
        faceCentres[f] = np.mean(triCentres, axis=0)

cellVolumes = np.zeros((len(cells)))
temp = np.sum(faceAreas*faceCentres, axis=1)/3
np.add.at(cellVolumes, faceOwners, temp)
np.subtract.at(cellVolumes, faceNeighbours, temp[:len(faceNeighbours)])

cellCentres = np.zeros((len(cells), 3))
mask = cellVolumes == 0
temp = np.expand_dims(np.sum(faceAreas*faceCentres, axis=1), axis=1)*faceCentres/4
np.add.at(cellCentres, faceOwners, temp)
np.subtract.at(cellCentres, faceNeighbours, temp[:len(faceNeighbours)])
if np.any(~mask):
    cellCentres[~mask] /= np.expand_dims(cellVolumes[~mask], axis=1)
temp = np.expand_dims(mag(faceAreas, axis=1), axis=1)
if np.any(mask):
    cellCentres[mask] = [ np.sum((temp*faceCentres)[cell], axis=0)/np.sum(temp[cell]) for cell in cells[mask] ]

#------------------------------------------------------------------------------#

def faceBehind(f, i, n):
    return (i + len(faces[f]) - n) % len(faces[f])

def faceAhead(f, i, n):
    return (i + n) % len(faces[f])

def isOwnerTet(t):
    return t >= 0

def isNeighbourTet(t):
    return not isOwnerTet(t)

def ownerTet(t):
    return t if isOwnerTet(t) else - 1 - t

def neighbourTet(t):
    return t if isNeighbourTet(t) else - 1 - t

def otherTet(t):
    return neighbourTet(t) if isOwnerTet(t) else ownerTet(t)

def getTetTopology(f, t):
    c = faceOwners[f] if isOwnerTet(t) else faceNeighbours[f]
    iB = faceTetPoints[f]
    i1 = faceAhead(f, iB, ownerTet(t))
    i2 = faceAhead(f, iB, ownerTet(t) + 1)
    pB = faces[f][iB]
    p1 = faces[f][i1 if isOwnerTet(t) else i2]
    p2 = faces[f][i2 if isOwnerTet(t) else i1]
    return c, pB, p1, p2

def getTetGeometry(f, t):
    c, pB, p1, p2 = getTetTopology(f, t)
    vO = cellCentres[c]
    vB = points[pB]
    v1 = points[p1]
    v2 = points[p2]
    return vO, vB, v1, v2

def getTetTransform(f, t):
    vO, vB, v1, v2 = getTetGeometry(f, t)
    return vO, np.array([vB - vO, v1 - vO, v2 - vO]).T

def getTetReverseTransform(f, t):
    vO, vB, v1, v2 = getTetGeometry(f, t)
    xB = np.cross(v1 - vO, v2 - vO)
    x1 = np.cross(v2 - vO, vB - vO)
    x2 = np.cross(vB - vO, v1 - vO)
    return vO, (vB - vO).dot(xB), np.array([xB, x1, x2])

def toBarycentricPosition(x):
    assert x.shape == (3, )
    return np.append(1 - np.sum(x), x)

def toBarycentricDisplacement(x):
    assert x.shape == (3, )
    return np.append(- np.sum(x), x)

def fromBarycentricPosition(x):
    assert x.shape == (4, )
    return (x/np.sum(x))[1:]

def fromBarycentricDisplacement(x):
    assert x.shape == (4, )
    return x[1:] # <-- ?

def toTetCoordinates(f, t, x):
    origin, detA, T = getTetReverseTransform(f, t)
    if detA == 0:
        return np.ones(4)*np.nan
    return toBarycentricPosition(T.dot(x - origin)/detA)

def fromTetCoordinates(f, t, x):
    origin, A = getTetTransform(f, t)
    return A.dot(fromBarycentricPosition(x)) + origin

def findTet(x):
    for f in range(len(faces)):
        for t in range(1, len(faces[f]) - 1):
            y = toTetCoordinates(f, t, x)
            if np.min(y) >= 0 and np.max(y) <= 1:
                return f, t, y
    assert False

#------------------------------------------------------------------------------#

def createPointsPlot(plotter, style='ko'):
    return plotter.plot([], [], [], style)[0]

def createFacesPlot(plotter, style='g-'):
    return plotter.plot([], [], [], style)[0]

def createCellsPlot(plotter, style='r-'):
    return plotter.plot([], [], [], style)[0]

def createTetPlot(plotter, style='m-'):
    return plotter.plot([], [], [], style)[0]

def updatePlot(plot, data):
    plot.set_data(data[:,0], data[:,1])
    plot.set_3d_properties(data[:,2])

def updatePointsPlot(plot):
    updatePlot(plot, points)

def updateFacesPlot(plot, shrink=0.8, centres=None):
    if centres is None:
        centres = faceCentres
    temp = points[flat(faces)]
    lines = np.array([
        temp,
        points[flat( np.roll(face, 1) for face in faces )],
        np.nan*np.ones(temp.shape)
        ]).swapaxes(0, 1)
    origins = np.expand_dims(flat(
        np.ones((len(face), 1))*centre
        for face, centre in zip(faces, centres) ), axis=1)
    temp = np.reshape(shrink*(lines - origins) + origins, (3*lines.shape[0], 3))
    updatePlot(plot, temp)

def updateCellsPlot(plot, shrink=0.8, centres=None):
    if centres is None:
        centres = cellCentres
    temp = points[flat(flat( faces[cell] for cell in cells ))]
    lines = np.array([
        temp,
        points[flat( np.roll(face, 1) for face in flat( faces[cell] for cell in cells ))],
        np.nan*np.ones(temp.shape)
        ]).swapaxes(0, 1)
    origins = np.expand_dims(flat(
        np.ones((len(flat(faces[cell])), 1))*centre
        for cell, centre in zip(cells, centres) ), axis=1)
    temp = np.reshape(shrink*(lines - origins) + origins, (3*lines.shape[0], 3))
    updatePlot(plot, temp)

def updateTetPlot(plot, f, t, shrink=0.8, centre=None):
    vO, vB, v1, v2 = getTetGeometry(f, t)
    vN = np.nan*np.ones(3)
    if centre is None:
        centre = (vO + vB + v1 + v2)/4
    lines = np.array([[vO, vB, vN], [vO, v1, vN], [vO, v2, vN], [vB, v1, vN], [vB, v2, vN], [v1, v2, vN]])
    temp = np.reshape(shrink*(lines - centre) + centre, (3*lines.shape[0], 3))
    updatePlot(plot, temp)

#------------------------------------------------------------------------------#

def changeCell(f, t, i, y0):
    '''
    Changing cell is quite simple. The tets match on either side of a face, so
    it's simply a matter of selecting the opposite tet index and reversing the
    flipped coordinates.
    '''
    return f, otherTet(t), y0[[0,1,3,2]]

def changeFace(f, t, i, y0):
    '''
    Changing face is the biggest pain. This code searches for the face in the
    same cell which shares the same edge. A new tet is selected using this edge.
    A strange transformation of the barycentric coordinates is then required
    whereby the coordinates are first rotated to place the base coordinate (1)
    opposite the edge. After this, the coordinates on the edge (2 and 3) are
    flipped, as the edge order is reversed on the other face. Finally, the
    coordinates are rotated so as to place the coordinate opposite the edge (1)
    at the base point of the new face.
    '''
    # Get the tet topology
    c, pB, p1, p2 = getTetTopology(f, t)
    isOwner = faceOwners[f] == c
    # Find the face with the shared edge
    if i == 1:
        edge = [p1, p2]
    elif i == 2:
        edge = [p2, pB]
    elif i == 3:
        edge = [pB, p1]
    else:
        assert False
    temp = [ g for g in cells[c] if g != f and set(edge).issubset(faces[g]) ]
    assert len(temp) == 1
    g = temp[0]
    # Get the tet point on the new face
    isNewOwner = faceOwners[g] == c
    u = list(faces[g]).index(edge[isNewOwner]) - faceTetPoints[g]
    u = (u + len(faces[g])) % len(faces[g])
    u = np.clip(u, 1, len(faces[g]) - 2)
    u = ownerTet(u) if isNewOwner else neighbourTet(u)
    # Get the new tet topology
    c, qB, q1, q2 = getTetTopology(g, u)
    # Transform the coordinates to the new tet
    preRotate = 0*(pB not in edge) + 1*(p1 not in edge) + 2*(p2 not in edge)
    postRotate = 0*(qB not in edge) + 1*(q1 not in edge) + 2*(q2 not in edge)
    temp = y0[1:]
    temp = np.roll(temp, - preRotate)
    temp = temp[[0, 2, 1]]
    temp = np.roll(temp, postRotate)
    return g, u, np.append(y0[0], temp)

def changeTet(f, t, i, y0):
    '''
    This is a mess of logic, but is relatively straighforward, at least in
    comparison with changing face. If the index is 0, we change cell. If it is
    1 we change face. If it is 2 or 3 we change tet, unless we are at the edge
    of the fan, in which case we change face.
    '''
    if i == 0:
        return changeCell(f, t, i, y0)
    elif i == 1:
        return changeFace(f, t, i, y0)
    elif i == 2:
        if isOwnerTet(t):
            if t == len(faces[f]) - 2:
                return changeFace(f, t, i, y0)
            else:
                return f, t + 1, y0[[0,1,3,2]]
        else:
            if ownerTet(t) == 1:
                return changeFace(f, t, i, y0)
            else:
                return f, neighbourTet(ownerTet(t) - 1), y0[[0,1,3,2]]
    elif i == 3:
        if isOwnerTet(t):
            if t == 1:
                return changeFace(f, t, i, y0)
            else:
                return f, t - 1, y0[[0,1,3,2]]
        else:
            if ownerTet(t) == len(faces[f]) - 2:
                return changeFace(f, t, i, y0)
            else:
                return f, neighbourTet(ownerTet(t) + 1), y0[[0,1,3,2]]
    assert False

def trackThroughTet(f, t, y0, x1, l):
    '''
    This performs the track through a tet, in the specified direction, either to
    the bounds of the tet, or until the end of the track. There are separate
    notes for how this works.
    '''
    # Get the tet geometry
    o, detA, T = getTetReverseTransform(f, t)
    # Get the local displacement
    detAY1 = toBarycentricDisplacement(T.dot((1 - l)*x1))
    iH = -1
    lHByDetA = np.finfo(float).max if detA == 0 else np.abs(1/detA)
    for i in range(4):
        if detAY1[i] != 0:
            lByDetA = - y0[i]/detAY1[i]
            if 0 < lByDetA and lByDetA < lHByDetA:
                iH = i
                lHByDetA = lByDetA
    # Set the new y
    yH = y0 + lHByDetA*detAY1
    # Remove tolerance issues in the event of a hit
    if iH != -1:
        yH[iH] = 0
    # Return the hit index, he new position, and the new tracking parameter
    return yH, iH, l + (1 - l)*lHByDetA*detA

#------------------------------------------------------------------------------#

import matplotlib.pyplot as pp
import matplotlib.animation as an
from mpl_toolkits.mplot3d import Axes3D

# Set up the figure and scale the axes
fg = pp.figure()
ax = fg.add_subplot(111, projection='3d')
minPoints = np.min(points, axis=0)
maxPoints = np.max(points, axis=0)
ax.set_xlim3d(minPoints[0], maxPoints[0])
ax.set_ylim3d(minPoints[1], maxPoints[1])
ax.set_zlim3d(minPoints[2], maxPoints[2])

# Create the plots
plots = [createPointsPlot(ax), createFacesPlot(ax), createCellsPlot(ax), createTetPlot(ax), ax.plot([], [], [], 'bo--')[0]]


#------------------------------------------------------------------------------#

## Update the mesh plots. These would need moving into the draw function if the mesh were moving.
#updatePointsPlot(plots[0])
#updateFacesPlot(plots[1])
#updateCellsPlot(plots[2])
#
## Define the frame drawing function
#def draw(data, *plots):
#    # If in the interior of a tet, then reset the track
#    if draw.i == -1:
#        draw.f, draw.t, draw.y = findTet(draw.x0)
#        draw.l = 0
#        draw.path = [draw.x0]
#        draw.moving = True
#        updateTetPlot(plots[3], draw.f, draw.t, shrink=1.0)
#    # If we are moving, then track through the tet and update the path plot
#    if draw.moving == True:
#        draw.y, draw.i, draw.l = trackThroughTet(draw.f, draw.t, draw.y, draw.x1, draw.l)
#        draw.path.append(fromTetCoordinates(draw.f, draw.t, draw.y))
#        updatePlot(plots[4], np.array(draw.path))
#    # If we are not moving, then transform to the next tet and update the tet plot
#    else:
#        draw.f, draw.t, draw.y = changeTet(draw.f, draw.t, draw.i, draw.y)
#        updateTetPlot(plots[3], draw.f, draw.t, shrink=1.0)
#    # Toggle the state
#    draw.moving = not draw.moving
#
## Initialise the track
#draw.x0 = np.array([ float(x) for x in sys.argv[2].split(',') ])
#draw.x1 = np.array([ float(x) for x in sys.argv[3].split(',') ])
#draw.i = -1
#
### Draw
##while True:
##    draw(None, *plots)
##    if draw.i == -1:
##        break
##pp.show()
#
## Animate
#animation = an.FuncAnimation(fg, draw, fargs=(plots), interval=500)
#pp.show()

#------------------------------------------------------------------------------#

path = __import__(sys.argv[2])
timeStep = float(sys.argv[3])

#motion = __import__('motion1')
#points0 = points.copy()

def track(data, *plots):
    # If the track is finished, set up the next one
    if track.i == -1:
        track.time += timeStep
        track.x0 = fromTetCoordinates(track.f, track.t, track.y)
        track.x1 = path.position(track.time) - track.x0
        track.l = 0
    # If on a face, change tet and update the plot
    else:
        track.f, track.t, track.y = changeTet(track.f, track.t, track.i, track.y)
        updateTetPlot(plots[3], track.f, track.t, shrink=1.0)
    # Move through the tet and plot the new position
    track.y, track.i, track.l = trackThroughTet(track.f, track.t, track.y, track.x1, track.l)
    track.path.append(fromTetCoordinates(track.f, track.t, track.y))
    updatePlot(plots[4], np.array(track.path))

    #track.time += timeStep
    #points[:] = motion.position(track.time, points0)
    #updatePointsPlot(plots[0])
    #updateFacesPlot(plots[1])
    #updateCellsPlot(plots[2])

# Initialise the track
track.time = 0
track.i = -1
track.f, track.t, track.y = findTet(path.position(0))
track.path = [path.position(0)]
track.hit = False

# Initialise the plots
updatePointsPlot(plots[0])
updateFacesPlot(plots[1])
updateCellsPlot(plots[2])
updatePlot(plots[4], np.array(track.path))
updateTetPlot(plots[3], track.f, track.t, shrink=1.0)

# Animate
animation = an.FuncAnimation(fg, track, fargs=(plots), interval=100)
pp.show()
