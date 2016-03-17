#!/usr/bin/python3

#------------------------------------------------------------------------------#

'''
Copyright (C) 2016 Will Bainbridge

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program.  If not, see <http://www.gnu.org/licenses/>.
'''

#------------------------------------------------------------------------------#

import numpy as np
import sys

np.set_printoptions(suppress=True)
np.random.seed(0)

#------------------------------------------------------------------------------#

def mag(x, axis=None):
    return np.sqrt(np.sum(np.square(x), axis=axis))

def flat(x):
    return [ item for sub in x for item in sub ]

#------------------------------------------------------------------------------#

class Mesh(object):

    def __init__(self, name):
        mesh = __import__(name)

        self.points = np.array(mesh.points, dtype=float)

        self.faces = np.array(mesh.faces)
        self.faceOwners = np.array(mesh.faceOwners)
        self.faceNeighbours = np.array(mesh.faceNeighbours)

        self.faceTetPoints = np.array([ np.random.randint(len(face)) for face in self.faces ])

        self.cells = [ [] for i in range(max(self.faceOwners) + 1) ]
        for f, c in enumerate(self.faceOwners):
            self.cells[c].append(f)
        for f, c in enumerate(self.faceNeighbours):
            self.cells[c].append(f)
        self.cells = np.array(self.cells)

        self.pointsOldest = self.points
        self.update()
        self.movePointsWithoutUpdate(self.points)

    def getFaceAreasAndCentres(self):
        faceAreas = np.zeros((len(self.faces), 3))
        faceCentres = np.zeros((len(self.faces), 3))
        for f in range(len(self.faces)):
            tri1 = self.points[self.faces[f]]
            tri2 = self.points[np.append(self.faces[f][1:], self.faces[f][0])]
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
        return faceAreas, faceCentres

    def getCellVolumesAndCentres(self):
        cellVolumes = np.zeros((len(self.cells)))
        temp = np.sum(self.faceAreas*self.faceCentres, axis=1)/3
        np.add.at(cellVolumes, self.faceOwners, temp)
        np.subtract.at(cellVolumes, self.faceNeighbours, temp[:len(self.faceNeighbours)])
        cellCentres = np.zeros((len(self.cells), 3))
        mask = cellVolumes == 0
        temp = np.expand_dims(np.sum(self.faceAreas*self.faceCentres, axis=1), axis=1)*self.faceCentres/4
        np.add.at(cellCentres, self.faceOwners, temp)
        np.subtract.at(cellCentres, self.faceNeighbours, temp[:len(self.faceNeighbours)])
        if np.any(~mask):
            cellCentres[~mask] /= np.expand_dims(cellVolumes[~mask], axis=1)
        temp = np.expand_dims(mag(self.faceAreas, axis=1), axis=1)
        if np.any(mask):
            cellCentres[mask] = [ np.sum((temp*self.faceCentres)[cell], axis=0)/np.sum(temp[cell]) for cell in self.cells[mask] ]
        return cellVolumes, cellCentres

    def getIntermediatePoints(self, l):
        return l*self.points + (1 - l)*self.pointsOld

    def update(self):
        self.faceAreas, self.faceCentres = self.getFaceAreasAndCentres()
        self.cellVolumes, self.cellCentres = self.getCellVolumesAndCentres()

    def movePointsWithoutUpdate(self, newPoints):
        self.pointsOld = self.points
        self.points = newPoints
        self.faceAreasOld, self.faceCentresOld = self.faceAreas, self.faceCentres
        self.cellVolumesOld, self.cellCentresOld = self.cellVolumes, self.cellCentres

    def movePoints(self, newPoints):
        self.movePointsWithoutUpdate(newPoints)
        self.update()

#------------------------------------------------------------------------------#

timeStep = float(sys.argv[1])
mesh = Mesh(sys.argv[2])
path = __import__(sys.argv[3])
motion = __import__(sys.argv[4]) if len(sys.argv) > 4 else None

#------------------------------------------------------------------------------#

def faceBehind(f, i, n):
    return (i + len(mesh.faces[f]) - n) % len(mesh.faces[f])

def faceAhead(f, i, n):
    return (i + n) % len(mesh.faces[f])

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
    c = mesh.faceOwners[f] if isOwnerTet(t) else mesh.faceNeighbours[f]
    iB = mesh.faceTetPoints[f]
    i1 = faceAhead(f, iB, ownerTet(t))
    i2 = faceAhead(f, iB, ownerTet(t) + 1)
    pB = mesh.faces[f][iB]
    p1 = mesh.faces[f][i1 if isOwnerTet(t) else i2]
    p2 = mesh.faces[f][i2 if isOwnerTet(t) else i1]
    return c, pB, p1, p2

def getTetGeometry(f, t):
    c, pB, p1, p2 = getTetTopology(f, t)
    vO = mesh.cellCentres[c]
    vB = mesh.points[pB]
    v1 = mesh.points[p1]
    v2 = mesh.points[p2]
    return vO, vB, v1, v2

def getTetTransform(f, t):
    vO, vB, v1, v2 = getTetGeometry(f, t)
    dB = vB - vO
    d1 = v1 - vO
    d2 = v2 - vO
    A = np.array([dB, d1, d2]).T
    return vO, A

def getTetReverseTransform(f, t):
    vO, A = getTetTransform(f, t)
    dB = A[:,0]
    d1 = A[:,1]
    d2 = A[:,2]
    detA = dB.dot(np.cross(d1, d2))
    T = np.array([
            np.cross(d1, d2),
            np.cross(d2, dB),
            np.cross(dB, d1),
            ])
    return vO, detA, T

def getMovingTetGeometry(f, t, l):
    c, pB, p1, p2 = getTetTopology(f, t)
    vO = np.array([mesh.cellCentresOld[c], mesh.cellCentres[c]])
    vB = np.array([mesh.pointsOld[pB], mesh.points[pB]])
    v1 = np.array([mesh.pointsOld[p1], mesh.points[p1]])
    v2 = np.array([mesh.pointsOld[p2], mesh.points[p2]])
    vO = np.array([vO[0] + l*(vO[1] - vO[0]), (1 - l)*(vO[1] - vO[0])])
    vB = np.array([vB[0] + l*(vB[1] - vB[0]), (1 - l)*(vB[1] - vB[0])])
    v1 = np.array([v1[0] + l*(v1[1] - v1[0]), (1 - l)*(v1[1] - v1[0])])
    v2 = np.array([v2[0] + l*(v2[1] - v2[0]), (1 - l)*(v2[1] - v2[0])])
    return vO, vB, v1, v2

def getMovingTetTransform(f, t, l):
    vO, vB, v1, v2 = getMovingTetGeometry(f, t, l)
    dB = vB - vO
    d1 = v1 - vO
    d2 = v2 - vO
    A0 = np.array([dB[0], d1[0], d2[0]]).T
    A1 = np.array([dB[1], d1[1], d2[1]]).T
    return vO, np.array([A0, A1])

def getMovingTetReverseTransform(f, t, l):
    vO, A = getMovingTetTransform(f, t, l)
    dB = A[:,:,0]
    d1 = A[:,:,1]
    d2 = A[:,:,2]
    detA0 = dB[0].dot(np.cross(d1[0], d2[0]))
    detA1 = dB[1].dot(np.cross(d1[0], d2[0])) + dB[0].dot(np.cross(d1[1], d2[0])) + dB[0].dot(np.cross(d1[0], d2[1]))
    detA2 = dB[0].dot(np.cross(d1[1], d2[1])) + dB[1].dot(np.cross(d1[0], d2[1])) + dB[1].dot(np.cross(d1[1], d2[0]))
    detA3 = dB[1].dot(np.cross(d1[1], d2[1]))
    T0 = np.array([
        np.cross(d1[0], d2[0]),
        np.cross(d2[0], dB[0]),
        np.cross(dB[0], d1[0]),
        ])
    T1 = np.array([
        np.cross(d1[0], d2[1]) + np.cross(d1[1], d2[0]),
        np.cross(d2[0], dB[1]) + np.cross(d2[1], dB[0]),
        np.cross(dB[0], d1[1]) + np.cross(dB[1], d1[0]),
        ])
    T2 = np.array([
        np.cross(d1[1], d2[1]),
        np.cross(d2[1], dB[1]),
        np.cross(dB[1], d1[1]),
        ])
    return vO, np.array([detA0, detA1, detA2, detA3]), np.array([T0, T1, T2])

def toBarycentricPosition(x):
    assert x.shape == (3, )
    return np.append(1 - np.sum(x), x)

def toBarycentricDisplacement(x):
    assert x.shape == (3, )
    return np.append(- np.sum(x), x)

def fromBarycentricPosition(x):
    assert x.shape == (4, )
    # <-- Assumes barycentrics function like homogeneous coordinates.
    return (x/np.sum(x))[1:]

#def fromBarycentricDisplacement(x):
#    assert x.shape == (4, )
#    # <-- Is this sufficient? How could it be normalised? The sum is zero for a displacement.
#    return x[1:]

def fromTetCoordinates(f, t, y, l):
    origin, A = getTetTransform(f, t)
    return A.dot(fromBarycentricPosition(y)) + origin

def toTetCoordinates(f, t, x, l):
    origin, detA, T = getTetReverseTransform(f, t)
    if detA == 0:
        return np.ones(4)*np.nan
    return toBarycentricPosition(T.dot(x - origin)/detA)

def fromMovingTetCoordinates(f, t, y, l):
    origin, A = getMovingTetTransform(f, t, l)
    return A[0].dot(fromBarycentricPosition(y)) + origin[0]

def toMovingTetCoordinates(f, t, x, l):
    origin, detA, T = getTetMovingTransform(f, t)
    if detA[0] == 0:
        return np.ones(4)*np.nan
    return toBarycentricPosition(T[0].dot(x - origin[0])/detA[0])

def findTet(x):
    for f in range(len(mesh.faces)):
        for t in range(1, len(mesh.faces[f]) - 1):
            y = toTetCoordinates(f, t, x, 0)
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

def updatePointsPlot(plot, l):
    updatePlot(plot, mesh.getIntermediatePoints(l))

def updateFacesPlot(plot, l, shrink=0.8, centres=None):
    if centres is None:
        centres = mesh.faceCentres
    points = mesh.getIntermediatePoints(l)
    temp = points[flat(mesh.faces)]
    lines = np.array([
        temp,
        points[flat( np.roll(face, 1) for face in mesh.faces )],
        np.nan*np.ones(temp.shape)
        ]).swapaxes(0, 1)
    origins = np.expand_dims(flat(
        np.ones((len(face), 1))*centre
        for face, centre in zip(mesh.faces, centres) ), axis=1)
    temp = np.reshape(shrink*(lines - origins) + origins, (3*lines.shape[0], 3))
    updatePlot(plot, temp)

def updateCellsPlot(plot, l, shrink=0.8, centres=None):
    if centres is None:
        centres = mesh.cellCentres
    points = mesh.getIntermediatePoints(l)
    temp = points[flat(flat( mesh.faces[cell] for cell in mesh.cells ))]
    lines = np.array([
        temp,
        points[flat( np.roll(face, 1) for face in flat( mesh.faces[cell] for cell in mesh.cells ))],
        np.nan*np.ones(temp.shape)
        ]).swapaxes(0, 1)
    origins = np.expand_dims(flat(
        np.ones((len(flat(mesh.faces[cell])), 1))*centre
        for cell, centre in zip(mesh.cells, centres) ), axis=1)
    temp = np.reshape(shrink*(lines - origins) + origins, (3*lines.shape[0], 3))
    updatePlot(plot, temp)

def updateTetPlot(plot, f, t, l, shrink=0.8, centre=None):
    vO, vB, v1, v2 = [ v[0] for v in getMovingTetGeometry(f, t, l) ]
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
    isOwner = mesh.faceOwners[f] == c
    # Find the face with the shared edge
    if i == 1:
        edge = [p1, p2]
    elif i == 2:
        edge = [p2, pB]
    elif i == 3:
        edge = [pB, p1]
    else:
        assert False
    temp = [ g for g in mesh.cells[c] if g != f and set(edge).issubset(mesh.faces[g]) ]
    assert len(temp) == 1
    g = temp[0]
    # Get the tet point on the new face
    isNewOwner = mesh.faceOwners[g] == c
    u = list(mesh.faces[g]).index(edge[isNewOwner]) - mesh.faceTetPoints[g]
    u = (u + len(mesh.faces[g])) % len(mesh.faces[g])
    u = np.clip(u, 1, len(mesh.faces[g]) - 2)
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
            if t == len(mesh.faces[f]) - 2:
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
            if ownerTet(t) == len(mesh.faces[f]) - 2:
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
    origin, detA, T = getTetReverseTransform(f, t)
    # Get the local displacement
    TX1 = toBarycentricDisplacement(T.dot((1 - l)*x1))
    # Calculate the hit fraction
    iH = -1
    muH = np.finfo(float).max if detA == 0 else np.abs(1/detA)
    for i in range(4):
        if TX1[i] < 0:
            mu = - y0[i]/TX1[i]
            if 0 <= mu and mu < muH:
                iH = i
                muH = mu
    # Set the new y
    yH = y0 + muH*TX1
    # Remove tolerance issues in the event of a hit
    if iH != -1:
        yH[iH] = 0
    # Return the hit index, the new position, and the new tracking parameter
    return yH, iH, l + (1 - l)*muH*detA

def trackThroughMovingTet(f, t, y0, x1, l):
    '''
    This performs the track in the same way as the function above, but for a
    moving tet.
    '''
    # Get the tet geometry
    origin, detA, T = getMovingTetReverseTransform(f, t, l)
    # Form the determinant and hit polynomials
    x0 = fromMovingTetCoordinates(f, t, y0, l)
    x0Rel = x0 - origin[0]
    x1Rel = (1 - l)*(x1 - origin[1])
    detAPoly = np.array([1, detA[1], detA[2]*detA[0], detA[3]*detA[0]**2])
    hitPolyB12 = np.array([
        y0[1:],
        T[0].dot(x1Rel) + T[1].dot(x0Rel),
        (T[1].dot(x1Rel) + T[2].dot(x0Rel))*detA[0],
        (T[2].dot(x1Rel))*detA[0]**2,
        ]).T
    hitPolyO = np.expand_dims(detAPoly - np.sum(hitPolyB12, axis=0), axis=0)
    hitPoly = np.append(hitPolyO, hitPolyB12, axis=0)
    dHitPoly = hitPoly[:,1:]*np.arange(1, 4)
    # Calculate the hit fraction
    iH = -1
    muH = np.finfo(float).max if detA[0] == 0 else np.abs(1/detA[0])
    for i in range(4):
        # <-- This root finding step is likely to have to deal with unpleasant cases
        #     Cases where the high order coefficients are zero or near-zero are difficult to handle
        #     Numpy does fine
        muRoots = np.roots(hitPoly[i][::-1])
        # <-- This method compares the imaginary part to zero, so this wouldn't work in the event of a small imaginary error
        #     Numpy doesn't seem to generate such errors, though, so this is OK
        muRoots = np.real(muRoots[np.isreal(muRoots)])
        dHitPolyVal = dHitPoly[i][0] + muRoots*dHitPoly[i][1] + muRoots*muRoots*dHitPoly[i][2]
        muRoots = muRoots[(dHitPolyVal < 0) & (0 <= muRoots) & (muRoots < muH)]
        if muRoots.size:
            iH = i
            muH = np.min(muRoots)
    # Set the new y
    detAPolyVal = detAPoly[0] + muH*detAPoly[1] + muH*muH*detAPoly[2] + muH*muH*muH*detAPoly[3]
    hitPolyVal = hitPoly[:,0] + muH*hitPoly[:,1] + muH*muH*hitPoly[:,2] + muH*muH*muH*hitPoly[:,3]
    # <-- This calculation will need special handling for cases when both polynomials approach zero
    #     Testing for zero will be accomplished by comparing the values to an estimate of the truncation error
    #     If both polynomials are zero, the limit will be evaluated by differentiating
    #     If the denominator is zero, but the numerator is not, an error will result
    yH = hitPolyVal/detAPolyVal
    # <-- In the event where both polynomials approach zero, the calculated yH will not be on a face
    #     A second track step will be required, through the degenerate tet
    #     This can be of the simpler, static mesh form
    # Remove tolerance issues in the event of a hit
    if iH != -1:
        yH[iH] = 0
    # Return the hit index, the new position, and the new tracking parameter
    return yH, iH, l + (1 - l)*muH*detA[0]

#------------------------------------------------------------------------------#

import matplotlib.pyplot as pp
import matplotlib.animation as an
from mpl_toolkits.mplot3d import Axes3D

# Set up the figure and scale the axes
fg = pp.figure()
ax = fg.add_subplot(111, projection='3d')
minPoints = np.min(mesh.points, axis=0)
maxPoints = np.max(mesh.points, axis=0)
ax.set_xlim3d(minPoints[0], maxPoints[0])
ax.set_ylim3d(minPoints[1], maxPoints[1])
ax.set_zlim3d(minPoints[2], maxPoints[2])

# Create the plots
plots = [createPointsPlot(ax), createFacesPlot(ax), createCellsPlot(ax), createTetPlot(ax), ax.plot([], [], [], 'bo--')[0]]

def track(data, *plots):
    # If the track is finished, set up the next one
    if track.i == -1:
        track.time += timeStep
        if motion is not None:
            mesh.movePoints(motion.position(track.time, mesh.pointsOldest))
        track.l = 0
        track.x0 = (fromTetCoordinates if motion is None else fromMovingTetCoordinates)(track.f, track.t, track.y, track.l)
        track.x1 = path.position(track.time) - track.x0
    # If on a face, change tet and update the plot
    else:
        track.f, track.t, track.y = changeTet(track.f, track.t, track.i, track.y)
    # Move through the tet and plot the new position
    track.y, track.i, track.l = (trackThroughTet if motion is None else trackThroughMovingTet)(track.f, track.t, track.y, track.x1, track.l)
    if np.any(0 > track.y) or np.any(track.y > 1):
        raise NameError('The track is outside the tet')
    track.path.append((fromTetCoordinates if motion is None else fromMovingTetCoordinates)(track.f, track.t, track.y, track.l))
    track.path = track.path[max(0, len(track.path) - 10):]
    # Update the plots
    updatePointsPlot(plots[0], track.l)
    updateFacesPlot(plots[1], track.l)
    updateCellsPlot(plots[2], track.l)
    updatePlot(plots[4], np.array(track.path))
    updateTetPlot(plots[3], track.f, track.t, track.l, shrink=1.0)

# Initialise the track
track.time = 0
track.i = -1
track.f, track.t, track.y = findTet(path.position(0))
track.path = [path.position(0)]
track.hit = False

# Initialise the plots
updatePointsPlot(plots[0], 0)
updateFacesPlot(plots[1], 0)
updateCellsPlot(plots[2], 0)
updatePlot(plots[4], np.array(track.path))
updateTetPlot(plots[3], track.f, track.t, 0, shrink=1.0)

# Animate
animation = an.FuncAnimation(fg, track, fargs=(plots), interval=100)
pp.show()

## Draw
#for i in range(15):
#    track(None, *plots)
#temp = np.array([track.x0, track.x0 + track.x1])
#ax.plot(temp[:,0], temp[:,1], temp[:,2], 'kx--')
#pp.show()
