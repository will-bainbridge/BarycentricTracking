#!/usr/bin/python3

# This file describes a mesh which looks like this...
# 
#   (3)o------o------o------o
#     / \    / \    / \    / \
# (2)o------o------o------o   \
#    |    \ |    \ |    \ |    \
#    | (4)o-|----o-|----o-|----o
# (1)o------o------o------o   /
#     \ /    \ /    \ /    \ /
#   (0)o------o------o------o

# Points
def translate(points, displacement):
    return [ [ x + d for x, d in zip(point, displacement) ] for point in points ]
pointsAtX0 = [
        [ 0,  0, -2],
        [ 0, -1, -1],
        [ 0, -1,  1],
        [ 0,  1,  2],
        [ 0,  2, -1]]
points = (
        pointsAtX0 +
        translate(pointsAtX0, [3, 0, 0]) +
        translate(pointsAtX0, [3, 0, 0]) +
        translate(pointsAtX0, [6, 0, 0]))

# Faces
def offset(faces, offset):
    return [ [ p + offset for p in face ] for face in faces ]
facesNormalToXAtX0 = [
        [0, 1, 2, 3, 4]]
facesNotNormalToXAtX0 = [
        [0, 5, 6, 1],
        [1, 6, 7, 2],
        [2, 7, 8, 3],
        [3, 8, 9, 4],
        [4, 9, 5, 0]]
faces = (
        offset(facesNormalToXAtX0, 5) + 
        offset(facesNormalToXAtX0, 10) + 
        facesNormalToXAtX0 +
        offset([ facesNormalToXAtX0[0][::-1] ], 15) + 
        facesNotNormalToXAtX0 + 
        offset(facesNotNormalToXAtX0, 5) + 
        offset(facesNotNormalToXAtX0, 10))

# Cells (indirectly)
faceOwners = [1, 2, 0, 2] + [0]*5 + [1]*5 + [2]*5
faceNeighbours = [0, 1]
