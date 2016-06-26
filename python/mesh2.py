#!/usr/bin/python3

# This file describes a mesh which looks like this...
# 
# (2)              (3)
#   o--------------o
#   |\             |
#   | \            |
#   |  \           |
#   |   \          |
#   |    \ (4)     |
#   |     o-__     |
#   |         --__ |
#   o--------------o
# (0)              (1)

# Points
def translate(points, displacement):
    return [ [ x + d for x, d in zip(point, displacement) ] for point in points ]
pointsAtZ0 = [
        [0   ,0   ,0   ],
        [1   ,0   ,0   ],
        [0   ,1   ,0   ],
        [1   ,1   ,0   ],
        [0.25,0.25,0   ]]
points = (
        pointsAtZ0 +
        translate(pointsAtZ0, [0, 0, 1]))

# Faces
faces = [
        [1,4,9,6],
        [4,2,7,9],
        [0,2,4,1],
        [1,4,2,3],
        [5,6,9,7],
        [6,8,7,9],
        [0,1,6,5],
        [1,3,8,6],
        [3,2,7,8],
        [0,5,7,2]]

# Cells (indirectly)
faceOwners = [0, 0, 0, 1, 0, 1, 0, 1, 1, 0]
faceNeighbours = [1, 1]
