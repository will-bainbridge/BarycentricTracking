# BarycentricTracking
This is a minimal proof of concept of three-dimensional tracking in barycentric coordinates.
Tracking is the process of following a particle as it progresses through a mesh.
This algorithm tracks a particle though a polyhedral mesh.
The mesh is decomposed into tetrahedra, and a particle is at any time associated with the tetrahedron which contains it.
The particle position is stored in barycentric coordinates, local to its tetrahedron.

Traditional tracking algorithms calculate the intersection between particle and tetrahedron face in global coordinates.
The problem with this approach is that the relative floating point error associated with calculating the location of the intersection can be very large.
If the error means the intersection cannot be found, the particle is effectively lost.

This tracking algorithm is formulated purely in terms of these barycentric coordinates.
Floating point error occurs as a result of converting the desired track direction to these coordinates.
Calculating the location of the intersection, however, is now exact, as it occurs when one local coordinate equals zero.
There is no danger of being unable to find an intersection or losing the particle.
In addition, the algorithm is capable of tracking through inverted or degenerate tetrahedra.

# Use
For an example of tracking through a stationary mesh, run...

    ./track.py 0.05 mesh1 path1
    
For an example of tracking through a moving mesh (warning, work in progress), run...

    ./track.py 0.05 mesh1 path1 motion1
