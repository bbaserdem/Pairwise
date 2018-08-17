# Pairwise

I'm trying to run a custom SPH simulation.
SPH is a fluid dynamic simulation where particles do short-range interactions
Thus, I have this specific case, where;

* I need to calculate pairwise distances below a certain threshold.
* Needs to generalize easily to arbitrary dimensions. (2D to 3D)
* Needs to run fast when threshold is low and particle density is large.
* Needs to be able to do this on a torus. (Simplify boundary conditions)

This description of problem requires either wasted computational power,
or building of an algorithm.

For this to work, I have come up with the algorithm of the next section.
Implementation wise, the particles will be moving around.
The tradeoff with this algorithm is rapid-insertion data structures are needed,
which MATLAB does not have.
So this neccessitates the building of a specialized linked-list.

This problem refers to the following;

* N: Number of points
* V: Volume of simulation
* D: Dimensionality of the system

# Algorithm

While a bulk all-to-all approach can work pretty well, it is a horrible idea.
In the case where the threshold distance is a small fraction of the geometry,
there will be a lot of operations wasted.
Whatmore, if there are a lot of points, either vectorization has to be given up
or allocation has to go up to N^2.
There will be double calculation of distances in an all-to-all approach;
double lookup of both points are uneccassary tradeoff.

To avoid this problem; I propose the following data structure.

* The points live in subspaces, to be referred as voxels.
Points belong in each of these voxels.
Each point's displacement is to be referred from a voxel origin.
* Each voxel is a D-dimensional, and referred by integer indices.

To find nearest neighbours;

* Points within individual voxels are candidates for being below threshold,
and can have their pairwise distances calculated.
This is `O((NR^D/V)^2)`
* For each voxel, compare each point and the points of a neighbour voxel.
Pixel coordinates need to be modulated with the displacement vector.
This also takes care of symmetric borders.
Each neighbour pair needs to be considered only once.
Check neighbours on a voxel with Lﯢ norm of 1, only if displacement vector is
positive. (First non-zero component in ordered vector is positive.)
This is `O( N^2 ((3R)^D/V)^2 )`, much better than `O(N^2)`.

## Positive vectors

For N-dimensions, positive vectors of Lﯢ norm can be found with the iteration;

* Grab all the positive vectors in N-1 dimensions, and set N'th index to 0.
* Grab all 1 Lﯢ norm N-1 vectors, and set N'th index to 1.
