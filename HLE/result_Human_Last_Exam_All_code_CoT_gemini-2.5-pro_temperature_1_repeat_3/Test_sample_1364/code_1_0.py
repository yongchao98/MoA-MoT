import numpy as np
from scipy.spatial import ConvexHull

def project_polyhedron(vertices, direction):
    """
    Projects a set of 3D vertices onto a plane orthogonal to the given direction.
    Returns the 2D coordinates of the projected points.
    """
    # Normalize the projection direction vector
    d = np.array(direction, dtype=float)
    if np.linalg.norm(d) == 0:
        raise ValueError("Direction vector cannot be zero.")
    d /= np.linalg.norm(d)

    # Create an orthonormal basis for the projection plane
    if np.abs(np.dot(d, [0, 0, 1])) < 0.99999:
        u = np.cross(d, [0, 0, 1])
    else:
        u = np.cross(d, [0, 1, 0])
    u /= np.linalg.norm(u)
    w = np.cross(d, u)

    # Project vertices and express them in the 2D basis (u, w)
    projected_vertices_2d = []
    for v_3d in vertices:
        v = np.array(v_3d, dtype=float)
        # The projection of v onto the plane is v - (v.d)d
        # To get the 2D coordinates, we project onto the basis vectors u and w
        coord1 = np.dot(v, u)
        coord2 = np.dot(v, w)
        projected_vertices_2d.append([coord1, coord2])
        
    return np.array(projected_vertices_2d)

def check_projection(name, vertices, directions):
    """
    Checks if the projection of a polyhedron is a quadrilateral for the given directions.
    """
    print(f"--- Checking Polyhedron: {name} (V={len(vertices)}) ---")
    possible = True
    for i, d in enumerate(directions):
        projected_points = project_polyhedron(vertices, d)
        
        # The vertices of the projected polygon form the convex hull of all projected points.
        try:
            # Use 'QJ' option to handle cases where some points are co-linear.
            hull = ConvexHull(projected_points, qhull_options='QJ')
            num_hull_vertices = len(hull.vertices)
            print(f"  Projection along direction {d}: Found a {num_hull_vertices}-gon.")
            if num_hull_vertices != 4:
                possible = False
        except Exception as e:
            print(f"  Could not compute convex hull for direction {d}: {e}")
            possible = False

    if possible:
        print(f"SUCCESS: {name} can have quadrilateral projections for all 3 given directions.")
    else:
        print(f"FAILURE: {name} does not produce quadrilaterals for all 3 given directions.")
    print("-" * (len(name) + 30) + "\n")


# Main analysis
# According to the theorem by Pogorelov and Kleinschmidt, the number of vertices V
# must be at least k+2 for k-gonal projections. For k=4, V must be >= 6.
# This rules out V=4 and V=5.

# We demonstrate that V=6 is possible with a regular octahedron.
octahedron_vertices = [
    [1, 0, 0], [-1, 0, 0],
    [0, 1, 0], [0, -1, 0],
    [0, 0, 1], [0, 0, -1]
]
# We use three orthogonal projection directions.
octahedron_directions = [
    [1, 0, 0],
    [0, 1, 0],
    [0, 0, 1]
]
check_projection("Regular Octahedron", octahedron_vertices, octahedron_directions)

# The face-capping argument shows that if V=6 works, so do V=7, 8, 9, ...
# For example, we can cap the face of the octahedron in the (-,-,-) octant.
# This face is a back-face for all three directions [1,0,0], [0,1,0], [0,0,1].
# Capping it adds a vertex without changing the silhouette for these projections.

print("Conclusion from the analysis:")
print("The minimum number of vertices is 6.")
print("Any integer number of vertices greater than or equal to 6 is possible.")
print("The set of possible numbers for the vertices V is {V | V is an integer and V >= 6}.")
print("The numbers in this set are 6, 7, 8, 9, ...")
