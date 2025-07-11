import numpy as np
from scipy.spatial import ConvexHull

def get_projection_vertices_count(points, direction):
    """
    Projects 3D points onto a plane orthogonal to the given direction
    and returns the number of vertices in their 2D convex hull.
    """
    # Normalize the direction vector
    d = np.array(direction) / np.linalg.norm(direction)

    # Create an orthonormal basis for the projection plane
    # Find a vector 't' not parallel to 'd'
    t = np.array([1.0, 0.0, 0.0])
    if np.allclose(np.cross(d, t), 0):
        t = np.array([0.0, 1.0, 0.0])
    
    u = np.cross(d, t)
    u /= np.linalg.norm(u)
    v = np.cross(d, u)

    # Project the 3D points to 2D
    projected_points = []
    for p in points:
        projected_points.append([np.dot(p, u), np.dot(p, v)])
    
    # Compute the convex hull of the 2D points
    hull = ConvexHull(np.array(projected_points))
    
    # Return the number of vertices of the convex hull
    return len(hull.vertices)

def check_polyhedron(name, vertices, directions):
    """
    Checks if a polyhedron satisfies the condition for a given set of directions.
    """
    num_vertices = len(vertices)
    print(f"Checking V={num_vertices} ({name})")
    is_solution = True
    for i, d in enumerate(directions):
        try:
            proj_v_count = get_projection_vertices_count(vertices, d)
            print(f"  Projection onto plane with normal {np.round(d, 2).tolist()} has {proj_v_count} vertices.")
            if proj_v_count != 4:
                is_solution = False
        except Exception as e:
            print(f"  Could not compute projection {i+1}: {e}")
            is_solution = False

    if is_solution:
        print(f"CONCLUSION: V={num_vertices} is a possible number of vertices.\n")
    else:
        print(f"CONCLUSION: V={num_vertices} is NOT a possible number of vertices for this configuration.\n")

def main():
    # Case V=5: Trigonal Bipyramid
    # It works for directions slightly perturbed from its equatorial plane.
    h = 1.0  # height of poles
    r = 1.0  # radius of equator
    eps = 0.1 # perturbation
    trigonal_bipyramid_v = [
        [0, 0, h], [0, 0, -h], # Poles
        [r, 0, 0], # Equatorial vertices
        [r * np.cos(2*np.pi/3), r * np.sin(2*np.pi/3), 0],
        [r * np.cos(4*np.pi/3), r * np.sin(4*np.pi/3), 0]
    ]
    # Three linearly independent directions
    trigonal_bipyramid_d = [
        [1, 0, eps],
        [np.cos(2*np.pi/3), np.sin(2*np.pi/3), eps],
        [np.cos(4*np.pi/3), np.sin(4*np.pi/3), eps]
    ]
    check_polyhedron("Trigonal Bipyramid", trigonal_bipyramid_v, trigonal_bipyramid_d)
    
    # Case V=6: Octahedron
    # It works for projections along the coordinate axes.
    octahedron_v = [
        [1, 0, 0], [-1, 0, 0],
        [0, 1, 0], [0, -1, 0],
        [0, 0, 1], [0, 0, -1]
    ]
    octahedron_d = [
        [1, 0, 0],
        [0, 1, 0],
        [0, 0, 1]
    ]
    check_polyhedron("Octahedron", octahedron_v, octahedron_d)

    # Case V=7: Capped Octahedron
    # We add a small cap on one face. The projections along the same axes should still be quadrilaterals.
    # We choose a cap height 'c' such that 1/3 < c <= 1/2. Let c = 0.4
    capped_octahedron_v = octahedron_v + [[0.4, 0.4, 0.4]]
    check_polyhedron("Singly Capped Octahedron", capped_octahedron_v, octahedron_d)
    
    print("Based on these examples and the face-capping argument (explained in the text),")
    print("we can generate a valid polyhedron for any number of vertices V >= 5.")
    print("The set of possible numbers of vertices is {5, 6, 7, ...}.")


if __name__ == '__main__':
    main()
