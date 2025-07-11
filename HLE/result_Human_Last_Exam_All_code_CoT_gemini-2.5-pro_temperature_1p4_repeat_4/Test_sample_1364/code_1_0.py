import numpy as np
from scipy.spatial import ConvexHull

def count_silhouette_vertices(vertices_3d, proj_direction):
    """
    Calculates the number of silhouette vertices of a polyhedron when projected along a given direction.
    
    Args:
        vertices_3d (np.ndarray): A NumPy array of shape (n, 3) representing the polyhedron's vertices.
        proj_direction (np.ndarray): A NumPy array of shape (3,) representing the projection direction.
        
    Returns:
        int: The number of vertices in the 2D projection (silhouette).
    """
    # 1. Normalize the projection direction vector.
    d = np.asarray(proj_direction) / np.linalg.norm(proj_direction)

    # 2. Create an orthonormal basis for the 2D projection plane.
    # Find a vector not parallel to d.
    u = np.array([1, 0, 0])
    if np.allclose(np.abs(np.dot(d, u)), 1.0):
        u = np.array([0, 1, 0])
        
    basis_v1 = np.cross(d, u)
    basis_v1 /= np.linalg.norm(basis_v1)
    basis_v2 = np.cross(d, basis_v1)
    basis_v2 /= np.linalg.norm(basis_v2)

    # 3. Project the 3D vertices onto the 2D plane.
    vertices_2d = np.dot(vertices_3d, np.stack([basis_v1, basis_v2], axis=1))

    # 4. Compute the 2D convex hull of the projected points.
    hull = ConvexHull(vertices_2d)

    # 5. The number of vertices of the hull is the number of silhouette vertices.
    return len(hull.vertices)

def check_polyhedron(name, vertices, directions):
    """
    Checks if a polyhedron meets the problem's criteria and prints the results.
    """
    print(f"--- Checking: {name} (V={len(vertices)}) ---")
    is_solution = True
    for i, d in enumerate(directions):
        try:
            k = count_silhouette_vertices(np.array(vertices), d)
            print(f"  Projection {i+1} direction {d}: silhouette has {k} vertices.")
            if k != 4:
                is_solution = False
        except Exception as e:
            print(f"  An error occurred during calculation for direction {d}: {e}")
            is_solution = False
            
    if is_solution:
        print(f"Result: {name} is a valid solution.\n")
    else:
        print(f"Result: {name} is NOT a valid solution.\n")

# Define vertices for candidate polyhedra
# 1. Tetrahedron (V=4)
tetra_vertices = [
    [1, 1, 1], [1, -1, -1], [-1, 1, -1], [-1, -1, 1]
]
tetra_directions = [
    [1, 0, 0], [0, 1, 0], [0, 0, 1]
]
check_polyhedron("Tetrahedron", tetra_vertices, tetra_directions)

# 2. Square Pyramid (V=5)
sq_pyramid_vertices = [
    [1, 1, 0], [1, -1, 0], [-1, -1, 0], [-1, 1, 0], [0, 0, 1]
]
# Directions from the open cone "looking at the base"
sq_pyramid_directions = [
    [0, 0, -1], [0.1, 0, -1], [0, 0.1, -1]
]
check_polyhedron("Square Pyramid", sq_pyramid_vertices, sq_pyramid_directions)

# 3. Octahedron (V=6)
octa_vertices = [
    [1, 0, 0], [-1, 0, 0], [0, 1, 0], [0, -1, 0], [0, 0, 1], [0, 0, -1]
]
octa_directions = [
    [1, 0, 0], [0, 1, 0], [0, 0, 1]
]
check_polyhedron("Octahedron", octa_vertices, octa_directions)

# 4. Cube (V=8)
cube_vertices = [
    [1, 1, 1], [1, 1, -1], [1, -1, 1], [1, -1, -1],
    [-1, 1, 1], [-1, 1, -1], [-1, -1, 1], [-1, -1, -1]
]
cube_directions = [
    [1, 0, 0], [0, 1, 0], [0, 0, 1]
]
check_polyhedron("Cube", cube_vertices, cube_directions)

# The final answer is the set of these possible vertex counts.
final_set = {4, 5, 6, 8}
print(f"The set of possible numbers of vertices is: {final_set}")