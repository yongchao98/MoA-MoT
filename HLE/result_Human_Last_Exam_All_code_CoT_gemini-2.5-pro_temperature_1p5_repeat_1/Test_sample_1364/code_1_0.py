import numpy as np
from scipy.spatial import ConvexHull

def check_projections(vertices, N):
    """
    Checks if the projections of a set of vertices onto the coordinate planes
    are all quadrilaterals. The polyhedron is the convex hull of these vertices.
    """
    print(f"Checking for a polyhedron with {N} vertices...")
    vertex_set = np.array(vertices)
    
    # We choose three orthogonal planes (xy, xz, yz), which are in a general position.
    projections = {
        "XY": {"plane_axes": [0, 1]},
        "XZ": {"plane_axes": [0, 2]},
        "YZ": {"plane_axes": [1, 2]},
    }
    
    is_possible = True
    for name, proj_info in projections.items():
        # Project the vertices onto the plane
        projected_points = vertex_set[:, proj_info["plane_axes"]]
        
        # We need to handle the case where all projected points are collinear
        if len(np.unique(projected_points, axis=0)) < 3:
             num_hull_vertices = 0 # Not a polygon
        else:
            try:
                # Compute the convex hull of the projected points
                hull = ConvexHull(projected_points)
                num_hull_vertices = len(hull.vertices)
            except: # QJ_ERR, happens with collinear points
                 num_hull_vertices = 0
        
        print(f"  Projection on {name} plane forms a polygon with {num_hull_vertices} vertices.")
        if num_hull_vertices != 4:
            is_possible = False
    
    if is_possible:
        print(f"  SUCCESS: N = {N} is a possible number of vertices.\n")
        return True
    else:
        print(f"  FAILURE: This vertex set for N = {N} does not work.\n")
        return False

# --- Main Program ---
possible_n_values = []

# N = 4: A specific tetrahedron
vertices_4 = [[0, 0, 0], [1, 1, 0], [1, 0, 1], [0, 1, 1]]
if check_projections(vertices_4, 4):
    possible_n_values.append(4)

# N = 5: A specific 5-vertex polyhedron
vertices_5 = [[0, 0, 0], [1, 0, 0], [0, 1, 0], [0, 0, 1], [1, 1, 1]]
if check_projections(vertices_5, 5):
    possible_n_values.append(5)

# N = 6: Vertices of a cube missing two opposite corners
vertices_6 = [[1, 0, 0], [0, 1, 0], [0, 0, 1], [1, 1, 0], [1, 0, 1], [0, 1, 1]]
if check_projections(vertices_6, 6):
    possible_n_values.append(6)

# N = 7: Vertices of a cube missing one corner
vertices_7 = [[0, 0, 0], [1, 0, 0], [0, 1, 0], [0, 0, 1], [1, 1, 0], [1, 0, 1], [0, 1, 1]]
if check_projections(vertices_7, 7):
    possible_n_values.append(7)

# N = 8: The 8 vertices of a cube
vertices_8 = [[0, 0, 0], [1, 0, 0], [0, 1, 0], [0, 0, 1], [1, 1, 0], [1, 0, 1], [0, 1, 1], [1, 1, 1]]
if check_projections(vertices_8, 8):
    possible_n_values.append(8)

print("----------------------------------------------------")
print("The complete set of possible numbers of vertices is:")
print(set(possible_n_values))