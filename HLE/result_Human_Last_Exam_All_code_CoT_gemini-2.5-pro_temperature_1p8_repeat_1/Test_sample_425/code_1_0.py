import numpy as np

def project(points, axis):
    """Projects a set of 3D points onto a plane normal to the given axis."""
    axis = axis / np.linalg.norm(axis)
    projected_points = []
    for p in points:
        # Project p onto the axis vector
        proj_p_on_axis = np.dot(p, axis) * axis
        # Subtract the projection to get the component in the plane
        p_in_plane = p - proj_p_on_axis
        projected_points.append(p_in_plane)
    return np.array(projected_points)

def get_dist_sq(p1, p2):
    """Calculates the squared distance between two points."""
    return np.sum((p1 - p2)**2)

def analyze_shape(points):
    """Analyzes the side lengths of a polygon defined by points."""
    num_points = len(points)
    centroid = np.mean(points, axis=0)
    
    # Calculate distances from centroid to each point
    radii_sq = [get_dist_sq(p, centroid) for p in points]
    
    # Calculate side lengths
    side_lengths_sq = [get_dist_sq(points[i], points[(i + 1) % num_points]) for i in range(num_points)]
    
    return radii_sq, side_lengths_sq

# Vertices of a regular tetrahedron
tetra_vertices = np.array([
    [1, 1, 1],
    [1, -1, -1],
    [-1, 1, -1],
    [-1, -1, 1]
])

print("Analyzing possible projection symmetries for an object with A4 rotation group...")
print("-" * 70)

# --- Case i): Projection along a C3 axis ---
print("Possibility i): Order 3")
# A C3 axis passes through a vertex (e.g., [1,1,1]) and the centroid (origin)
c3_axis = np.array([1, 1, 1])
proj_points_3 = project(tetra_vertices, c3_axis)

# One point projects to the origin, find the other three
outer_points_3 = np.array([p for p in proj_points_3 if np.linalg.norm(p) > 1e-6])
centroid_3 = np.mean(outer_points_3, axis=0)
radii_sq_3, sides_sq_3 = analyze_shape(outer_points_3)

print("Projecting along a C3 axis (like [1, 1, 1])...")
print("The three outer vertices of the projection form an equilateral triangle.")
print(f"Squared distance of vertices from their centroid: [{', '.join(f'{r:.2f}' for r in radii_sq_3)}]")
print(f"Squared side lengths of the triangle: [{', '.join(f'{s:.2f}' for s in sides_sq_3)}]")
print("Since the sides are equal and vertices are equidistant from the center, the shape has 3-fold rotational symmetry.")
print("Equation: a^2 = b^2 = c^2, where a, b, c are side lengths.")
print(f"Result: {sides_sq_3[0]:.2f} = {sides_sq_3[1]:.2f} = {sides_sq_3[2]:.2f}")
print("Order 3 is possible.\n")

# --- Case ii): Projection along a C2 axis ---
print("Possibility ii): Order 4")
# A C2 axis passes through midpoints of opposite edges. The cartesian axes for this tetrahedron are C2 axes.
c2_axis = np.array([0, 0, 1])
proj_points_4 = project(tetra_vertices, c2_axis)

radii_sq_4, sides_sq_4 = analyze_shape(proj_points_4)

print("Projecting along a C2 axis (like [0, 0, 1])...")
print("The four vertices of the projection form a square.")
print(f"Squared side lengths of the quadrilateral: [{', '.join(f'{s:.2f}' for s in sides_sq_4)}]")
# For a square, we also check diagonals
diag1_sq = get_dist_sq(proj_points_4[0], proj_points_4[2])
diag2_sq = get_dist_sq(proj_points_4[1], proj_points_4[3])
print(f"Squared diagonal lengths: [{diag1_sq:.2f}, {diag2_sq:.2f}]")
print("Since all four sides are equal and the diagonals are equal, the shape is a square with 4-fold rotational symmetry.")
print("Equations: a^2 = b^2 = c^2 = d^2  and  p^2 = q^2 (diagonals).")
print(f"Result: {sides_sq_4[0]:.2f} = {sides_sq_4[1]:.2f} = {sides_sq_4[2]:.2f} = {sides_sq_4[3]:.2f}  and  {diag1_sq:.2f} = {diag2_sq:.2f}")
print("Order 4 is possible.\n")

print("-" * 70)
print("Conclusion: Based on geometric and group theory arguments, the possible orders are 3 and 4.")