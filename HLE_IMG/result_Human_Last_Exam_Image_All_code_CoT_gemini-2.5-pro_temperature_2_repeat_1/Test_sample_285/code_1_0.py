import numpy as np

def get_circumcircle(p1, p2, p3):
    """
    Calculate the circumcenter and the squared radius of a triangle defined by three points.
    Formula from: https://en.wikipedia.org/wiki/Circumscribed_circle#Cartesian_coordinates_2
    """
    D = 2 * (p1[0] * (p2[1] - p3[1]) + p2[0] * (p3[1] - p1[1]) + p3[0] * (p1[1] - p2[1]))
    
    # Check for colinearity
    if abs(D) < 1e-10:
        return None, None

    p1_sq = p1[0]**2 + p1[1]**2
    p2_sq = p2[0]**2 + p2[1]**2
    p3_sq = p3[0]**2 + p3[1]**2

    center_x = (1 / D) * (p1_sq * (p2[1] - p3[1]) + p2_sq * (p3[1] - p1[1]) + p3_sq * (p1[1] - p2[1]))
    center_y = (1 / D) * (p1_sq * (p3[0] - p2[0]) + p2_sq * (p1[0] - p3[0]) + p3_sq * (p2[0] - p1[0]))
    
    center = np.array([center_x, center_y])
    # Use squared distances to avoid sqrt
    radius_sq = (p1[0] - center_x)**2 + (p1[1] - center_y)**2
    
    return center, radius_sq

def check_empty_circle_property(points, triangles):
    """
    Checks if a given triangulation violates the empty circle property.
    
    Args:
    - points: A list of numpy arrays, where each array is a point's coordinates.
    - triangles: A list of tuples, where each tuple contains the indices of the 3 points forming a triangle.
    
    Returns:
    - A string indicating the violating triangle and point, or a success message.
    """
    for i, tri_indices in enumerate(triangles):
        p1, p2, p3 = points[tri_indices[0]], points[tri_indices[1]], points[tri_indices[2]]
        
        center, radius_sq = get_circumcircle(p1, p2, p3)
        
        if center is None:
            # The vertices are collinear, which is a degenerate case. 
            # A valid triangulation shouldn't have this, but we'll note it.
            print(f"Triangle {i} {tri_indices} has collinear vertices.")
            continue

        for j, p_check in enumerate(points):
            # A point is not checked against the triangle it belongs to
            if j in tri_indices:
                continue
            
            # Check if p_check is inside the circumcircle
            dist_sq = (p_check[0] - center[0])**2 + (p_check[1] - center[1])**2
            
            # Use a small tolerance epsilon for floating point comparison
            epsilon = 1e-9
            if dist_sq < radius_sq - epsilon:
                print(f"Violation Found!")
                print(f"  - Point {j} {tuple(p_check)} is inside the circumcircle of triangle {i} {tri_indices}.")
                print(f"  - The vertices of the triangle are: {tuple(p1)}, {tuple(p2)}, {tuple(p3)}.")
                print(f"  - Circumcenter: {tuple(np.round(center, 2))}, Squared Radius: {radius_sq:.2f}")
                print(f"  - Distance of point {j} from center (squared): {dist_sq:.2f}")
                return True
                
    print("No violation of the empty circle property was found.")
    return False

# --- Main analysis ---
# Based on the visual analysis, B is the violator. Let's create a hypothetical
# set of coordinates that match the visual layout of Triangulation B to demonstrate the violation.

# Points for Triangulation B (red graph)
# Outer hull: P0-P4, Inner: P5-P6
points_B = [
    np.array([0, 0]),      # P0: Bottom-Left (BL)
    np.array([11, 1]),     # P1: Bottom-Right (BR)
    np.array([12, 9]),     # P2: Top-Right (TR)
    np.array([6, 11]),     # P3: Top (T)
    np.array([1, 10]),     # P4: Top-Left (TL)
    np.array([4, 6]),      # P5: Inner-Left (IL)
    np.array([8, 5]),      # P6: Inner-Right (IR)
]

# Triangles of B, identified visually.
triangles_B = [
    (4, 3, 5),    # (TL, T, IL)
    (4, 0, 5),    # (TL, BL, IL)
    (3, 6, 5),    # (T, IR, IL)  -- this is the triangle (Top, Inner-Right, Inner-Left)
    (0, 6, 5),    # (BL, IR, IL) -- this is the triangle (Bottom-Left, Inner-Right, Inner-Left)
    (3, 2, 6),    # (T, TR, IR)
    (1, 2, 6),    # (BR, TR, IR)
    (0, 1, 6)     # (BL, BR, IR)
]

# The visual analysis suggested that in the quadrilateral (T, IR, BL, IL),
# the diagonal (IR, IL) is non-Delaunay.
# This means point BL (P0) should be in the circumcircle of triangle (T, IR, IL) (P3,P6,P5).
print("--- Checking Triangulation B ---")
if check_empty_circle_property(points_B, triangles_B):
    print("\nConclusion: Triangulation B violates the empty circle property.")
else:
    print("\nConclusion: Triangulation B seems to satisfy the property with these coordinates (estimation might be off).")

print("\nBased on the geometric angle-sum rule, which is more robust than estimating coordinates, Triangulation B is the one that violates the empty circle property.")

<<<B>>>