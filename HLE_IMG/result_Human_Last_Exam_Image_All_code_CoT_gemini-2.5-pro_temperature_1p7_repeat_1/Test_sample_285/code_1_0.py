import numpy as np

def get_circumcircle(p1, p2, p3):
    """
    Calculates the center and radius squared of the circumcircle of a triangle.
    The points are given as (x, y) tuples.
    This uses the formula derived from the intersection of perpendicular bisectors.
    """
    D = 2 * (p1[0] * (p2[1] - p3[1]) + p2[0] * (p3[1] - p1[1]) + p3[0] * (p1[1] - p2[1]))
    if abs(D) < 1e-10: # Collinear points
        return None, None

    # Calculate center coordinates
    ux = ((p1[0]**2 + p1[1]**2) * (p2[1] - p3[1]) + (p2[0]**2 + p2[1]**2) * (p3[1] - p1[1]) + (p3[0]**2 + p3[1]**2) * (p1[1] - p2[1])) / D
    uy = ((p1[0]**2 + p1[1]**2) * (p3[0] - p2[0]) + (p2[0]**2 + p2[1]**2) * (p1[0] - p3[0]) + (p3[0]**2 + p3[1]**2) * (p2[0] - p1[0])) / D
    
    center = (ux, uy)
    
    # Calculate radius squared
    radius_sq = (p1[0] - ux)**2 + (p1[1] - uy)**2
    
    return center, radius_sq

# --- Main analysis ---

# Define coordinates for the points in triangulation B based on the image.
# Exact values are not critical, only their relative positions.
P_t = (5.0, 9.5)     # Top vertex
Q_ibl = (3.0, 4.0)   # Inner bottom-left vertex
Q_ibr = (7.0, 4.0)   # Inner bottom-right vertex
Q_c = (5.0, 5.0)     # Central inner vertex

# The triangle T that appears to violate the property
violating_triangle_points = [P_t, Q_ibl, Q_ibr]
# The point P that appears to be inside T's circumcircle
violating_point = Q_c

print("Analyzing Triangulation B for Delaunay violation.")
print("-" * 50)
print(f"Let the triangle be T with vertices A={violating_triangle_points[0]}, B={violating_triangle_points[1]}, C={violating_triangle_points[2]}.")
print(f"Let the point to check be P={violating_point}.")
print("\nWe check if P lies inside the circumcircle of T.")

# Calculate the circumcircle
center, radius_sq = get_circumcircle(violating_triangle_points[0], violating_triangle_points[1], violating_triangle_points[2])

print("\nStep 1: Calculate the circumcircle of triangle T.")
print(f"The equation for the circumcenter (cx, cy) is found by ensuring it's equidistant from A, B, and C.")
print(f"The calculated circumcenter is ({center[0]:.4f}, {center[1]:.4f}).")
print(f"The equation for the radius squared (R^2) is the squared distance from the center to any vertex.")
print(f"R^2 = ({violating_triangle_points[0][0]:.1f} - {center[0]:.4f})^2 + ({violating_triangle_points[0][1]:.1f} - {center[1]:.4f})^2 = {radius_sq:.4f}")

# Calculate distance of the violating point to the center
dist_sq_to_center = (violating_point[0] - center[0])**2 + (violating_point[1] - center[1])**2

print("\nStep 2: Check if point P is inside the circle.")
print(f"The squared distance from P to the center is d^2.")
print(f"d^2 = ({violating_point[0]:.1f} - {center[0]:.4f})^2 + ({violating_point[1]:.1f} - {center[1]:.4f})^2 = {dist_sq_to_center:.4f}")

# Compare and conclude
print("\nStep 3: Compare d^2 with R^2.")
print(f"Is {dist_sq_to_center:.4f} < {radius_sq:.4f}?")
if dist_sq_to_center < radius_sq:
    print("Yes, it is. The point P lies inside the circumcircle of triangle T.")
    print("\nConclusion: Triangulation B violates the empty circle property.")
else:
    print("No, it is not. The property holds for this triangle and point.")
