import numpy as np

def get_circumcircle(p1, p2, p3):
    """Calculates the circumcenter (h, k) and radius squared (r^2) of a triangle."""
    D = 2 * (p1[0] * (p2[1] - p3[1]) + p2[0] * (p3[1] - p1[1]) + p3[0] * (p1[1] - p2[1]))
    if abs(D) < 1e-10: return None, None
    
    p1_sq = p1[0]**2 + p1[1]**2
    p2_sq = p2[0]**2 + p2[1]**2
    p3_sq = p3[0]**2 + p3[1]**2
    
    h = (p1_sq * (p2[1] - p3[1]) + p2_sq * (p3[1] - p1[1]) + p3_sq * (p1[1] - p2[1])) / D
    k = (p1_sq * (p3[0] - p2[0]) + p2_sq * (p1[0] - p3[0]) + p3_sq * (p2[0] - p1[0])) / D
    
    center = np.array([h, k])
    radius_sq = (p1[0] - h)**2 + (p1[1] - k)**2
    
    return center, radius_sq

# Assign representative coordinates to the 7 points based on the image.
points = {
    'P_bl': np.array([-2.0, 0.0]),   # Bottom-left
    'P_br': np.array([2.0, 0.0]),    # Bottom-right
    'P_il': np.array([0.0, 0.8]),   # Inner-lower
    'P_iu': np.array([0.0, 2.2]),   # Inner-upper
    'P_tl': np.array([-1.5, 3.5]), # Top-left
    'P_tr': np.array([1.5, 3.5]), # Top-right
    'P_top': np.array([0.0, 4.0])     # Top
}

# In triangulation D, consider the triangle T = (P_top, P_bl, P_br).
triangle_vertices = [points['P_top'], points['P_bl'], points['P_br']]

# Calculate the circumcircle of T.
center, radius_sq = get_circumcircle(*triangle_vertices)
h, k = center[0], center[1]
r_sq = radius_sq

print("Analysis of Triangulation D")
print("Consider the triangle with vertices: P_top(0.0, 4.0), P_bl(-2.0, 0.0), P_br(2.0, 0.0)")
print(f"Its circumcircle has the equation: (x - {h:.2f})^2 + (y - {k:.2f})^2 = {r_sq:.2f}")

# Check if another point, e.g., P_il, lies inside this circle.
point_to_test = points['P_il']
px, py = point_to_test[0], point_to_test[1]

print(f"\nTesting if point P_il({px}, {py}) is inside the circle.")
dist_sq = (px - h)**2 + (py - k)**2

# The condition for a point (px, py) to be inside is (px - h)^2 + (py - k)^2 < r^2.
print(f"We check the inequality: ({px:.2f} - {h:.2f})^2 + ({py:.2f} - {k:.2f})^2 < {r_sq:.2f}")
print(f"Plugging in the numbers: ({px - h:.2f})^2 + ({py - k:.2f})^2 < {r_sq:.2f}")
print(f"Calculating the terms: {((px - h)**2):.2f} + {((py - k)**2):.2f} < {r_sq:.2f}")
print(f"Final check: {dist_sq:.2f} < {r_sq:.2f}")

if dist_sq < r_sq:
    print("\nThe inequality is TRUE.")
    print("Conclusion: The point P_il is inside the circumcircle of triangle (P_top, P_bl, P_br).")
    print("Therefore, Triangulation D violates the empty circle property.")
else:
    print("\nThe inequality is FALSE.")
