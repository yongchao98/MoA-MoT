import numpy as np

# Define the points from graph B
p_t = np.array([0, 12])
i_l = np.array([-4, 10])
i_r = np.array([4, 10])
i_b = np.array([0, 9])

# The triangle in question is (p_t, i_l, i_r)
# The point to check is i_b

# Function to find the circumcenter of a triangle
def find_circumcenter(p1, p2, p3):
    """Calculates the circumcenter of a triangle defined by three points."""
    D = 2 * (p1[0] * (p2[1] - p3[1]) + p2[0] * (p3[1] - p1[1]) + p3[0] * (p1[1] - p2[1]))
    
    # Check for colinearity
    if D == 0:
        return None, float('inf')

    p1_sq = p1[0]**2 + p1[1]**2
    p2_sq = p2[0]**2 + p2[1]**2
    p3_sq = p3[0]**2 + p3[1]**2
    
    ux = (p1_sq * (p2[1] - p3[1]) + p2_sq * (p3[1] - p1[1]) + p3_sq * (p1[1] - p2[1])) / D
    uy = (p1_sq * (p3[0] - p2[0]) + p2_sq * (p1[0] - p3[0]) + p3_sq * (p2[0] - p1[0])) / D
    
    center = np.array([ux, uy])
    radius_sq = (p1[0] - ux)**2 + (p1[1] - uy)**2
    
    return center, radius_sq

# Calculate circumcircle for triangle (p_t, i_l, i_r)
triangle_vertices = [p_t, i_l, i_r]
center, radius_sq = find_circumcenter(*triangle_vertices)

# Check if point i_b is inside the circumcircle
# A point is inside if its squared distance to the center is less than the squared radius
dist_sq_to_center = (i_b[0] - center[0])**2 + (i_b[1] - center[1])**2

print(f"Triangle vertices: P_T={tuple(p_t)}, I_L={tuple(i_l)}, I_R={tuple(i_r)}")
print(f"Point to check: I_B={tuple(i_b)}")
print(f"Circumcenter of the triangle: ({center[0]:.2f}, {center[1]:.2f})")
print(f"Squared radius of the circumcircle: {radius_sq:.2f}")
print(f"Squared distance of I_B from the center: {dist_sq_to_center:.2f}")

if dist_sq_to_center < radius_sq:
    print("\nResult: Point I_B is INSIDE the circumcircle of triangle (P_T, I_L, I_R).")
    print("Therefore, triangulation B violates the empty circle property.")
else:
    print("\nResult: Point I_B is OUTSIDE the circumcircle.")
    print("This specific test does not show a violation for triangulation B.")
