import numpy as np

def find_circle_from_three_points(p1, p2, p3):
    """
    Finds the center and radius of a circle passing through three points.
    Returns (h, k, r) for circle (x-h)^2 + (y-k)^2 = r^2
    """
    D = 2 * (p1[0] * (p2[1] - p3[1]) + p2[0] * (p3[1] - p1[1]) + p3[0] * (p1[1] - p2[1]))
    if abs(D) < 1e-9:
        # Points are collinear, cannot form a circle
        return None, None, None
        
    p1_sq = p1[0]**2 + p1[1]**2
    p2_sq = p2[0]**2 + p2[1]**2
    p3_sq = p3[0]**2 + p3[1]**2
    
    h = (p1_sq * (p2[1] - p3[1]) + p2_sq * (p3[1] - p1[1]) + p3_sq * (p1[1] - p2[1])) / D
    k = (p1_sq * (p3[0] - p2[0]) + p2_sq * (p1[0] - p3[0]) + p3_sq * (p2[0] - p1[0])) / D
    
    r = np.sqrt((p1[0] - h)**2 + (p1[1] - k)**2)
    return h, k, r

def find_line_circle_intersections(m, c, h, k, r):
    """
    Finds intersections of line y = mx + c and circle (x-h)^2 + (y-k)^2 = r^2
    """
    # (x-h)^2 + (mx+c-k)^2 = r^2
    # x^2 - 2hx + h^2 + (mx + (c-k))^2 = r^2
    # x^2 - 2hx + h^2 + m^2x^2 + 2m(c-k)x + (c-k)^2 - r^2 = 0
    # (1+m^2)x^2 + (2m(c-k)-2h)x + (h^2+(c-k)^2-r^2) = 0
    A = 1 + m**2
    B = 2 * m * (c - k) - 2 * h
    C = h**2 + (c - k)**2 - r**2
    
    delta = B**2 - 4 * A * C
    if delta < 0:
        return []
        
    x1 = (-B + np.sqrt(delta)) / (2 * A)
    y1 = m * x1 + c
    x2 = (-B - np.sqrt(delta)) / (2 * A)
    y2 = m * x2 + c
    return [np.array([x1, y1]), np.array([x2, y2])]
    
def find_intersections_with_horizontal_line(y0, h, k, r, M_x):
    """
    Finds intersection of circle with y=y0, returning the point that isn't M.
    """
    # (x-h)^2 + (y0-k)^2 = r^2
    # (x-h)^2 = r^2 - (y0-k)^2
    val = r**2 - (y0-k)**2
    if val < 0:
        return None
    
    x_offset = np.sqrt(val)
    x1 = h + x_offset
    x2 = h - x_offset

    # We return the point that is not M
    if abs(x1 - M_x) < 1e-9:
        return np.array([x2, y0])
    else:
        return np.array([x1, y0])

# --- Main logic ---

# 1. Define Circle O
O_h, O_k, O_r = 1, 2, 10
print(f"Circle O is defined by center ({O_h}, {O_k}) and radius {O_r}")

# 2. Define a horizontal chord AB on Circle O
y_AB = 8  # must be between O_k - O_r and O_k + O_r
# Find x coordinates for A and B
x_offset_AB = np.sqrt(O_r**2 - (y_AB - O_k)**2)
A = np.array([O_h - x_offset_AB, y_AB])
B = np.array([O_h + x_offset_AB, y_AB])
print(f"Chord AB is on the line y={y_AB}, with A={A}, B={B}")

# 3. Define point M on segment AB
# Let's place M at a point dividing AB in ratio 1:3
M = 0.75 * A + 0.25 * B
print(f"Point M is on AB at {M}")

MA = np.linalg.norm(M - A)
MB = np.linalg.norm(M - B)
print(f"Lengths MA = {MA:.4f}, MB = {MB:.4f}")

# 4. Define two chords CD and EF of circle O passing through M
# Chord CD (e.g., slope 1)
m_CD = 1
c_CD = M[1] - m_CD * M[0]
intersections_CD = find_line_circle_intersections(m_CD, c_CD, O_h, O_k, O_r)
C, D = intersections_CD[0], intersections_CD[1]
print(f"Chord CD passes through M, with C={C}, D={D}")

# Chord EF (e.g., slope -0.5)
m_EF = -0.5
c_EF = M[1] - m_EF * M[0]
intersections_EF = find_line_circle_intersections(m_EF, c_EF, O_h, O_k, O_r)
E, F = intersections_EF[0], intersections_EF[1]
print(f"Chord EF passes through M, with E={E}, F={F}")

# 5. Find circumcircle of triangle EMD
h_EMD, k_EMD, r_EMD = find_circle_from_three_points(E, M, D)
if h_EMD is None:
    print("Points E, M, D are collinear, cannot continue.")
else:
    # 6. Find G
    G = find_intersections_with_horizontal_line(y_AB, h_EMD, k_EMD, r_EMD, M[0])
    print(f"Circumcircle of EMD intersects AB at G={G}")
    MG = np.linalg.norm(M - G)
    print(f"Length MG = {MG:.4f}")

    # 7. Find circumcircle of triangle CMF
    h_CMF, k_CMF, r_CMF = find_circle_from_three_points(C, M, F)
    if h_CMF is None:
        print("Points C, M, F are collinear, cannot continue.")
    else:
        # 8. Find H
        H = find_intersections_with_horizontal_line(y_AB, h_CMF, k_CMF, r_CMF, M[0])
        print(f"Circumcircle of CMF intersects AB at H={H}")
        MH = np.linalg.norm(M - H)
        print(f"Length MH = {MH:.4f}")
        
        # 9. Compare the results
        result1 = MG - MH
        result2 = MB - MA

        print("\n--- Final Comparison ---")
        print(f"The equation we want to express is: MG - MH")
        print(f"Let's test the hypothesis that MG - MH = MB - MA.")
        print(f"Calculated value of (MG - MH) is: {result1:.6f}")
        print(f"Calculated value of (MB - MA) is: {result2:.6f}")
        print(f"The difference between these two results is: {abs(result1 - result2):.6f}")
        
        # Final expression
        print(f"\nThus, we can express MG - MH in terms of MA and MB as:")
        print(f"The final equation is: {MG:.2f} - {MH:.2f} = {MB:.2f} - {MA:.2f}")
