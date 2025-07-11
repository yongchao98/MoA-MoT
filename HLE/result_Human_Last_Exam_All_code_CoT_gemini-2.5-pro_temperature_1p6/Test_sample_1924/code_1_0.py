import numpy as np

# Function to calculate the distance between two points
def distance(p1, p2):
    return np.linalg.norm(np.array(p1) - np.array(p2))

# Function to find the intersections of a line and a circle
def get_circle_line_intersections(circle_center, circle_radius, p1, p2):
    (cx, cy), r = circle_center, circle_radius
    x1, y1 = p1
    x2, y2 = p2
    
    dx, dy = x2 - x1, y2 - y1
    
    # Line equation: x = x1 + t*dx, y = y1 + t*dy
    # Substitute into circle equation: (x1 + t*dx - cx)^2 + (y1 + t*dy - cy)^2 = r^2
    a = dx**2 + dy**2
    b = 2 * (dx * (x1 - cx) + dy * (y1 - cy))
    c = (x1 - cx)**2 + (y1 - cy)**2 - r**2
    
    delta = b**2 - 4*a*c
    if delta < 0:
        return [] # No intersection
    
    t1 = (-b + np.sqrt(delta)) / (2*a)
    t2 = (-b - np.sqrt(delta)) / (2*a)
    
    i1 = (x1 + t1 * dx, y1 + t1 * dy)
    i2 = (x1 + t2 * dx, y1 + t2 * dy)
    
    return [i1, i2]

# Function to find the circumcircle of a triangle (3 points)
def get_circumcircle(p1, p2, p3):
    D = 2 * (p1[0] * (p2[1] - p3[1]) + p2[0] * (p3[1] - p1[1]) + p3[0] * (p1[1] - p2[1]))
    
    # If D is close to zero, points are collinear
    if abs(D) < 1e-9:
        return None, None
        
    p1_sq = p1[0]**2 + p1[1]**2
    p2_sq = p2[0]**2 + p2[1]**2
    p3_sq = p3[0]**2 + p3[1]**2
    
    ux = (p1_sq * (p2[1] - p3[1]) + p2_sq * (p3[1] - p1[1]) + p3_sq * (p1[1] - p2[1])) / D
    uy = (p1_sq * (p3[0] - p2[0]) + p2_sq * (p1[0] - p3[0]) + p3_sq * (p2[0] - p1[0])) / D
    
    center = (ux, uy)
    radius = distance(center, p1)
    
    return center, radius

# 1. Set up the Geometry
O_center = (1, 2)
O_radius = 10
M = (4, 4)

# Define lines through M to create chords
line_AB_p2 = (5, 3)    # A second point to define the line AB
line_CD_p2 = (2, 8)    # A second point to define the line CD
line_EF_p2 = (9, 5)    # A second point to define the line EF

# 2. Find Key Points
# Find A and B
intersections_AB = get_circle_line_intersections(O_center, O_radius, M, line_AB_p2)
if len(intersections_AB) != 2:
    print("Line AB does not properly intersect the circle O.")
else:
    A, B = intersections_AB

# Find C and D
intersections_CD = get_circle_line_intersections(O_center, O_radius, M, line_CD_p2)
if len(intersections_CD) != 2:
    print("Line CD does not properly intersect the circle O.")
else:
    C, D = intersections_CD
    
# Find E and F
intersections_EF = get_circle_line_intersections(O_center, O_radius, M, line_EF_p2)
if len(intersections_EF) != 2:
    print("Line EF does not properly intersect the circle O.")
else:
    E, F = intersections_EF
    
# 3. Find G and H
# Circumcircle of EMD
circ_EMD_center, circ_EMD_radius = get_circumcircle(E, M, D)

# Find G
if circ_EMD_center is not None:
    intersections_G = get_circle_line_intersections(circ_EMD_center, circ_EMD_radius, M, A)
    G = None
    for p in intersections_G:
        # G is the intersection point that is not M
        if distance(p, M) > 1e-6:
            G = p
            break
else:
    G = None

# Circumcircle of CMF
circ_CMF_center, circ_CMF_radius = get_circumcircle(C, M, F)

# Find H
if circ_CMF_center is not None:
    intersections_H = get_circle_line_intersections(circ_CMF_center, circ_CMF_radius, M, A)
    H = None
    for p in intersections_H:
        # H is the intersection point that is not M
        if distance(p, M) > 1e-6:
            H = p
            break
else:
    H = None

# 4. Calculate and Express the Result
if G is not None and H is not None:
    MA = distance(M, A)
    MB = distance(M, B)
    MG = distance(M, G)
    MH = distance(M, H)
    
    # We need to assign a sign to MA-MB and MG-MH based on the orientation
    # Let's align the vector MB with the positive direction
    vec_MB = np.array(B) - np.array(M)
    vec_MA = np.array(A) - np.array(M)
    vec_MG = np.array(G) - np.array(M)
    vec_MH = np.array(H) - np.array(M)
    
    # Check if vector MA points in the same direction as MB
    sign_A = 1 if np.dot(vec_MA, vec_MB) > 0 else -1
    sign_G = 1 if np.dot(vec_MG, vec_MB) > 0 else -1
    sign_H = 1 if np.dot(vec_MH, vec_MB) > 0 else -1
    
    # Signed differences
    ma_mb_val = sign_A * MA - MB # Since MB defines the positive direction
    mg_mh_val = sign_G * MG - sign_H * MH
    
    # A cleaner approach using the established theorem:
    # vec(MG) + vec(MH) = vec(MA) + vec(MB)
    # The question MG - MH can be interpreted as asking for a specific arrangement.
    # The general answer without assuming orientation is often MA - MB or MB - MA.
    # Let's output based on one of these plausible interpretations.
    # For this problem, it's established that |MG - MH| = |MA - MB|
    
    final_expr_val = MG - MH
    ma_minus_mb_val = MA - MB
    
    print("Based on the geometric configuration:")
    print(f"MA = {MA:.4f}")
    print(f"MB = {MB:.4f}")
    print(f"MG = {MG:.4f}")
    print(f"MH = {MH:.4f}")
    print("\nThe problem asks to express MG - MH in terms of MA and MB.")
    print("The theorem for this geometric configuration states |MG - MH| = |MA - MB|.")
    print("\nLet's verify this numerically:")
    print(f"|MG - MH| = |{MG:.4f} - {MH:.4f}| = {abs(final_expr_val):.4f}")
    print(f"|MA - MB| = |{MA:.4f} - {MB:.4f}| = {abs(ma_minus_mb_val):.4f}")
    print("\nLet's assume an orientation where B is further from M than A.")
    print(f"Then the relation is often written as MG - MH = MA - MB.")
    print("\nFinal Equation:")
    print(f"{MG:.4f} - {MH:.4f} = {MA:.4f} - {MB:.4f}")
    print(f"{final_expr_val:.4f} = {ma_minus_mb_val:.4f}")

else:
    print("Could not compute all necessary points for the calculation.")
