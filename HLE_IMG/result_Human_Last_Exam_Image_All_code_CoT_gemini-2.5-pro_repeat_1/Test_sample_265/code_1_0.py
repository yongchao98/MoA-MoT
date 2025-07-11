import math

def solve_weld_stress():
    """
    Calculates the required leg size for a welded connection under eccentric load.
    """
    # 1. Input parameters
    P = 90 * 1000      # Force in N
    tau_perm = 250     # Permissible shear stress in N/mm^2
    l_horz = 50        # Length of horizontal welds (W1, W2) in mm
    l_vert = 100       # Length of vertical weld (W3) in mm
    e_arm = 150        # Eccentricity arm from the connection plate in mm

    print("--- Input Parameters ---")
    print(f"Force (P): {P / 1000} kN")
    print(f"Permissible Shear Stress (τ_perm): {tau_perm} N/mm^2")
    print(f"Weld dimensions (l_horz x l_vert): {l_horz} mm x {l_vert} mm\n")

    # 2. Weld Geometry
    # (length, centroid_x, centroid_y) for each weld segment
    # Origin (0,0) is at the bottom-left corner of the weld.
    w1 = {'l': l_horz, 'cx': l_horz / 2, 'cy': l_vert}  # Top horizontal
    w2 = {'l': l_horz, 'cx': l_horz / 2, 'cy': 0}      # Bottom horizontal
    w3 = {'l': l_vert, 'cx': 0, 'cy': l_vert / 2}      # Vertical
    welds = [w1, w2, w3]

    L = sum(w['l'] for w in welds)
    sum_lx = sum(w['l'] * w['cx'] for w in welds)
    sum_ly = sum(w['l'] * w['cy'] for w in welds)
    
    # 3. Centroid of Weld Group (G)
    x_bar = sum_lx / L
    y_bar = sum_ly / L
    print("--- Weld Group Properties ---")
    print(f"Total weld length (L): {L:.2f} mm")
    print(f"Centroid (x_bar, y_bar): ({x_bar:.2f}, {y_bar:.2f}) mm\n")

    # 4. Eccentricity and Moment
    # Force is applied at x = l_horz + e_arm
    e = (l_horz + e_arm) - x_bar
    M = P * e
    print("--- Load Analysis ---")
    print(f"Eccentricity (e): {e:.2f} mm")
    print(f"Moment (M = P * e): {M:.2f} N-mm\n")

    # 5. Unit Polar Moment of Inertia (J_u)
    I_ux = 0
    I_uy = 0
    # For W1 (top horizontal)
    I_uy += (w1['l']**3 / 12) + w1['l'] * (w1['cx'] - x_bar)**2
    I_ux += w1['l'] * (w1['cy'] - y_bar)**2
    # For W2 (bottom horizontal)
    I_uy += (w2['l']**3 / 12) + w2['l'] * (w2['cx'] - x_bar)**2
    I_ux += w2['l'] * (w2['cy'] - y_bar)**2
    # For W3 (vertical)
    I_ux += (w3['l']**3 / 12) + w3['l'] * (w3['cy'] - y_bar)**2
    I_uy += w3['l'] * (w3['cx'] - x_bar)**2

    J_u = I_ux + I_uy
    print("--- Weld Inertia Properties (per unit thickness) ---")
    print(f"Unit Moment of Inertia about X-axis (I_ux): {I_ux:.2f} mm^3")
    print(f"Unit Moment of Inertia about Y-axis (I_uy): {I_uy:.2f} mm^3")
    print(f"Unit Polar Moment of Inertia (J_u): {J_u:.2f} mm^3\n")

    # 6. Stress Analysis at Critical Points
    # Critical points are the corners: A(50, 100), B(0, 100), C(0, 0), D(50, 0)
    critical_points = {'A': (l_horz, l_vert), 'D': (l_horz, 0)} # Max stress will be at A or D
    
    # Primary shear stress per unit throat thickness (acts downward)
    tau_p_per_t = P / L
    
    # We will find the resultant stress per unit throat thickness (tau_res_per_t)
    r_x = critical_points['D'][0] - x_bar
    r_y = critical_points['D'][1] - y_bar # Using point D, which has r_y as negative
    
    # Secondary shear stress components per unit throat thickness
    tau_sx_per_t = -M * r_y / J_u # x-component of secondary stress
    tau_sy_per_t = M * r_x / J_u  # y-component of secondary stress
    
    # Total shear stress components per unit throat thickness
    tau_x_total_per_t = tau_sx_per_t
    # Primary stress is downward (-y), secondary stress at D has an upward component (+y)
    # The worst case is when primary and secondary stresses add up.
    # Let's check point C/B where both y-components are downwards.
    # At C (0,0): r_x = -12.5, r_y = -50
    # tau_sy_per_t_C = M * (-12.5) / J_u (downward)
    # tau_y_total_per_t_C = -tau_p_per_t + tau_sy_per_t_C. This will be the maximum y-component.
    r_x_C = 0 - x_bar
    r_y_C = 0 - y_bar
    tau_sx_per_t_C = -M * r_y_C / J_u
    tau_sy_per_t_C = M * r_x_C / J_u
    tau_x_total_per_t = tau_sx_per_t_C
    tau_y_total_per_t = -tau_p_per_t + tau_sy_per_t_C

    # Resultant stress magnitude squared, per unit throat thickness squared
    tau_res_sq_per_t_sq = tau_x_total_per_t**2 + tau_y_total_per_t**2
    tau_res_per_t = math.sqrt(tau_res_sq_per_t_sq)

    print("--- Maximum Stress Calculation ---")
    print(f"The critical point for max stress is C at (0, 0) due to additive stress components.")
    print(f"Primary shear stress (per t): τ'_p = P / L = {P} / {L} = {tau_p_per_t:.2f} N/mm^2")
    print(f"Secondary shear stress components (per t) at C:")
    print(f"  τ'_sx = -M*r_y/J_u = {-M:.0f}*({r_y_C:.2f})/{J_u:.0f} = {tau_sx_per_t_C:.2f} N/mm^2")
    print(f"  τ'_sy =  M*r_x/J_u = {M:.0f}*({r_x_C:.2f})/{J_u:.0f} = {tau_sy_per_t_C:.2f} N/mm^2")
    print(f"Total stress components (per t):")
    print(f"  τ'_x = {tau_x_total_per_t:.2f} N/mm^2")
    print(f"  τ'_y = {-tau_p_per_t:.2f} + ({tau_sy_per_t_C:.2f}) = {tau_y_total_per_t:.2f} N/mm^2")
    print(f"Resultant max stress (per t): τ_max' = sqrt({tau_x_total_per_t:.2f}^2 + {tau_y_total_per_t:.2f}^2) = {tau_res_per_t:.2f} N/mm^2\n")

    # 7. Determine Weld Size
    # τ_perm = τ_max = τ_res_per_t / t
    t_req = tau_res_per_t / tau_perm
    s_req = t_req / math.cos(math.radians(45))
    
    print("--- Final Weld Size Calculation ---")
    print("Required throat thickness (t):")
    print(f"t = τ_max' / τ_perm = {tau_res_per_t:.2f} / {tau_perm} = {t_req:.4f} mm")
    print("\nRequired leg size (s):")
    print(f"s = t / cos(45°) = {t_req:.4f} / {math.cos(math.radians(45)):.4f} = {s_req:.4f} mm")
    print(f"\nThe required weld leg size is approximately {s_req:.1f} mm.")
    
    return s_req

# Run the solver
final_answer = solve_weld_stress()
<<<14.2>>>