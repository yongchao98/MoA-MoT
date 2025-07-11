import math

def solve_weld_size():
    """
    Calculates the required weld leg size for an eccentrically loaded connection.
    """
    # --- Input Data from the problem statement ---
    P = 90 * 1000       # Force in N
    tau_permissible = 250 # Permissible shear stress in N/mm^2
    l_horiz = 50        # Length of horizontal welds (W1, W2) in mm
    l_vert = 100        # Length of vertical weld (W3) in mm
    e_arm = 150         # Eccentricity arm from the end of the weld in mm

    # --- Step 1: Find Centroid of Weld Group ---
    # We assume a unit throat thickness (t=1) for geometric calculations.
    # The origin (0,0) is set at the bottom-left corner of the vertical weld W3.
    # Segments are defined as [length, centroid_x, centroid_y]
    w1 = [l_horiz, l_horiz / 2, l_vert]   # Top weld W1
    w2 = [l_horiz, l_horiz / 2, 0]        # Bottom weld W2
    w3 = [l_vert, 0, l_vert / 2]          # Vertical weld W3
    welds = [w1, w2, w3]

    # Calculate total length and locate the centroid (x_bar, y_bar)
    L_total = 0
    sum_lx = 0
    sum_ly = 0
    for length, cx, cy in welds:
        L_total += length
        sum_lx += length * cx
        sum_ly += length * cy

    x_bar = sum_lx / L_total
    y_bar = sum_ly / L_total

    print("--- Step 1: Find Centroid of Weld Group ---")
    print(f"The weld group consists of two horizontal welds of {l_horiz} mm and one vertical weld of {l_vert} mm.")
    print(f"Total weld length L = {w1[0]} + {w2[0]} + {w3[0]} = {L_total:.2f} mm")
    print(f"The centroid G of the weld group is at ({x_bar:.2f} mm, {y_bar:.2f} mm).")

    # --- Step 2: Calculate Eccentricity and Moment ---
    # The force is applied at a horizontal position of l_horiz + e_arm
    force_pos_x = l_horiz + e_arm
    e = force_pos_x - x_bar
    M = P * e

    print("\n--- Step 2: Calculate Eccentricity and Moment ---")
    print(f"The eccentric force is P = {P/1000} kN = {P} N.")
    print(f"The eccentricity e = (distance from W3 to force) - x_bar = ({l_horiz} + {e_arm}) - {x_bar:.2f} = {e:.2f} mm.")
    print(f"This creates a moment M = P * e = {P} * {e:.2f} = {M:.2f} N-mm.")

    # --- Step 3: Calculate Polar Moment of Inertia (for unit thickness) ---
    I_xx_total = 0
    I_yy_total = 0
    # Using the parallel axis theorem: I = I_c + A*d^2, where Area A = length * (unit thickness=1)
    # W1 (top horizontal)
    I_xx1 = (l_horiz * 1**3 / 12) + (l_horiz * 1) * (w1[2] - y_bar)**2 # Term l*t^3/12 is negligible
    I_yy1 = (1 * l_horiz**3 / 12) + (l_horiz * 1) * (w1[1] - x_bar)**2
    # W2 (bottom horizontal)
    I_xx2 = (l_horiz * 1**3 / 12) + (l_horiz * 1) * (w2[2] - y_bar)**2 # Negligible
    I_yy2 = (1 * l_horiz**3 / 12) + (l_horiz * 1) * (w2[1] - x_bar)**2
    # W3 (vertical)
    I_xx3 = (1 * l_vert**3 / 12) + (l_vert * 1) * (w3[2] - y_bar)**2
    I_yy3 = (l_vert * 1**3 / 12) + (l_vert * 1) * (w3[1] - x_bar)**2 # Negligible
    
    I_xx_total = I_xx1 + I_xx2 + I_xx3
    I_yy_total = I_yy1 + I_yy2 + I_yy3
    J_unit = I_xx_total + I_yy_total

    print("\n--- Step 3: Calculate Polar Moment of Inertia (J_unit) ---")
    print("Using the parallel axis theorem for each segment (assuming unit throat thickness):")
    print(f"Second moment of area about x-axis, I_xx = {I_xx_total:.2f} mm^4")
    print(f"Second moment of area about y-axis, I_yy = {I_yy_total:.2f} mm^4")
    print(f"Polar moment of inertia J_unit = I_xx + I_yy = {J_unit:.2f} mm^4")

    # --- Step 4: Calculate Stresses at Critical Point ---
    # The primary shear stress is uniform and acts downwards (-y direction).
    tau_prime_per_t = P / L_total

    # The critical points are the corners furthest from the centroid, A(50, 100) and B(50, 0).
    # Let's analyze point A(50, 100), as it's a common point of failure.
    px_crit, py_crit = (l_horiz, l_vert)
    rx_crit = px_crit - x_bar
    ry_crit = py_crit - y_bar

    # Secondary shear stress components (per unit thickness 't')
    tau_sec_x_per_t = (M * ry_crit) / J_unit
    tau_sec_y_per_t = -(M * rx_crit) / J_unit

    # Total stress components
    tau_total_x_per_t = tau_sec_x_per_t
    tau_total_y_per_t = -tau_prime_per_t + tau_sec_y_per_t

    # Resultant stress magnitude (per unit thickness)
    max_res_stress_per_t = math.sqrt(tau_total_x_per_t**2 + tau_total_y_per_t**2)

    print("\n--- Step 4: Calculate Stresses at Critical Point ---")
    print("The maximum stress occurs at the corners. We analyze the top-right corner A(50, 100).")
    print(f"Primary shear stress (per unit 't'): τ'_y = -P / L = -{P}/{L_total:.2f} = {-tau_prime_per_t:.2f} N/mm per mm")
    print("Secondary (torsional) shear stress components (per unit 't'):")
    print(f"  τ''_x = (M * r_y) / J_unit = ({M:.0f} * {ry_crit:.2f}) / {J_unit:.0f} = {tau_sec_x_per_t:.2f}")
    print(f"  τ''_y = -(M * r_x) / J_unit = -({M:.0f} * {rx_crit:.2f}) / {J_unit:.0f} = {tau_sec_y_per_t:.2f}")
    print("Total stress components (per unit 't'):")
    print(f"  τ_total_x = 0 + {tau_sec_x_per_t:.2f} = {tau_total_x_per_t:.2f}")
    print(f"  τ_total_y = {-tau_prime_per_t:.2f} + {tau_sec_y_per_t:.2f} = {tau_total_y_per_t:.2f}")
    print("Maximum resultant shear stress (per unit 't'):")
    print(f"  τ_res / t = sqrt({tau_total_x_per_t:.2f}² + {tau_total_y_per_t:.2f}²) = {max_res_stress_per_t:.2f} N/mm per mm")
    
    # --- Step 5: Determine Required Weld Size ---
    # τ_res = max_res_stress_per_t / t <= τ_permissible
    t_required = max_res_stress_per_t / tau_permissible
    s_required = t_required / math.cos(math.radians(45))

    print("\n--- Step 5: Determine Required Weld Size ---")
    print("The maximum resultant stress must not exceed the permissible stress.")
    print("Final Equation (Throat Thickness 't'):")
    print(f"t = (Resultant Stress per unit 't') / (Permissible Stress) = {max_res_stress_per_t:.2f} / {tau_permissible} = {t_required:.4f} mm")
    print("\nFinal Equation (Leg Size 's'):")
    print("For a fillet weld, t = s * cos(45°). Therefore, s = t / cos(45°).")
    print(f"s = {t_required:.4f} / {math.cos(math.radians(45)):.4f} = {s_required:.4f} mm")
    print("\n-----------------------------------------------------")
    print(f"The required size of the weld (leg) is approximately {s_required:.2f} mm.")
    print("-----------------------------------------------------")
    
    return s_required

if __name__ == "__main__":
    final_answer = solve_weld_size()
    # The final answer is also printed inside the function.
    # To conform to the format, let's also prepare the final numerical value.
    # print(f"<<<{final_answer:.2f}>>>") # This would be the final output line