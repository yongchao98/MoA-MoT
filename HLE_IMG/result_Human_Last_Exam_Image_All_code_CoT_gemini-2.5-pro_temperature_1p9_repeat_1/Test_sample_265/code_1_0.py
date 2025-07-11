import math

def solve_weld_size():
    """
    Calculates the required weld size for an eccentric load on a C-shaped weld group.
    """
    # 1. Define System Geometry and Loads
    P = 90000  # Force in N
    tau_allowable = 250  # Permissible shear stress in N/mm^2
    l_horiz = 50  # Length of horizontal welds (W1, W2) in mm
    l_vert = 100  # Length of vertical weld (W3) in mm
    # The eccentricity is the distance from the vertical support plate (where W3 is) to the force.
    e_arm_dist = 150 # mm
    
    print("--- Step 1: Given Parameters ---")
    print(f"Force (P): {P} N")
    print(f"Permissible Shear Stress (τ_allowable): {tau_allowable} N/mm^2")
    print(f"Weld dimensions (l_horiz, l_vert): {l_horiz} mm, {l_vert} mm")
    print("-" * 30 + "\n")

    # 2. Calculate Weld Group Properties (assuming unit throat thickness t=1mm)
    print("--- Step 2: Weld Group Properties (for unit throat thickness) ---")
    
    # Areas of each weld segment
    A1 = l_horiz * 1
    A2 = l_horiz * 1
    A3 = l_vert * 1
    A_total = A1 + A2 + A3

    # Centroid calculation (origin at the center of the vertical weld W3)
    # W1: x=25, y=-50; W2: x=25, y=50; W3: x=0, y=0
    x_bar = (A1 * (l_horiz / 2) + A2 * (l_horiz / 2) + A3 * 0) / A_total
    y_bar = (A1 * (-l_vert / 2) + A2 * (l_vert / 2) + A3 * 0) / A_total # Should be 0 by symmetry

    print(f"Centroid location (x_bar, y_bar): ({x_bar:.2f} mm, {y_bar:.2f} mm)")

    # Moment of Inertia about the centroid (using Parallel Axis Theorem: I = I_c + A*d^2)
    # I_x_unit
    Ix1 = (A1 * (-l_vert / 2 - y_bar)**2)
    Ix2 = (A2 * (l_vert / 2 - y_bar)**2)
    Ix3 = (1 * l_vert**3 / 12) + A3 * (0 - y_bar)**2
    Ix_total_unit = Ix1 + Ix2 + Ix3

    # I_y_unit
    Iy1 = (1 * l_horiz**3 / 12) + A1 * (l_horiz / 2 - x_bar)**2
    Iy2 = (1 * l_horiz**3 / 12) + A2 * (l_horiz / 2 - x_bar)**2
    Iy3 = A3 * (0 - x_bar)**2
    Iy_total_unit = Iy1 + Iy2 + Iy3

    # Polar Moment of Inertia
    J_unit = Ix_total_unit + Iy_total_unit
    
    print(f"Total unit moment of inertia about x-axis (Ix): {Ix_total_unit:.2f} mm^4")
    print(f"Total unit moment of inertia about y-axis (Iy): {Iy_total_unit:.2f} mm^4")
    print(f"Total unit polar moment of inertia (J): {J_unit:.2f} mm^4\n")

    # 3. Determine Forces on the Weld
    print("--- Step 3: Resultant Forces at Centroid ---")
    # Eccentricity
    e = e_arm_dist - x_bar
    # Torque (Moment)
    T = P * e
    
    print(f"Eccentricity (e): {e_arm_dist:.2f} - {x_bar:.2f} = {e:.2f} mm")
    print(f"Torque (T = P * e): {P} * {e:.2f} = {T:.2f} N-mm\n")

    # 4. Calculate Stresses (in terms of throat thickness 't')
    print("--- Step 4: Stress Calculation at Critical Point ---")
    
    # Critical point is the corner furthest from the centroid, which is at the ends of W1 and W2.
    # Let's analyze the top-right corner of weld W2.
    r_x = l_horiz - x_bar
    r_y = l_vert / 2 - y_bar
    
    print(f"Critical point coordinates relative to centroid (r_x, r_y): ({r_x:.2f} mm, {r_y:.2f} mm)")
    
    # Primary shear stress (acts downward, in -y direction)
    tau_prime_y = -P / A_total
    
    # Secondary shear stress components (from Torque T)
    # The stress vector is perpendicular to the radius vector, for a CW torque
    # τ''_x = T * r_y / J_unit
    # τ''_y = -T * r_x / J_unit
    tau_dprime_x = (T * r_y) / J_unit
    tau_dprime_y = (-T * r_x) / J_unit
    
    print("Stress values (multiplied by throat thickness 't'):")
    print(f"Primary shear stress component (τ'_y * t): {tau_prime_y:.2f} N/mm")
    print(f"Secondary shear stress component (τ''_x * t): {tau_dprime_x:.2f} N/mm")
    print(f"Secondary shear stress component (τ''_y * t): {tau_dprime_y:.2f} N/mm\n")

    # 5. Find the Maximum Resultant Stress
    print("--- Step 5: Resultant Shear Stress ---")
    # Total stress components
    tau_total_x = tau_dprime_x
    tau_total_y = tau_prime_y + tau_dprime_y

    # Resultant stress (multiplied by 't')
    tau_res_t = math.sqrt(tau_total_x**2 + tau_total_y**2)

    print(f"Resultant Stress (τ_res) = sqrt( (τ'_x + τ''_x)^2 + (τ'_y + τ''_y)^2 )")
    print(f"τ_res * t = sqrt( (0 + {tau_dprime_x:.2f})^2 + ({tau_prime_y:.2f} + {tau_dprime_y:.2f})^2 )")
    print(f"τ_res * t = sqrt( ({tau_total_x:.2f})^2 + ({tau_total_y:.2f})^2 )")
    print(f"τ_res * t = {tau_res_t:.2f} N/mm\n")

    # 6. Determine the Required Weld Size
    print("--- Step 6: Determine Weld Size 's' ---")
    
    # τ_res = τ_allowable
    # (tau_res_t / t) = tau_allowable
    t = tau_res_t / tau_allowable
    
    # Weld leg size 's' relates to throat thickness 't' by: t = s * cos(45°) = s / sqrt(2)
    s = t / math.cos(math.radians(45))
    
    print("Equating resultant stress to permissible stress to find throat thickness 't':")
    print(f"τ_res = {tau_res_t:.2f} / t = {tau_allowable} N/mm^2")
    print(f"t = {tau_res_t:.2f} / {tau_allowable:.2f} = {t:.3f} mm")
    print("\nConverting throat thickness 't' to weld leg size 's':")
    print(f"s = t / cos(45°) = {t:.3f} / {math.cos(math.radians(45)):.3f} = {s:.3f} mm")
    
    return s

# --- Execute the calculation and print the final answer ---
final_weld_size = solve_weld_size()
print(f"\nThe determined size of the weld (leg) is approximately {final_weld_size:.2f} mm.")
print(f"<<<{final_weld_size:.2f}>>>")
