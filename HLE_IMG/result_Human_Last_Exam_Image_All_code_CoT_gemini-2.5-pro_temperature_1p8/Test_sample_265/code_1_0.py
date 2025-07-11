import math

def solve_weld_problem():
    """
    Calculates the required weld size for a connection subjected to an eccentric force.
    """
    # --- 1. Define Givens and Constants ---
    P = 90 * 1000      # Force in N
    tau_per = 250      # Permissible shear stress in N/mm^2
    b = 50.0           # Length of horizontal welds W1, W2 in mm
    d = 100.0          # Length of vertical weld W3 in mm
    l_offset = 150.0   # Offset distance from weld group edge to force application point

    # Throat size to leg size conversion factor for a fillet weld
    K = 0.707

    print("### WELD STRESS ANALYSIS ###")
    print("\n--- Step 1: Given Parameters ---")
    print(f"Applied eccentric force, P = {P/1000} kN = {P} N")
    print(f"Permissible shear stress, τ_per = {tau_per} N/mm²")
    print(f"Horizontal weld length (W1, W2), b = {b} mm")
    print(f"Vertical weld length (W3), d = {d} mm")

    # --- 2. Calculate Centroid of Weld Group ---
    # We model the weld as a line group.
    # Origin (0,0) is at the bottom-left corner of the weld.
    L_total = b + b + d
    # x-coordinate of centroid G
    x_g = (b * (b / 2) + b * (b / 2) + d * 0) / L_total
    # y-coordinate of centroid G (by symmetry)
    y_g = d / 2

    print("\n--- Step 2: Centroid of the Weld Group ---")
    print(f"Total weld length, L_w = {b} + {b} + {d} = {L_total} mm")
    print(f"x-coordinate of centroid, x_g = ({b} * {b/2} + {b} * {b/2} + {d} * 0) / {L_total} = {x_g:.2f} mm")
    print(f"y-coordinate of centroid, y_g = {d} / 2 = {y_g:.2f} mm")
    print(f"Centroid G is at ({x_g:.2f}, {y_g:.2f})")

    # --- 3. Calculate Stresses (per unit throat thickness 't') ---
    # The weld connection resists the external force, creating an upward primary shear force
    # and a counter-clockwise (CCW) moment M.

    # Primary shear stress (τ') acts UPWARDS. It's P divided by total throat area A_w = L_total * t.
    tau_prime_per_t = P / L_total
    tau_prime_x_per_t = 0
    tau_prime_y_per_t = tau_prime_per_t # Positive, as it acts upwards to resist the downward force

    print("\n--- Step 3: Stress Calculations ---")
    print("\n3a. Primary Shear Stress (due to direct force P)")
    print(f"The primary shear stress is uniform and acts upwards to balance P.")
    print(f"τ' = P / A_w = {P} / ({L_total}*t) = {tau_prime_per_t:.2f} / t  (N/mm²)")

    # Secondary shear stress (τ'') due to Moment M
    eccentricity = (b + l_offset) - x_g
    M = P * eccentricity  # Counter-clockwise moment resisted by the weld group

    print("\n3b. Secondary Shear Stress (due to moment M)")
    print(f"Eccentricity, e = ({b} + {l_offset}) - {x_g:.2f} = {eccentricity:.2f} mm")
    print(f"Resisting Moment, M = P * e = {P} N * {eccentricity:.2f} mm = {M:.2f} N-mm (CCW)")

    # Unit polar moment of inertia, J_u (moment of inertia per unit throat thickness t)
    # I_xx_u for welds about centroidal x-axis (y = y_g)
    I_xx_u_1 = b * (d - y_g)**2
    I_xx_u_2 = b * (0 - y_g)**2
    I_xx_u_3 = d**3 / 12
    I_xx_u = I_xx_u_1 + I_xx_u_2 + I_xx_u_3
    # I_yy_u for welds about centroidal y-axis (x = x_g)
    I_yy_u_1 = (b**3 / 12) + b * (b/2 - x_g)**2
    I_yy_u_2 = (b**3 / 12) + b * (b/2 - x_g)**2
    I_yy_u_3 = d * (0 - x_g)**2
    I_yy_u = I_yy_u_1 + I_yy_u_2 + I_yy_u_3

    J_u = I_xx_u + I_yy_u

    print(f"Polar Moment of Inertia (per t), J_u = I_xx_u + I_yy_u")
    print(f"J_u = {I_xx_u:.2f} + {I_yy_u:.2f} = {J_u:.2f} mm³")

    # --- 4. Find Maximum Resultant Stress at Critical Points ---
    # Critical points are corners farthest from the centroid. Point A (b,d) and B (b,0) are chosen.
    xA, yA = b, d
    rAx, rAy = xA - x_g, yA - y_g
    rA = math.sqrt(rAx**2 + rAy**2)

    # Secondary shear magnitude (per t) at A
    tau_double_prime_A_per_t = (M * rA) / J_u
    
    # Secondary shear vector components (per t) at A
    # For a CCW moment, the stress vector direction at a point (rx,ry) from the centroid is (-ry, rx)
    tau_double_prime_Ax_per_t = tau_double_prime_A_per_t * (-rAy / rA)
    tau_double_prime_Ay_per_t = tau_double_prime_A_per_t * (rAx / rA)
    
    # Total resultant stress (per t) at A
    tau_total_Ax_per_t = tau_prime_x_per_t + tau_double_prime_Ax_per_t
    tau_total_Ay_per_t = tau_prime_y_per_t + tau_double_prime_Ay_per_t
    tau_max_per_t = math.sqrt(tau_total_Ax_per_t**2 + tau_total_Ay_per_t**2)

    print("\n--- Step 4: Maximum Resultant Stress ---")
    print(f"Analyzing critical point A({xA:.0f},{yA:.0f}), which is furthest from G({x_g:.2f},{y_g:.2f})")
    print(f"Radius, r = sqrt(({rAx:.2f})² + ({rAy:.2f})²) = {rA:.2f} mm")
    print(f"\nResultant stress components (per unit thickness t):")
    print(f"τ'_x={tau_prime_x_per_t:.2f}, τ'_y={tau_prime_y_per_t:.2f}")
    print(f"τ''_x={tau_double_prime_Ax_per_t:.2f}, τ''_y={tau_double_prime_Ay_per_t:.2f}")
    print(f"τ_total_x = {tau_prime_x_per_t:.2f} + {tau_double_prime_Ax_per_t:.2f} = {tau_total_Ax_per_t:.2f}")
    print(f"τ_total_y = {tau_prime_y_per_t:.2f} + {tau_double_prime_Ay_per_t:.2f} = {tau_total_Ay_per_t:.2f}")

    print(f"\nMax resultant shear stress |τ_max| = sqrt(({tau_total_Ax_per_t:.2f})² + ({tau_total_Ay_per_t:.2f})²) / t")
    print(f"|τ_max| = {tau_max_per_t:.2f} / t (N/mm²)")
    
    # --- 5. Determine the required Weld Size ---
    # τ_max = τ_per => τ_max_per_t / t = tau_per
    throat_t = tau_max_per_t / tau_per
    leg_s = throat_t / K

    print("\n--- Step 5: Final Weld Size Calculation ---")
    print(f"Set max stress equal to permissible stress: τ_max = τ_per")
    print(f"{tau_max_per_t:.2f} / t = {tau_per}")
    print(f"Required throat thickness, t = {tau_max_per_t:.2f} / {tau_per} = {throat_t:.2f} mm")
    print(f"\nRequired weld size (leg), s = t / K")
    print(f"s = {throat_t:.2f} / {K} = {leg_s:.2f} mm")
    
    # Final answer for automated checking
    return leg_s

if __name__ == '__main__':
    result = solve_weld_problem()
    print(f"\n<<< {result:.2f} >>>")
