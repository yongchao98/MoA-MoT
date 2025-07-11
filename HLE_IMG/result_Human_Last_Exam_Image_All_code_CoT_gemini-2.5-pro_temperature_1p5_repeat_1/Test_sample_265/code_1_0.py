import math

def solve_weld_connection():
    """
    Calculates the required leg size for a welded connection under eccentric loading.
    """
    # Step 1: Define the given parameters
    P = 90 * 1000  # Force in N
    tau_perm = 250  # Permissible shear stress in N/mm^2
    l_hor = 50  # Length of horizontal welds (W1, W2) in mm
    l_ver = 100  # Length of vertical weld (W3) in mm
    dist1 = 50  # Horizontal distance from vertical plate to end of horizontal weld
    dist2 = 150 # Horizontal distance from end of horizontal weld to force application

    print("### Analysis of Welded Connection ###\n")

    print("--- Step 1: Given Parameters ---")
    print(f"Force, P = {P/1000} kN = {P} N")
    print(f"Permissible shear stress, τ_perm = {tau_perm} N/mm^2")
    print(f"Length of horizontal welds, L1=L2 = {l_hor} mm")
    print(f"Length of vertical weld, L3 = {l_ver} mm")
    print("-" * 40)

    # Step 2: Calculate properties of the weld group
    L_w = 2 * l_hor + l_ver
    print("\n--- Step 2: Properties of the Weld Group ---")
    print(f"Total weld length, L_w = 2 * {l_hor} + {l_ver} = {L_w} mm")

    # Centroid of the weld group (origin at the bottom-left corner of W3)
    x_bar = (l_hor * (l_hor / 2) + l_hor * (l_hor / 2) + l_ver * 0) / L_w
    y_bar = l_ver / 2
    print(f"Centroid of weld group (Gx, Gy) = ({x_bar:.2f} mm, {y_bar:.2f} mm)")
    print("-" * 40)

    # Step 3: Calculate Moment and Eccentricity
    e = (dist1 + dist2) - x_bar
    M = P * e
    print("\n--- Step 3: Moment Calculation ---")
    print(f"Eccentricity, e = ({dist1} + {dist2}) - {x_bar:.2f} = {e:.2f} mm")
    print(f"Moment, M = P * e = {P} N * {e:.2f} mm = {M:.2f} N-mm (Clockwise)")
    print("-" * 40)

    # Step 4: Calculate Polar Moment of Inertia (per unit throat thickness)
    # Moment of inertia about x-axis (I_xx_u)
    I_xx_u_1_2 = 2 * (l_hor * (y_bar)**2)
    I_xx_u_3 = (l_ver**3) / 12
    I_xx_u = I_xx_u_1_2 + I_xx_u_3

    # Moment of inertia about y-axis (I_yy_u)
    I_yy_u_1_2 = 2 * ((l_hor**3 / 12) + l_hor * (l_hor/2 - x_bar)**2)
    I_yy_u_3 = l_ver * x_bar**2
    I_yy_u = I_yy_u_1_2 + I_yy_u_3

    # Polar moment of inertia per unit thickness
    J_u = I_xx_u + I_yy_u
    print("\n--- Step 4: Polar Moment of Inertia (per unit throat thickness 't') ---")
    print(f"I_xx_u = {I_xx_u_1_2:.2f} + {I_xx_u_3:.2f} = {I_xx_u:.2f} mm^3")
    print(f"I_yy_u = {I_yy_u_1_2:.2f} + {I_yy_u_3:.2f} = {I_yy_u:.2f} mm^3")
    print(f"Unit polar moment of inertia, J_u = {I_xx_u:.2f} + {I_yy_u:.2f} = {J_u:.2f} mm^3")
    print("-" * 40)


    # Step 5: Calculate Stresses at the Critical Point
    print("\n--- Step 5: Stress Calculation ---")
    # Critical point is where stresses combine to a maximum. For this geometry, it's the
    # top-right or bottom-right corner. Let's analyze the bottom-right corner.
    x_crit_rel = l_hor - x_bar
    y_crit_rel = 0 - y_bar
    print(f"Analyzing critical point (bottom-right corner) relative to centroid:")
    print(f"(x_rel, y_rel) = ({x_crit_rel:.2f} mm, {y_crit_rel:.2f} mm)")

    # Primary shear stress per unit throat thickness 't' (acts downwards)
    tau_prime_y_per_t = P / L_w
    print(f"\nPrimary shear stress (vertical), τ'_y = -(P / L_w) / t = -({P} / {L_w}) / t = -{tau_prime_y_per_t:.2f} / t")

    # Secondary (torsional) shear stress components per unit throat thickness 't'
    # For a clockwise moment M: τ''_x = M*y/J, τ''_y = -M*x/J
    tau_second_x_per_t = M * y_crit_rel / J_u
    tau_second_y_per_t = -(M * x_crit_rel / J_u)
    print(f"Secondary shear stress (horizontal), τ''_x = (M*y_rel/J_u)/t = {tau_second_x_per_t:.2f} / t")
    print(f"Secondary shear stress (vertical), τ''_y = -(M*x_rel/J_u)/t = {tau_second_y_per_t:.2f} / t")

    # Total stress components per unit throat thickness 't'
    tau_total_x_per_t = tau_second_x_per_t
    tau_total_y_per_t = -tau_prime_y_per_t + tau_second_y_per_t

    print("\nResultant Stress Components at Critical Point:")
    print(f"τ_x_total = τ''_x = {tau_total_x_per_t:.2f} / t N/mm^2")
    print(f"τ_y_total = τ'_y + τ''_y = ( -{tau_prime_y_per_t:.2f} {tau_second_y_per_t:+.2f} ) / t = {tau_total_y_per_t:.2f} / t N/mm^2")

    # Resultant stress per unit throat thickness 't'
    C = math.sqrt(tau_total_x_per_t**2 + tau_total_y_per_t**2)
    print(f"\nResultant stress, τ_res = sqrt( (τ_x_total)² + (τ_y_total)² ) = sqrt( ({tau_total_x_per_t:.2f}/t)² + ({tau_total_y_per_t:.2f}/t)² )")
    print(f"τ_res = (1/t) * sqrt({tau_total_x_per_t**2:.2f} + {tau_total_y_per_t**2:.2f}) = {C:.2f} / t N/mm^2")
    print("-" * 40)

    # Step 6: Determine the Required Weld Size
    print("\n--- Step 6: Weld Size Calculation ---")
    t_req = C / tau_perm
    print(f"Set τ_res = τ_perm => {C:.2f} / t = {tau_perm}")
    print(f"Required throat thickness, t = {C:.2f} / {tau_perm} = {t_req:.4f} mm")

    s_req = t_req / 0.707
    print(f"\nFinal Equation: Weld leg size, s = t / 0.707 = {t_req:.4f} / 0.707")
    print(f"The required weld size (leg) is s = {s_req:.2f} mm.")

solve_weld_connection()