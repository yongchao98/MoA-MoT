import math

def solve_weld_size():
    """
    Calculates the required weld leg size for an eccentrically loaded connection.
    """
    # --- Step 1: Define Given Parameters ---
    P = 90 * 1000  # Force in N
    b = 50.0       # Length of horizontal welds in mm
    d = 100.0      # Length of vertical weld in mm
    e_prime = 150.0  # Distance from weld end to force in mm
    tau_perm = 250.0 # Permissible shear stress in N/mm^2

    print("--- Step 1: System Parameters ---")
    print(f"Force, P = {P/1000} kN = {P:.0f} N")
    print(f"Length of horizontal welds (W1, W2), b = {b:.0f} mm")
    print(f"Length of vertical weld (W3), d = {d:.0f} mm")
    print(f"Permissible shear stress, τ_perm = {tau_perm:.0f} N/mm^2\n")

    # --- Step 2: Weld Group Properties ---
    L = 2 * b + d
    # Centroid calculation (origin at the bottom of the vertical weld)
    x_bar = (b * (b / 2) + b * (b / 2) + d * 0) / L
    y_bar = (b * d + b * 0 + d * (d / 2)) / L
    print("--- Step 2: Weld Group Centroid (G) ---")
    print(f"Total weld length, L = 2*b + d = 2*{b:.0f} + {d:.0f} = {L:.0f} mm")
    print(f"Centroid x-coordinate, x_bar = (2 * {b:.0f} * {b/2:.1f}) / {L:.0f} = {x_bar:.2f} mm")
    print(f"Centroid y-coordinate, y_bar = ({b:.0f}*{d:.0f} + {d:.0f}*{d/2:.0f}) / {L:.0f} = {y_bar:.2f} mm\n")

    # --- Step 3: Moment Calculation ---
    e = (b + e_prime) - x_bar
    M = P * e
    print("--- Step 3: Moment about Centroid ---")
    print(f"Eccentricity, e = ({b:.0f} + {e_prime:.0f}) - {x_bar:.2f} = {e:.2f} mm")
    print(f"Moment, M = P * e = {P:.0f} N * {e:.2f} mm = {M:.2f} N-mm\n")

    # --- Step 4: Polar Moment of Inertia ---
    # Moments of inertia per unit throat thickness (I_u), using parallel axis theorem
    Ixx1 = b * (d - y_bar)**2
    Ixx2 = b * (0 - y_bar)**2
    Ixx3 = d**3 / 12
    I_xx_u = Ixx1 + Ixx2 + Ixx3
    Iyy1 = (b**3 / 12) + b * (b/2 - x_bar)**2
    Iyy2 = (b**3 / 12) + b * (b/2 - x_bar)**2
    Iyy3 = d * (0 - x_bar)**2
    I_yy_u = Iyy1 + Iyy2 + Iyy3
    J_u = I_xx_u + I_yy_u
    print("--- Step 4: Polar Moment of Inertia (per unit throat) ---")
    print(f"I_xx_u = {I_xx_u:.2f} mm^3")
    print(f"I_yy_u = {I_yy_u:.2f} mm^3")
    print(f"Polar Moment of Inertia, J_u = I_xx_u + I_yy_u = {J_u:.2f} mm^3\n")

    # --- Step 5, 6: Stresses at Critical Point (Top-Right Corner) ---
    # Coordinates of critical point relative to the centroid
    x_c = b - x_bar
    y_c = d - y_bar
    
    # Primary shear stress (downwards) per unit throat
    tau_prime_y = P / L
    # Secondary shear stress components (right and down) per unit throat
    tau_double_prime_x = (M * y_c) / J_u
    tau_double_prime_y = (M * x_c) / J_u
    
    # Total stress components
    tau_total_x = tau_double_prime_x
    tau_total_y = tau_prime_y + tau_double_prime_y
    tau_max_per_t = math.sqrt(tau_total_x**2 + tau_total_y**2)

    print("--- Step 5 & 6: Stress Analysis at Critical Point ---")
    print(f"Primary shear stress (downward), τ'_y = P / L = {P:.0f} / {L:.0f} = {tau_prime_y:.2f} / t")
    print(f"Secondary shear stress (horizontal), τ''_x = (M*y_c)/J_u = ({M:.0f}*{y_c:.2f})/{J_u:.2f} = {tau_double_prime_x:.2f} / t")
    print(f"Secondary shear stress (vertical), τ''_y = (M*x_c)/J_u = ({M:.0f}*{x_c:.2f})/{J_u:.2f} = {tau_double_prime_y:.2f} / t")
    print(f"Resultant stress, τ_max = sqrt((τ''_x)^2 + (τ'_y + τ''_y)^2)")
    print(f"τ_max = sqrt(({tau_double_prime_x:.2f}/t)^2 + (({tau_prime_y:.2f} + {tau_double_prime_y:.2f})/t)^2)")
    print(f"τ_max = {tau_max_per_t:.2f} / t\n")

    # --- Step 7: Determine Weld Size ---
    # τ_max = τ_perm => tau_max_per_t / t = tau_perm
    t = tau_max_per_t / tau_perm
    s = t / 0.7071
    print("--- Step 7: Final Weld Size Calculation ---")
    print(f"Required throat thickness, t = τ_max_per_t / τ_perm = {tau_max_per_t:.2f} / {tau_perm:.0f} = {t:.4f} mm")
    print("Required leg size, s = t / 0.7071")
    print(f"s = {t:.4f} / 0.7071 = {s:.4f} mm")
    print("\n-------------------------------------------")
    print(f"The required size of the weld (leg) is {s:.2f} mm.")
    print("-------------------------------------------")
    return s

# Run the solver
final_answer = solve_weld_size()
<<<17.13>>>