import math

def solve_weld_problem():
    """
    Calculates the required weld leg size for an eccentrically loaded connection.
    """
    # 1. Given Data
    P = 90000  # Force in N (90 kN)
    tau_perm = 250  # Permissible shear stress in N/mm^2
    l1 = 50  # Length of weld W1 in mm
    l2 = 50  # Length of weld W2 in mm
    l3 = 100 # Length of weld W3 in mm
    d1 = 50  # Horizontal length of the member at the weld
    d2 = 150 # Distance from the member end to the force application

    print("--- Given Data ---")
    print(f"Force (P) = {P / 1000} kN")
    print(f"Permissible Shear Stress (τ_perm) = {tau_perm} N/mm²")
    print(f"Weld dimensions (l1, l2, l3) = ({l1}, {l2}, {l3}) mm")
    print(f"Eccentricity dimensions = {d1} mm + {d2} mm")
    print("-" * 20 + "\n")

    # 2. Centroid Calculation
    # Set origin at the bottom of weld W3.
    # Centroid of W1: (l1/2, l3) = (25, 100)
    # Centroid of W2: (l2/2, 0)  = (25, 0)
    # Centroid of W3: (0, l3/2)   = (0, 50)
    L = l1 + l2 + l3
    x_bar = (l1 * (l1 / 2) + l2 * (l2 / 2) + l3 * 0) / L
    y_bar = (l1 * l3 + l2 * 0 + l3 * (l3 / 2)) / L
    
    print("--- Step 1: Centroid of the Weld Group (G) ---")
    print(f"Total weld length (L) = {l1} + {l2} + {l3} = {L} mm")
    print(f"x_bar = ({l1}*{l1/2} + {l2}*{l2/2} + {l3}*0) / {L} = {x_bar:.2f} mm")
    print(f"y_bar = ({l1}*{l3} + {l2}*0 + {l3}*{l3/2}) / {L} = {y_bar:.2f} mm")
    print(f"Centroid G is at ({x_bar:.2f}, {y_bar:.2f}) mm from the bottom-left corner.")
    print("-" * 20 + "\n")
    
    # 3. Eccentricity and Moment
    e = (d1 + d2) - x_bar
    M = P * e

    print("--- Step 2: Eccentricity and Moment (M) ---")
    print(f"Eccentricity (e) = ({d1} + {d2}) - x_bar = {e:.2f} mm")
    print(f"Moment (M) = P * e = {P} N * {e:.2f} mm = {M:.2f} N-mm")
    print("-" * 20 + "\n")

    # 4. Polar Moment of Inertia (per unit throat thickness)
    # Using Parallel Axis Theorem: I = I_c + A*d^2 (where A is length l for unit thickness)
    # For I_ux
    I_ux1 = (l1 * (l3 - y_bar)**2)
    I_ux2 = (l2 * (0 - y_bar)**2)
    I_ux3 = l3**3 / 12 + l3 * (l3/2 - y_bar)**2
    I_ux = I_ux1 + I_ux2 + I_ux3
    
    # For I_uy
    I_uy1 = l1**3 / 12 + l1 * (l1/2 - x_bar)**2
    I_uy2 = l2**3 / 12 + l2 * (l2/2 - x_bar)**2
    I_uy3 = l3 * (0 - x_bar)**2
    I_uy = I_uy1 + I_uy2 + I_uy3

    J_u = I_ux + I_uy

    print("--- Step 3: Unit Polar Moment of Inertia (Ju) ---")
    print(f"I_ux = {I_ux:.2f} mm³")
    print(f"I_uy = {I_uy:.2f} mm³")
    print(f"Ju = I_ux + I_uy = {I_ux:.2f} + {I_uy:.2f} = {J_u:.2f} mm³")
    print("-" * 20 + "\n")

    # 5. Stresses at Critical Point
    # The critical point is the top-right corner, furthest from the centroid.
    # Coordinates of critical point A: (l1, l3) = (50, 100)
    # Coordinates relative to centroid G:
    r_x = l1 - x_bar
    r_y = l3 - y_bar
    r = math.sqrt(r_x**2 + r_y**2)
    
    # Stresses are calculated per unit of throat thickness 't'
    tau_prime_y_per_t = -P / L  # Direct shear (downwards)
    tau_second_x_per_t = (M * r_y) / J_u # Torsional shear (horizontal component)
    tau_second_y_per_t = (-M * r_x) / J_u # Torsional shear (vertical component, clockwise moment)

    tau_total_x_per_t = tau_second_x_per_t
    tau_total_y_per_t = tau_prime_y_per_t + tau_second_y_per_t
    
    tau_res_per_t = math.sqrt(tau_total_x_per_t**2 + tau_total_y_per_t**2)

    print("--- Step 4: Stresses at Critical Point ---")
    print(f"Critical point A is at ({l1}, {l3}), which is ({r_x:.2f}, {r_y:.2f}) mm relative to the centroid.")
    print("Calculating stresses per unit throat thickness (N/mm per mm):")
    print(f"Primary shear stress component τ'_y = -P/L = {-P}/{L} = {tau_prime_y_per_t:.2f}")
    print(f"Secondary shear stress component τ''_x = (M*ry)/Ju = ({M:.0f}*{r_y:.2f})/{J_u:.0f} = {tau_second_x_per_t:.2f}")
    print(f"Secondary shear stress component τ''_y = (-M*rx)/Ju = (-{M:.0f}*{r_x:.2f})/{J_u:.0f} = {tau_second_y_per_t:.2f}")
    print(f"Total stress component τ_x_total = {tau_total_x_per_t:.2f}")
    print(f"Total stress component τ_y_total = {tau_prime_y_per_t:.2f} + {tau_second_y_per_t:.2f} = {tau_total_y_per_t:.2f}")
    print(f"Resultant stress (τ_res * t) = sqrt({tau_total_x_per_t:.2f}² + {tau_total_y_per_t:.2f}²) = {tau_res_per_t:.2f} N/mm")
    print("-" * 20 + "\n")
    
    # 6. Required Weld Size
    t = tau_res_per_t / tau_perm
    s = t / 0.707

    print("--- Step 5: Final Calculation for Weld Size ---")
    print("The final equation for the required leg size 's' is:")
    print("s = [ sqrt( ( (M*r_y)/J_u )^2 + ( -P/L - (M*r_x)/J_u )^2 ) ] / (τ_perm * 0.707)")
    print("\nSubstituting the calculated values:")
    # Using format to make it clear
    eq_str = (
        f"s = [ sqrt( ( ({M:.0f}*{r_y:.2f})/{J_u:.0f} )^2 + "
        f"( {-P/L:.0f} - ({M:.0f}*{r_x:.2f})/{J_u:.0f} )^2 ) ] / ( {tau_perm} * 0.707 )\n"
        f"s = [ sqrt( ( {tau_second_x_per_t:.1f} )^2 + ( {tau_prime_y_per_t:.0f} + {tau_second_y_per_t:.1f} )^2 ) ] / ( {tau_perm * 0.707:.1f} )\n"
        f"s = [ sqrt( {tau_second_x_per_t:.1f}**2 + {tau_total_y_per_t:.1f}**2 ) ] / ( {tau_perm * 0.707:.1f} )\n"
        f"s = {tau_res_per_t:.1f} / {tau_perm * 0.707:.1f}"
    )
    print(eq_str)
    
    print(f"\nRequired throat thickness (t) = {tau_res_per_t:.2f} / {tau_perm} = {t:.2f} mm")
    print(f"Required leg size (s) = t / 0.707 = {t:.2f} / 0.707 = {s:.2f} mm")
    print("\n" + "="*30)
    print(f"The final required size of the weld (leg) is {s:.2f} mm.")
    print("="*30)
    return s

# Execute the function and print the final answer in the required format
final_answer = solve_weld_problem()
print(f"\n<<<{final_answer:.2f}>>>")
