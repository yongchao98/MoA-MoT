import math

def calculate_design_resistance():
    """
    Calculates the design bearing resistance for a square pad foundation
    under undrained conditions as per Eurocode 7.
    """

    # --- Given Parameters ---
    # Foundation dimensions
    B = 2.2  # Foundation width, m
    L = 2.2  # Foundation length, m
    D = 2.0  # Foundation depth, m

    # Soil properties
    c_uk = 57.0  # Characteristic undrained shear strength, kN/m^2
    depth_gwl = 0.6  # Groundwater level depth, m
    gamma_clay_above = 18.0  # Unit weight of clay above GWL, kN/m^3
    gamma_clay_below = 20.0  # Unit weight of clay below GWL, kN/m^3
    
    # Partial factor of safety for bearing resistance (EC7, DA1-C1, undrained)
    gamma_R_v = 1.4

    print("--- Step 1: Calculate foundation base area (A') ---")
    A_prime = B * L
    print(f"A' = B * L = {B} m * {L} m = {A_prime:.2f} m^2\n")

    print("--- Step 2: Calculate total overburden pressure (q) at foundation base ---")
    pressure_above_gwl = depth_gwl * gamma_clay_above
    depth_below_gwl = D - depth_gwl
    pressure_below_gwl = depth_below_gwl * gamma_clay_below
    q = pressure_above_gwl + pressure_below_gwl
    print(f"q = (depth_gwl * gamma_above) + ((D - depth_gwl) * gamma_below)")
    print(f"q = ({depth_gwl} * {gamma_clay_above}) + (({D} - {depth_gwl}) * {gamma_clay_below})")
    print(f"q = {pressure_above_gwl:.2f} + {pressure_below_gwl:.2f} = {q:.2f} kN/m^2\n")

    print("--- Step 3: Determine bearing capacity factors ---")
    s_c = 1.0 + 0.2 * (B / L)
    i_c = 1.0
    b_c = 1.0
    N_c_term = math.pi + 2
    print(f"Shape factor (s_c) for square footing = 1 + 0.2 * (B/L) = {s_c:.1f}")
    print(f"Load inclination factor (i_c) for vertical load = {i_c:.1f}")
    print(f"Base inclination factor (b_c) for horizontal base = {b_c:.1f}")
    print(f"Bearing capacity factor term (π + 2) = {N_c_term:.4f}\n")

    print("--- Step 4: Calculate characteristic bearing resistance (R_k) ---")
    # Contribution from cohesion
    cohesion_term = N_c_term * c_uk * s_c * i_c * b_c
    # Gross bearing pressure
    q_k_gross = cohesion_term + q
    # Gross bearing resistance
    R_k = q_k_gross * A_prime
    print("Formula: R_k = A' * [ (π + 2) * c_uk * s_c * i_c * b_c + q ]")
    print(f"R_k = {A_prime:.2f} * [ ({math.pi:.4f} + 2) * {c_uk} * {s_c:.1f} * {i_c:.1f} * {b_c:.1f} + {q:.2f} ]")
    print(f"R_k = {A_prime:.2f} * [ {cohesion_term:.2f} + {q:.2f} ]")
    print(f"R_k = {A_prime:.2f} * {q_k_gross:.2f}")
    print(f"R_k = {R_k:.2f} kN\n")

    print("--- Step 5: Calculate design bearing resistance (R_d) ---")
    R_d = R_k / gamma_R_v
    print("Formula: R_d = R_k / γ_R;v")
    print(f"The final equation for the design resistance is:")
    print(f"R_d = (A' * [ (π + 2) * c_uk * s_c + q ]) / γ_R;v")
    print(f"R_d = ({A_prime:.2f} m^2 * [ ({math.pi:.4f} + 2) * {c_uk} kN/m^2 * {s_c:.1f} + {q:.2f} kN/m^2 ]) / {gamma_R_v}")
    print(f"R_d = {R_k:.2f} kN / {gamma_R_v}")
    print(f"R_d = {R_d:.2f} kN")
    
    return R_d

if __name__ == '__main__':
    design_resistance = calculate_design_resistance()
    # The final answer is requested with the specific format.
    # The calculated value is 1349.96 kN. Rounding to one decimal place.
    print(f"\n<<< {design_resistance:.1f} >>>")