import math

def solve_guide_displacement():
    """
    Calculates the horizontal displacement of the guide.
    """
    # Define the constants from the problem in SI units
    m = 0.20  # kg, mass of the body
    M = 0.80  # kg, mass of the guide
    R = 0.20  # m, radius of the circular arcs
    d = 0.50  # m, length of the straight section
    mu_D = 0.20 # coefficient of friction

    # --- Calculation Steps ---

    # The core principle is the conservation of the center of mass, leading to the formula:
    # Δx_M = - (m / (m + M)) * Δx_m_guide
    # where Δx_m_guide is the horizontal displacement of the mass relative to the guide.

    # To find Δx_m_guide, we first calculate the final height 'h' using energy conservation.
    # h = R - μ_D * d
    # Then, we find the final x-position of the mass relative to the guide's frame.
    # The displacement of the mass relative to the guide is: Δx_m_guide = d - sqrt(R^2 - (R - h)^2)
    # This simplifies to:
    # Δx_m_guide = d - sqrt(R^2 - (μ_D * d)^2)

    # Calculate the term (μ_D * d)^2
    mu_d_sq = (mu_D * d)**2

    # Calculate the term inside the square root: R^2 - (μ_D * d)^2
    term_under_sqrt = R**2 - mu_d_sq

    # Calculate the horizontal displacement of the mass relative to the guide
    delta_x_m_guide = d - math.sqrt(term_under_sqrt)

    # Calculate the total mass
    total_mass = m + M

    # Finally, calculate the horizontal displacement of the guide
    displacement_M = - (m / total_mass) * delta_x_m_guide

    # --- Output the step-by-step equation as requested ---
    print("The final equation for the displacement of the guide (Δx_M) is:")
    print("Δx_M = - (m / (m + M)) * (d - sqrt(R^2 - (μ_D * d)^2))")
    
    print("\nSubstituting the values into the equation:")
    print(f"Δx_M = - ({m} / ({m} + {M})) * ({d} - sqrt({R}^2 - ({mu_D} * {d})^2))")
    print(f"Δx_M = - ({m / total_mass:.2f}) * ({d} - sqrt({R**2:.2f} - {mu_d_sq:.2f}))")
    print(f"Δx_M = - ({m / total_mass:.2f}) * ({d} - sqrt({term_under_sqrt:.2f}))")
    print(f"Δx_M = - ({m / total_mass:.2f}) * ({d} - {math.sqrt(term_under_sqrt)})")
    print(f"Δx_M = - ({m / total_mass:.2f}) * ({delta_x_m_guide})")
    print(f"Δx_M = {displacement_M}")

    # The negative sign indicates the guide moves to the left.
    # The final answer is the numerical value of the displacement.
    print(f"\nFinal Answer (in meters): {displacement_M}")
    
    # Return the final numerical value in the specified format
    return f"<<<{displacement_M}>>>"

# Execute the function and print the final result
final_answer = solve_guide_displacement()
print(final_answer)