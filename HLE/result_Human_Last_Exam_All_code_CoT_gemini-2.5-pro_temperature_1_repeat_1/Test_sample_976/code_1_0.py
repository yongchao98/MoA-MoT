def solve_and_print_equations():
    """
    This function presents the final equations for the electric potential and field
    outside the conducting sphere, as derived from the principles of electrostatics
    in a steady state.
    """

    print("The problem is to find the electric potential and field in the region outside the sphere (r > R).")
    print("Based on solving Laplace's equation with the appropriate boundary conditions, the correct expressions are found in Answer B.")
    print("-" * 50)

    # Define the individual terms of the equations to meet the "output each number" instruction.
    # The potential is composed of the uniform field potential and the induced dipole potential.
    potential_uniform_term = "-E_0 * r * cos(θ)"
    dipole_coefficient = "(σ₁ - σ₂)*R³ / (σ₁ + 2σ₂)"
    potential_dipole_term = f"{dipole_coefficient} / r² * cos(θ)"

    # The field is composed of radial and polar components.
    # We define the key factor appearing in the field equations.
    field_factor = f"2*(σ₁ - σ₂)*R³ / ((σ₁ + 2σ₂)*r³)"
    field_factor_theta = f"(σ₁ - σ₂)*R³ / ((σ₁ + 2σ₂)*r³)"


    print("\nFinal Expression for Electric Potential Φ(r, θ) for r > R:")
    # We print the complete, properly formatted equation.
    # The numbers in the equation are 1 (implicit), 1, -1, 3, 1, 2, 2.
    final_potential_expression = r"-E_0 * (r - ( (σ₁ - σ₂) * R³ ) / ( (σ₁ + 2*σ₂) * r² ) ) * cos(θ)"
    print(final_potential_expression)
    print("\nBroken down into its parts:")
    print(f"  Φ(r, θ) = ({potential_uniform_term}) + ({potential_dipole_term})")


    print("\nFinal Expression for Electric Field E(r, θ) for r > R:")
    # We print the complete, properly formatted equation.
    # The numbers in the E_r part are 1, 2, 1, -1, 3, 1, 2, 3.
    # The numbers in the E_theta part are -1, 1, 1, -1, 3, 1, 2, 3.
    radial_component = r"E₀ * [1 + 2*(σ₁ - σ₂)*R³ / ((σ₁ + 2σ₂)*r³)] * cos(θ) * r̂"
    polar_component = r"- E₀ * [1 - (σ₁ - σ₂)*R³ / ((σ₁ + 2σ₂)*r³)] * sin(θ) * θ̂"
    print(f"E(r, θ) = ({radial_component}) + ({polar_component})")
    print("-" * 50)


solve_and_print_equations()