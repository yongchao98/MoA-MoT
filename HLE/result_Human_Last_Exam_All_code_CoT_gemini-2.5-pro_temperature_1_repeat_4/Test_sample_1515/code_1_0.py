def print_nsvz_beta_function():
    """
    This function prints the formula for the NSVZ beta function.
    The formula relates the beta function of the gauge coupling 'α' to the anomalous
    dimensions 'γ' of the matter fields.
    """

    # Define the components of the equation as string variables for clarity.
    beta_alpha = "β(α)"
    alpha = "α"
    numerator = "3*T(G) - Σ_i T(r_i) * (1 - γ_i)"
    denominator = "1 - (T(G) * g**2) / (8 * π**2)"

    # Assemble the final equation string.
    # We represent the fraction with a placeholder to be explained below.
    # The NSVZ relation is often written in terms of g or α = g^2 / (4π).
    # Here we present the form involving the anomalous dimension explicitly.
    equation_part1 = f"{beta_alpha} = - (g**3 / (16 * π**2))"
    equation_part2 = f"[({numerator}) / ({denominator})]"

    print("The exact NSVZ beta function is given by the following relation:")
    print("-" * 60)
    # We print each part of the final equation.
    print(f"{equation_part1} * {equation_part2}")
    print("-" * 60)
    print("Where:")
    print("  β(α) is the beta function for the gauge coupling α = g^2 / (4π).")
    print("  g is the gauge coupling constant.")
    print("  π is the mathematical constant Pi.")
    print("  T(G) is the Dynkin index for the adjoint representation of the gauge group.")
    print("  T(r_i) is the Dynkin index for the matter field representation r_i.")
    print("  γ_i is the anomalous dimension of the i-th matter superfield.")
    print("\nThis relation holds in a specific regularization scheme (the 'NSVZ scheme')")
    print("that preserves the holomorphy required by non-renormalization theorems.")

# Execute the function to print the output.
print_nsvz_beta_function()