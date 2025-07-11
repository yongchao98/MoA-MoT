def solve_field_theory_inner_product():
    """
    This function constructs and prints the expression for the inner product
    (ϕ, D_ϕ) for a neutral scalar field in finite-temperature field theory.
    """

    # Define the numbers that appear in the final equation
    spacetime_dimension = 4
    mass_term_exponent = 2

    # Define the mathematical symbols as strings
    integral_symbol = "∫"
    field_phi = "ϕ(x)"
    euclidean_laplacian = "□_E"
    mass_symbol = "m"

    # Construct the components of the equation using the defined numbers and symbols
    differential_measure = f"d^{spacetime_dimension}x_E"
    mass_term = f"{mass_symbol}^{mass_term_exponent}"
    
    # The operator D = (-□_E + m²)
    operator_D_acting_on_phi = f"(-{euclidean_laplacian} + {mass_term}) {field_phi}"

    # The full integrand is ϕ(x) * [D * ϕ(x)]
    integrand = f"{field_phi} {operator_D_acting_on_phi}"

    # Assemble the final equation for the inner product
    final_equation = f"{integral_symbol} {differential_measure} {integrand}"

    print("The inner product (ϕ, D_ϕ) for a free, neutral scalar field is given by the expression:")
    
    # The final print statement outputs the full equation, including the numbers 4 and 2.
    print(final_equation)

solve_field_theory_inner_product()