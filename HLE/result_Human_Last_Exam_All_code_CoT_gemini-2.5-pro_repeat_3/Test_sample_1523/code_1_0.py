import sympy

def get_inner_product_expression():
    """
    Symbolically computes and displays the inner product (phi, D*phi)
    in finite-temperature scalar field theory.
    """
    # Define symbols and functions
    m = sympy.Symbol('m')
    tau = sympy.Symbol('τ')
    x1, x2, x3 = sympy.symbols('x1 x2 x3')
    phi = sympy.Function('ϕ')(tau, x1, x2, x3)

    # Define the Euclidean d'Alembertian operator acting on phi
    # Note: ∂_μ² = ∂_τ² + ∇²
    d2_dtau2 = sympy.Derivative(phi, tau, 2)
    d2_dx1_2 = sympy.Derivative(phi, x1, 2)
    d2_dx2_2 = sympy.Derivative(phi, x2, 2)
    d2_dx3_2 = sympy.Derivative(phi, x3, 2)
    box_E_phi = d2_dtau2 + d2_dx1_2 + d2_dx2_2 + d2_dx3_2

    # Define the Euclidean Klein-Gordon operator D acting on phi
    # D = -∂_μ² + m²
    D_phi = -box_E_phi + m**2 * phi

    # The inner product (ϕ, Dϕ) is the integral of ϕ * (Dϕ) over spacetime.
    # We will print the full expression symbolically.
    
    # Using unicode characters for a clearer mathematical representation
    integral_symbol = "∫"
    partial_symbol = "∂"
    
    # Constructing the string for the final equation
    # The term (-∂_μ² + m²)ϕ
    operator_term = f"(-{partial_symbol}_τ² - ∇² + m²)ϕ"
    
    # The full inner product equation
    final_equation = f"(ϕ, Dϕ) = {integral_symbol} dτ d³x  [ ϕ * ({operator_term}) ]"
    
    print("In finite-temperature field theory for a neutral scalar field, the operator D in the quadratic action S = ½(ϕ, Dϕ) is the Euclidean Klein-Gordon operator.")
    print("The inner product (ϕ, Dϕ) is therefore given by the following expression:")
    print("-" * 50)
    print(final_equation)
    print("-" * 50)
    
    # As requested by the prompt, we also output each component of the equation
    # from the symbolic calculation.
    # The term inside the brackets is the integrand ϕ * Dϕ
    integrand = phi * D_phi
    print("\nSymbolic representation of the integrand ϕ * (-∂_μ² + m²)ϕ:")
    sympy.pprint(integrand, use_unicode=True)


if __name__ == '__main__':
    get_inner_product_expression()
