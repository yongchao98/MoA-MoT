def print_inner_product_expression():
    """
    This function prints the mathematical expression for the inner product
    (φ, Dφ) in finite-temperature field theory for a neutral scalar field.
    The expression includes all relevant numbers (exponents and dimensions)
    as requested.
    """
    
    # Define the components of the equation using strings with unicode for clarity
    lhs = "(φ, Dφ)"
    integral_sign = "∫"
    differential = "d⁴x_E"
    
    # Integrand terms: (∂_μ φ)² + m²φ²
    partial_mu_phi_sq = "(∂_μ φ)²"
    m_sq_phi_sq = "m²φ²"
    
    integrand = f"[{partial_mu_phi_sq} + {m_sq_phi_sq}]"
    
    # Assemble the final equation string
    full_equation = f"{lhs} = {integral_sign} {differential} {integrand}"
    
    print("In the functional integral representation for a neutral scalar field at finite temperature,")
    print("the inner product (φ, Dφ) is given by the following expression:")
    
    # Print a formatted box around the equation
    print("\n" + "=" * (len(full_equation) + 4))
    print(f"  {full_equation}  ")
    print("=" * (len(full_equation) + 4) + "\n")
    
    print("Where:")
    print("  - φ: The neutral scalar field.")
    print("  - D: The Euclidean Klein-Gordon operator, defined as D = -∂_μ² + m².")
    print("  - ∂_μ: The four-gradient in Euclidean spacetime.")
    print("  - m: The mass of the scalar field.")
    print("  - ∫ d⁴x_E: The integral over 4-dimensional Euclidean spacetime.")
    print("  - The numbers in the equation are 4 (for the spacetime dimension) and 2 (for the exponents).")

# Execute the function to print the result
print_inner_product_expression()