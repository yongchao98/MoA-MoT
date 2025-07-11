def solve_infinite_product():
    """
    This function provides the symbolic solution to the given infinite product.
    The solution is derived using the Weierstrass product representation of the Gamma function.
    """
    # The infinite product from n=1 to infinity is 1 / (Gamma(1-z) * Gamma(1-z*omega) * Gamma(1-z*omega^2))
    # where omega = exp(i*2*pi/3).
    # The requested product starts from n=3, so we divide by the terms for n=1 and n=2.
    
    # Construct the denominator parts
    n1_term = "(1 - z^3/1^3)"
    n2_term = "(1 - z^3/2^3)"
    gamma_z = "Gamma(1 - z)"
    gamma_z_omega = "Gamma(1 - z*exp(i*2*pi/3))"
    gamma_z_omega2 = "Gamma(1 - z*exp(i*4*pi/3))"
    
    # Assemble the final expression
    # Note that all numbers from the original expression (3, 1, 2) are present.
    final_expression = f"1 / ({n1_term} * {n2_term} * {gamma_z} * {gamma_z_omega} * {gamma_z_omega2})"
    
    print(final_expression)

solve_infinite_product()