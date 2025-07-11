def solve_infinite_product():
    """
    This function prints the closed-form expression for the given infinite product.
    The derivation relies on the Weierstrass product representation of the Gamma function.
    """
    
    # Define the components of the formula to emphasize the numbers involved.
    start_index = 3
    power = 3
    term_n1_denom = 1**power
    term_n2_denom = 2**power # which is 8
    
    # Construct the final formula string using unicode characters for clarity.
    lhs = f"Π_{{n={start_index}..∞}} (1 - z³ / n³)"
    
    # The denominator consists of the terms for n=1 and n=2, and the product of Gamma functions.
    # We write the explicit simplified form.
    term1 = f"(1 - z³)"
    term2 = f"(1 - z³/{term_n2_denom})"
    gamma_terms = "Γ(1-z) Γ(1-ωz) Γ(1-ω²z)"
    
    rhs = f"1 / ({term1} * {term2} * {gamma_terms})"

    print(f"{lhs} = {rhs}")
    print("\nwhere:")
    print("  • ω = exp(2πi/3) is the principal cube root of unity.")
    print("  • Γ(z) is the Gamma function.")

solve_infinite_product()