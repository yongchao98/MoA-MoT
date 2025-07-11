def display_infinite_product_formula():
    """
    This function generates and displays the closed-form expression for the infinite product:
    Π_{n=3 to ∞} (1 - z³/n³)
    """
    # The derivation relies on the Weierstrass product for the Gamma function, which gives:
    # Π_{n=1 to ∞} (1 - z³/n³) = 1 / [Γ(1-z)Γ(1-ωz)Γ(1-ω²z)]
    # where ω = exp(2πi/3) is a primitive cube root of unity.

    # The product starting from n=3 is obtained by dividing the product from n=1
    # by the terms for n=1 and n=2.
    # Term n=1: (1 - z³/1³) = 1 - z³
    # Term n=2: (1 - z³/2³) = 1 - z³/8

    # Combining these gives the final expression.

    # Using unicode for mathematical symbols for better readability
    # Π = \u03a0, ∞ = \u221e, Γ = \u0393, ω = \u03c9, ² = \u00b2, ³ = \u00b3
    
    product_str = "\u03a0_{n=3}^\u221e (1 - z\u00b3/n\u00b3)"
    
    numerator = "1"
    
    # The denominator shows all the numbers from the formula, including 3 (power),
    # 8 (from 2^3), and implies the starting n=3 by the context.
    denominator = "(1 - z\u00b3) * (1 - z\u00b3/8) * \u0393(1-z)\u0393(1-\u03c9z)\u0393(1-\u03c9\u00b2z)"
    
    print(f"The closed-form expression for {product_str} is:")
    print()
    
    # Pretty print the fraction
    max_len = len(denominator)
    print(f" {numerator.center(max_len)} ")
    print(f" {'-' * (max_len + 2)} ")
    print(f" {denominator} ")
    print()
    print("where:")
    print("  \u0393(x) is the Gamma function.")
    print("  \u03c9 = exp(2\u03c0i/3) is a primitive cube root of unity.")

display_infinite_product_formula()