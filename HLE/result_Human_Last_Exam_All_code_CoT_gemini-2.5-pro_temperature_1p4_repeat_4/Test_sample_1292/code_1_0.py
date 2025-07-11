def print_field_equation():
    """
    This function prints the derived field equation for the given theory of gravity.
    Each term is represented symbolically.
    """
    kappa = "8pi G / c^4"
    T_munu = "T_{\u03BC\u03BD}"

    # Term 1: Divergence of the superpotential
    term1_coeff = "-2/\u221A(-g)"
    term1_expr = "\u2202\u2090(\u221A(-g)P\u207D\u2090\u2098\u2099)"
    term1 = f"{term1_coeff}{term1_expr}"

    # Term 2: Quadratic in non-metricity
    term2 = "- P_{\u03BC\u03B1\u03B2} Q\u2099\u207D\u1D43\u1D45"

    # Term 3: Quadratic in non-metricity
    term3 = "+ 2Q\u207D\u1D43\u1D45\u2098 P_{\u03B1\u03B2\u03BD}"
    
    # Term 4: The non-metricity scalar term
    term4 = "+ (1/2)Q g_{\u03BC\u03BD}"

    # Combine all terms for the left-hand side (LHS)
    lhs = f"{term1} {term2} {term3} {term4}"

    # Right-hand side (RHS)
    rhs = f"({kappa}) {T_munu}"
    
    # Print the full equation
    print(f"The derived field equation is:")
    print(f"{lhs} = {rhs}")
    
    # Print each term individually as requested
    print("\nBreaking down the equation:")
    print(f"Term 1: {term1}")
    print(f"Term 2: {term2}")
    print(f"Term 3: {term3}")
    print(f"Term 4: {term4}")
    print(f"RHS: {rhs}")
    
    print("\nFinal equation in plain text:")
    final_eq_text = f"(-2/sqrt(-g)) * d_alpha(sqrt(-g) * P^alpha_munu) - P_mualphabeta * Q_nu^alphabeta + 2 * Q^alphabeta_mu * P_alphabetanu + (1/2) * Q * g_munu = (8*pi*G/c^4) * T_munu"
    print(final_eq_text)
    
print_field_equation()