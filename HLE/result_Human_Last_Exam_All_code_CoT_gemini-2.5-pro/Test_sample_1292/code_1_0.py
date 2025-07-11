def print_field_equation():
    """
    This function prints the derived field equation for the given theory of gravity.
    Each term is printed individually to construct the full equation.
    """
    
    # Define the constant kappa^2 = 8 pi G / c^4
    kappa_squared = "\\frac{8\\pi G}{c^4}"
    
    # Define the terms of the equation
    term1 = "-\\frac{2}{\\sqrt{-g}}\\partial_{\\alpha}(\\sqrt{-g}P^\\alpha_{\\mu\\nu})"
    term2 = "- 2P_{\\mu\\alpha\\beta} Q_\\nu^{\\alpha\\beta}"
    term3 = "+ Q^{\\alpha\\beta}_\\mu P_{\\alpha\\beta\\nu}"
    term4 = "- \\frac{1}{2}Qg_{\\mu\\nu}"
    rhs = f"= {kappa_squared} T_{{\\mu\\nu}}"
    
    # Print the equation
    print("The derived field equation is:")
    print(f"{term1} {term2} {term3} {term4} {rhs}")

print_field_equation()