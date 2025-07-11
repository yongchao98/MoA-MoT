def solve_lie_algebra_poincare():
    """
    Computes and prints the Poincaré polynomial for the given Lie algebra.
    """
    # The Betti numbers b_k for k=0 to 6
    b = [1, 3, 6, 8, 6, 3, 1]
    
    poly_terms = []
    for k, coeff in enumerate(b):
        if coeff == 0:
            continue
        
        # Format coefficient
        if coeff == 1 and k > 0:
            coeff_str = ""
        else:
            coeff_str = str(coeff)
            
        # Format variable and exponent
        if k == 0:
            var_str = ""
        elif k == 1:
            var_str = "x"
        else:
            var_str = f"x^{k}"
        
        # Combine coefficient and variable
        if k > 0 and coeff > 1:
             term = coeff_str + var_str
        elif k > 0 and coeff == 1:
             term = var_str
        else:
             term = coeff_str

        poly_terms.append(term)
        
    polynomial_string = " + ".join(poly_terms)
    
    print("The Poincaré polynomial is:")
    print(f"P(x) = {polynomial_string}")

solve_lie_algebra_poincare()