def solve_poincare_polynomial():
    """
    This function constructs and prints the Poincaré polynomial for the given Lie algebra.
    The Betti numbers are pre-calculated based on the mathematical derivation.
    b_k = dim H^k(g, R)
    b_0 = 1
    b_1 = 3
    b_2 = 6
    b_3 = 8
    b_4 = 6 (by Poincare duality b_4 = b_2)
    b_5 = 3 (by Poincare duality b_5 = b_1)
    b_6 = 1 (by Poincare duality b_6 = b_0)
    """
    
    # The Betti numbers for the Lie algebra
    betti_numbers = [1, 3, 6, 8, 6, 3, 1]
    
    # The variable for the polynomial
    var = 'x'
    
    # Construct the polynomial string
    poly_terms = []
    
    # Zeroth term
    if betti_numbers[0] != 0:
        poly_terms.append(f"{betti_numbers[0]}")
    
    # First degree term
    if betti_numbers[1] != 0:
        if betti_numbers[1] == 1:
            poly_terms.append(f"{var}")
        else:
            poly_terms.append(f"{betti_numbers[1]}{var}")

    # Higher degree terms
    for k in range(2, len(betti_numbers)):
        if betti_numbers[k] != 0:
            if betti_numbers[k] == 1:
                 poly_terms.append(f"{var}^{k}")
            else:
                 poly_terms.append(f"{betti_numbers[k]}{var}^{k}")
    
    polynomial_string = " + ".join(poly_terms)
    
    print("The Poincaré polynomial is:")
    print(f"P(x) = {polynomial_string}")

solve_poincare_polynomial()