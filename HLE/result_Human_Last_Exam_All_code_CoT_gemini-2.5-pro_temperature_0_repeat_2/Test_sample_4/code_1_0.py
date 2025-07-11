def solve_poincare_polynomial():
    """
    This function calculates and prints the Poincaré polynomial for the given Lie algebra.
    The Betti numbers are pre-calculated based on the analysis of the
    Chevalley-Eilenberg complex associated with the Lie algebra.
    """
    # Betti numbers b_k = dim(H_k(g)) for k = 0, ..., 6
    betti_numbers = [1, 3, 6, 8, 6, 3, 1]
    
    polynomial_terms = []
    
    # Handle the constant term (k=0)
    if betti_numbers[0] != 0:
        polynomial_terms.append(str(betti_numbers[0]))
        
    # Handle the term for k=1
    if betti_numbers[1] != 0:
        if betti_numbers[1] == 1:
            polynomial_terms.append("x")
        else:
            polynomial_terms.append(f"{betti_numbers[1]}*x")
            
    # Handle terms for k > 1
    for k in range(2, len(betti_numbers)):
        if betti_numbers[k] != 0:
            if betti_numbers[k] == 1:
                polynomial_terms.append(f"x^{k}")
            else:
                polynomial_terms.append(f"{betti_numbers[k]}*x^{k}")
                
    # Join the terms to form the polynomial string
    polynomial_string = " + ".join(polynomial_terms)
    
    print("The Poincaré polynomial is:")
    print(f"P(x) = {polynomial_string}")

solve_poincare_polynomial()