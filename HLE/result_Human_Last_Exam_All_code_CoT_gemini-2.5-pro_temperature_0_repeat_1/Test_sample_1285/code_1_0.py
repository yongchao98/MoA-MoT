def solve_wronskian_roots():
    """
    This function calculates the maximum number of roots for the Wronskian R_t.
    """
    
    # Part (a): Find the expression for the maximum number of roots for a general t.
    # The maximum number of roots of R_t in the interval ]0, 1[ is given by
    # the degree of the polynomial part of the Wronskian of the functions f_i.
    # For functions of the form x^k * (1-x)^l, this degree is known to be at most t(t-1)/2.
    # This bound is sharp, meaning we can find exponents k_i and l_i for which
    # this number of roots is achieved.
    # The formula is the binomial coefficient "t choose 2".
    
    formula_a = "t * (t - 1) / 2"
    print(f"(a) The expression for the maximum number of roots is: {formula_a}")
    
    # Part (b): Calculate the maximum number of roots for t = 5.
    t = 5
    
    # Calculate the result using the formula.
    numerator = t * (t - 1)
    result_b = numerator / 2
    
    # Print the calculation step-by-step as requested.
    print(f"\n(b) For t = {t}, we substitute t into the formula:")
    print(f"    Maximum roots = {t} * ({t} - 1) / 2")
    print(f"                  = {t} * {t - 1} / 2")
    print(f"                  = {numerator} / 2")
    print(f"                  = {int(result_b)}")
    print(f"\nTherefore, for t = 5, the maximum number of roots is {int(result_b)}.")

solve_wronskian_roots()