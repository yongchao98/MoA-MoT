import mpmath

def find_machin_like_coefficients():
    """
    This function finds the integer coefficients n, c1, ..., c6 for the given
    Machin-like formula for pi.
    """
    # Set the precision for the calculations. A high precision is necessary
    # for the PSLQ algorithm to work correctly.
    mpmath.mp.dps = 100

    # The denominators of the arguments to the arctan functions
    denominators = [122, 239, 682, 1252, 2855, 12943]

    # The problem can be formulated as finding an integer relation for a vector of real numbers.
    # The equation is:
    # n*(pi/4) - c1*arctan(1/122) - c2*arctan(1/239) - ... = 0
    # The vector of values is [pi/4, arctan(1/122), arctan(1/239), ...].
    values = [mpmath.pi / 4]
    for d in denominators:
        values.append(mpmath.atan(mpmath.mpf(1) / d))

    # Use the PSLQ algorithm to find a list of integer coefficients [a0, a1, ..., a6]
    # such that a0*v0 + a1*v1 + ... + a6*v6 = 0.
    coeffs = mpmath.pslq(values)

    # The problem requires the smallest positive n. If PSLQ returns a set of
    # coefficients where the first one (our n) is negative, we can multiply
    # the entire set by -1 to get an equivalent relation with a positive n.
    if coeffs and coeffs[0] < 0:
        coeffs = [-c for c in coeffs]

    # Extract n and the c_k coefficients from the result of PSLQ.
    # n = a0
    # ck = -a_{k+1}
    n = coeffs[0]
    c_coeffs = [-k for k in coeffs[1:]]

    # Print the full equation as requested.
    print("The discovered equation is:")
    print(f"{n} * (pi/4) = ", end="")
    
    is_first_term = True
    for i in range(len(c_coeffs)):
        coeff = c_coeffs[i]
        denom = denominators[i]

        if coeff == 0:
            continue

        # Determine the sign to print (+ or -).
        sign = ""
        if not is_first_term:
            sign = " + " if coeff > 0 else " - "
        elif coeff < 0:
            sign = "-"
        
        print(sign, end="")

        # Print the coefficient if it's not 1 or -1.
        abs_coeff = abs(coeff)
        if abs_coeff != 1:
            print(f"{abs_coeff}*arctan(1/{denom})", end="")
        else:
            print(f"arctan(1/{denom})", end="")
        
        is_first_term = False
    print("\n")

    # Print the final answer in the required comma-separated format.
    print("The solution for (n, c1, c2, c3, c4, c5, c6) is:")
    result_string = f"{n},{','.join(map(str, c_coeffs))}"
    print(result_string)

find_machin_like_coefficients()