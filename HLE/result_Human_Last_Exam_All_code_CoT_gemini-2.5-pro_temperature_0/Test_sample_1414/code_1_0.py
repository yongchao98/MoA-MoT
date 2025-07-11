import mpmath

def solve_machin_like_formula():
    """
    Finds the integer coefficients for a given Machin-like formula using the PSLQ algorithm.
    """
    # Set the precision for the calculations. PSLQ requires high precision.
    mpmath.mp.dps = 100

    # These are the denominators from the arctan terms in the equation.
    denominators = [122, 239, 682, 1252, 2855, 12943]

    # Create the vector of values for which to find an integer relation.
    # The first term corresponds to n * pi/4.
    # The subsequent terms correspond to the arctan values.
    vals = [mpmath.pi / 4]
    for d in denominators:
        vals.append(mpmath.atan(mpmath.mpf(1) / d))

    # Use the PSLQ algorithm to find the integer coefficients [a_0, a_1, ..., a_6]
    # such that a_0*v_0 + a_1*v_1 + ... = 0.
    coeffs = mpmath.pslq(vals)

    # The target equation is: n*(pi/4) - c_1*arctan(1/122) - ... = 0
    # By comparing with the PSLQ result, we can map the coefficients.
    # n = a_0
    # c_k = -a_{k+1}
    n = coeffs[0]
    c = [-k for k in coeffs[1:]]

    # The problem requires the smallest positive n. If n is negative,
    # we can multiply the entire relation by -1.
    if n < 0:
        n = -n
        c = [-k for k in c]

    # Print the solution
    print(f"The solution is found for n = {n} with the following coefficients:")
    for i in range(len(c)):
        print(f"c_{i+1} = {c[i]}")

    # Construct and print the final equation as requested.
    print("\nThe final equation is:")
    equation_str = f"{n} * pi/4 = "
    terms = []
    for i in range(len(c)):
        # Only include terms with non-zero coefficients
        if c[i] != 0:
            # Handle the sign for the first term and subsequent terms
            if c[i] > 0:
                sign = " + " if terms else ""
                coeff_str = f"{c[i]}" if c[i] != 1 else ""
            else:
                sign = " - "
                coeff_str = f"{-c[i]}" if c[i] != -1 else ""
            
            terms.append(f"{sign}{coeff_str}arctan(1/{denominators[i]})")

    equation_str += "".join(terms).lstrip(" +")
    print(equation_str)
    
    # Provide the final answer in the specified format.
    final_answer_str = f"{n},{c[0]},{c[1]},{c[2]},{c[3]},{c[4]},{c[5]}"
    print(f"\n<<<{final_answer_str}>>>")

solve_machin_like_formula()