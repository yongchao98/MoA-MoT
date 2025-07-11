import mpmath
import math

def solve_machin_like_formula():
    """
    Finds the integer coefficients for the given Machin-like formula for pi.
    """
    # Set the precision for the calculations. 100 digits is usually sufficient for PSLQ.
    mpmath.mp.dps = 100

    # These are the denominators from the problem statement.
    denominators = [122, 239, 682, 1252, 2855, 12943]

    # Create a list of the high-precision values for the integer relation search.
    # The list contains the arctan values and the value pi/4.
    values = [mpmath.atan(mpmath.mpf(1)/d) for d in denominators]
    values.append(mpmath.pi() / 4)

    # Use the PSLQ algorithm to find a list of integer coefficients `coeffs`
    # such that: coeffs[0]*values[0] + ... + coeffs[6]*values[6] = 0
    try:
        coeffs = mpmath.pslq(values)
    except Exception as e:
        print(f"An error occurred during PSLQ calculation: {e}")
        print("This might happen if the precision (dps) is too low.")
        return

    # The equation is given as:
    # c1*atan(1/x1) + ... = n * pi/4
    # Our relation found by PSLQ is:
    # coeffs[0]*atan(1/x1) + ... + coeffs[5]*atan(1/x6) + coeffs[6]*(pi/4) = 0
    # So, we can set c_k = coeffs[k] and n = -coeffs[6].
    c_candidates = coeffs[:-1]
    n_candidate = -coeffs[-1]

    # The problem asks for the smallest POSITIVE n.
    # If the n found is negative, we can multiply all coefficients by -1.
    if n_candidate < 0:
        n = -n_candidate
        c = [-x for x in c_candidates]
    else:
        n = n_candidate
        c = c_candidates

    # To ensure we have the smallest integers, we can divide by the GCD.
    # PSLQ usually returns a primitive relation, but this is a safeguard.
    all_coeffs = [n] + c
    
    # Python's math.gcd only takes two arguments. We can reduce the list.
    common_divisor = all_coeffs[0]
    for i in range(1, len(all_coeffs)):
        common_divisor = math.gcd(common_divisor, all_coeffs[i])

    if common_divisor > 1:
        n = n // common_divisor
        c = [x // common_divisor for x in c]

    # Now, we print the results in the required format.
    # First, print the full equation.
    equation_str = f"{n} * pi / 4 = "
    first_term = True
    for i, coeff in enumerate(c):
        if coeff == 0:
            continue

        # Determine the sign for the term.
        sign = ""
        if not first_term:
            if coeff > 0:
                sign = " + "
            else:
                sign = " - "
        elif coeff < 0:
            sign = "-"
        
        abs_coeff = abs(coeff)

        if abs_coeff == 1:
            term_str = f"arctan(1/{denominators[i]})"
        else:
            term_str = f"{abs_coeff}*arctan(1/{denominators[i]})"
        
        equation_str += sign + term_str
        first_term = False

    print("The final equation is:")
    print(equation_str)
    print("-" * 20)

    # Then, provide the answer as requested.
    result_list = [n] + c
    print("The coefficients (n, c1, c2, c3, c4, c5, c6) are:")
    print(",".join(map(str, result_list)))

# Execute the solver function.
solve_machin_like_formula()