import mpmath

def find_machin_like_coefficients():
    """
    Solves for the integer coefficients n, c1, ..., c6 in the given
    Machin-like formula using the PSLQ algorithm.
    """
    # Set the precision for mpmath. PSLQ requires high precision.
    mpmath.mp.dps = 100

    # The denominators from the arctan terms in the formula
    denominators = [122, 239, 682, 1252, 2855, 12943]

    # We are looking for an integer relation for the equation:
    # n*(pi/4) - c1*arctan(1/122) - c2*arctan(1/239) - ... = 0
    # The input vector for PSLQ will be [pi/4, arctan(1/122), arctan(1/239), ...].
    terms_to_relate = [mpmath.pi / 4] + [mpmath.atan(mpmath.mpf(1) / d) for d in denominators]

    # Use the PSLQ algorithm to find the integer coefficients [k0, k1, k2, ...]
    # such that k0*x0 + k1*x1 + k2*x2 + ... = 0.
    found_coeffs = mpmath.pslq(terms_to_relate)

    # The problem asks for the smallest positive n.
    # The first coefficient, k0, corresponds to n. If it's negative,
    # we flip the signs of all coefficients to make n positive.
    if found_coeffs[0] < 0:
        found_coeffs = [-c for c in found_coeffs]

    # Map the PSLQ result to the variables in the problem.
    # n = k0
    # ci = -ki
    n = found_coeffs[0]
    c = [-k for k in found_coeffs[1:]]

    # Print the full equation with the found coefficients for clarity.
    print("The solved equation is:")
    print(f"{n} * pi/4 = ", end="")
    
    is_first_term = True
    for i in range(len(c)):
        coeff = c[i]
        denom = denominators[i]
        
        # Add a sign before the term, unless it's the first one.
        if not is_first_term:
            if coeff >= 0:
                print(" + ", end="")
            else:
                print(" - ", end="")
                coeff = abs(coeff)
        # Handle a negative coefficient for the very first term.
        elif coeff < 0:
            print("-", end="")
            coeff = abs(coeff)
            
        print(f"{coeff}*arctan(1/{denom})", end="")
        is_first_term = False
    print("\n")

    # Print the coefficients individually.
    print("The coefficients are:")
    print(f"n = {n}")
    for i, val in enumerate(c, 1):
        print(f"c{i} = {val}")
        
    # Prepare the final answer in the specified format for the platform.
    final_answer_list = [n] + c
    final_answer_str = ", ".join(map(str, final_answer_list))
    # print(f"\n<<<{[n] + c}>>>") # For verification

if __name__ == '__main__':
    find_machin_like_coefficients()