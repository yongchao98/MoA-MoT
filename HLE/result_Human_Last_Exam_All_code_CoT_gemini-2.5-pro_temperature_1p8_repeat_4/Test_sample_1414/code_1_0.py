from mpmath import mp

def find_machin_like_coefficients():
    """
    This function finds the integer coefficients n, c1, ..., c6 for the given
    Machin-like formula for pi using the PSLQ algorithm.
    """
    # PSLQ requires high precision arithmetic to find the correct integer relation.
    # We set the precision to 100 decimal places.
    mp.dps = 100

    # These are the denominators from the problem statement.
    denominators = [122, 239, 682, 1252, 2855, 12943]

    # Create the vector of arctangent values for PSLQ.
    # The first term corresponds to pi/4, which is arctan(1).
    vec = [mp.atan(1)] + [mp.atan(mp.mpf(1)/x) for x in denominators]

    # Use the PSLQ algorithm to find the integer relation vector [a0, a1, ...].
    # This vector satisfies a0*vec[0] + a1*vec[1] + ... = 0.
    raw_coeffs = mp.pslq(vec)

    # From the problem's structure, n*atan(1) - c1*atan(1/122) - ... = 0,
    # we can map the PSLQ result to our coefficients.
    # n = a0, c1 = -a1, c2 = -a2, and so on.
    n = raw_coeffs[0]
    coeffs = [-c for c in raw_coeffs[1:]]

    # The problem asks for the smallest positive n. If PSLQ returns a negative n,
    # we can flip the sign of all coefficients to get an equivalent relation.
    if n < 0:
        n = -n
        coeffs = [-c for c in coeffs]

    # To fulfill the requirement of outputting each number in the final equation,
    # we first reconstruct the equation in a readable format.
    equation_terms = []
    for i, coeff in enumerate(coeffs):
        if coeff != 0:
            term = f"{coeff}*arctan(1/{denominators[i]})"
            equation_terms.append(term)
    
    # Clean up the signs for printing.
    equation_str = " + ".join(equation_terms).replace(" + -", " - ")

    print("The final equation is:")
    print(f"{n}*pi/4 = {equation_str}")
    
    # Print the final list of coefficients as the result.
    c1, c2, c3, c4, c5, c6 = coeffs
    print("\nThe coefficients (n, c1, c2, c3, c4, c5, c6) are:")
    print(f"{n},{c1},{c2},{c3},{c4},{c5},{c6}")

if __name__ == "__main__":
    find_machin_like_coefficients()