import mpmath

def solve_arctan_relation():
    """
    Finds the integer coefficients for the given Machin-like formula for pi
    using the PSLQ algorithm.
    """
    # Set a high precision for the calculations to ensure PSLQ works correctly.
    mpmath.mp.dps = 100

    # The arguments for the arctan functions in the equation.
    # The first one is 1, corresponding to arctan(1) = pi/4.
    x_inv = [1, 1/122, 1/239, 1/682, 1/1252, 1/2855, 1/12943]

    # Create the vector of high-precision arctan values.
    try:
        atan_values = [mpmath.atan(mpmath.mpf(v)) for v in x_inv]
    except Exception as e:
        print(f"An error occurred during mpmath calculation: {e}")
        return

    # Use the PSLQ algorithm to find an integer relation among the arctan values.
    # It finds a list of integers `coeffs` such that sum(coeffs[i] * atan_values[i]) is close to 0.
    coeffs = mpmath.pslq(atan_values)

    if not coeffs:
        print("PSLQ algorithm did not find a relation.")
        return

    # The relation is of the form:
    # k0*atan(1) + k1*atan(1/122) + k2*atan(1/239) + ... = 0
    # We want to match it to:
    # n*atan(1) - c1*atan(1/122) - c2*atan(1/239) - ... = 0
    # So, n = k0 and ci = -ki for i > 0.
    n = coeffs[0]
    c_coeffs = [-k for k in coeffs[1:]]

    # We need the smallest positive n. If n is negative, flip all signs.
    if n < 0:
        n = -n
        c_coeffs = [-c for c in c_coeffs]

    c1, c2, c3, c4, c5, c6 = c_coeffs
    x_vals = [122, 239, 682, 1252, 2855, 12943]

    # Print the resulting coefficients
    print(f"Found coefficients:")
    print(f"n = {n}")
    print(f"c1 = {c1}")
    print(f"c2 = {c2}")
    print(f"c3 = {c3}")
    print(f"c4 = {c4}")
    print(f"c5 = {c5}")
    print(f"c6 = {c6}")
    
    # Print the final equation with each number explicitly shown
    print("\nThe final equation is:")
    lhs = f"{n} * pi/4"
    rhs_parts = []
    for c, x in zip(c_coeffs, x_vals):
        # We check c!=0 to not include terms with a zero coefficient
        if c != 0:
            rhs_parts.append(f"({c})*arctan(1/{x})")

    rhs = " + ".join(rhs_parts)
    print(f"{lhs} = {rhs}")
    
    # Print the final answer in the requested format for parsing
    final_answer = f"<<<{n},{c1},{c2},{c3},{c4},{c5},{c6}>>>"
    # This line will not be printed to the user console, but it's part of the standard output.
    # We add an extra print statement to ensure it's captured.
    print(f"\nFinal answer string for submission:")
    print(final_answer)


if __name__ == '__main__':
    solve_arctan_relation()