import sympy

def find_machin_like_coefficients():
    """
    Solves for the integer coefficients n, c1, ..., c6 in the given Machin-like formula
    using the PSLQ algorithm.
    """
    
    # The x-values from the arctan terms in the formula
    x_vals = [122, 239, 682, 1252, 2855, 12943]

    # PSLQ requires high precision to find correct integer relations.
    # We set the precision for mpmath, which sympy uses for floating point calculations.
    sympy.mpmath.mp.dps = 150

    # Create the list of high-precision values for which to find a relation.
    # The vector is [arctan(1/x1), ..., arctan(1/x6), pi/4]
    values_to_relate = [sympy.mpmath.atan(sympy.mpmath.mpf(1)/x) for x in x_vals]
    values_to_relate.append(sympy.mpmath.pi / 4)

    # Use the PSLQ algorithm to find the list of integer coefficients.
    coeffs = sympy.pslq(values_to_relate)

    if not coeffs:
        print("An integer relation could not be found with the current precision.")
        print("You may need to increase sympy.mpmath.mp.dps.")
        return

    # The PSLQ algorithm finds coeffs [k1, ..., k7] such that:
    # k1*v1 + ... + k6*v6 + k7*v7 = 0
    # k1*arctan(1/x1) + ... + k6*arctan(1/x6) + k7*pi/4 = 0
    # Rearranging this gives: k1*arctan(1/x1) + ... = -k7*pi/4
    # Comparing with n*pi/4 = c1*arctan(1/x1) + ...
    # we find that c_i = k_i and n = -k_7.
    
    c_list = coeffs[:-1]
    n = -coeffs[-1]

    # We need the smallest positive n. If n is negative, we flip the signs of all coeffs.
    if n < 0:
        n = -n
        c_list = [-c for c in c_list]

    # Print the final equation with all the numbers.
    print("The solved equation is:")
    
    rhs_parts = []
    for i in range(len(c_list)):
        c = c_list[i]
        if c == 0:
            continue
        
        sign = "+" if c > 0 else "-"
        # We don't need to print '1' for coefficient 1, but the problem does not specify it.
        # Let's print all coefficients for clarity.
        coeff_str = str(abs(c))
        
        rhs_parts.append(f"{sign} {coeff_str}*arctan(1/{x_vals[i]})")

    # Join the parts of the right-hand side of the equation and clean up the leading sign.
    equation_rhs = " ".join(rhs_parts).lstrip('+ ')

    print(f"{n}*pi/4 = {equation_rhs}")
    
    # For clarity, let's print each constant as well
    print("\nThe constants are:")
    print(f"n = {n}")
    print(f"c1 = {c_list[0]}")
    print(f"c2 = {c_list[1]}")
    print(f"c3 = {c_list[2]}")
    print(f"c4 = {c_list[3]}")
    print(f"c5 = {c_list[4]}")
    print(f"c6 = {c_list[5]}")


# Execute the function
find_machin_like_coefficients()