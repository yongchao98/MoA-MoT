import mpmath

def solve_machin_like_formula():
    """
    Finds the integer coefficients for the given Machin-like formula
    using the PSLQ algorithm.
    """
    # Set a high precision for the calculations, as required by PSLQ.
    mpmath.mp.dps = 150

    # The denominators of the arctan arguments.
    xs = [122, 239, 682, 1252, 2855, 12943]

    # Create the vector of real values for which we want to find an integer relation.
    # The relation has the form: n*(pi/4) - c1*atan(1/x1) - c2*atan(1/x2) - ... = 0
    vals = [mpmath.pi / 4] + [mpmath.atan(mpmath.mpf(1) / x) for x in xs]

    # Use the PSLQ algorithm to find the integer coefficients.
    # The result from pslq corresponds to [n, -c1, -c2, ...]
    coeffs = mpmath.pslq(vals)

    # PSLQ returns a primitive relation (gcd of coefficients is 1).
    # We need to ensure that 'n' (the first coefficient) is positive.
    # If it's negative, we flip the signs of all coefficients.
    if coeffs[0] < 0:
        coeffs = [-c for c in coeffs]

    n = coeffs[0]
    # The remaining coefficients from pslq are -c_i, so we flip their signs back.
    cs = [-c for c in coeffs[1:]]

    # Output the discovered equation
    print("The discovered identity is:")
    
    equation_parts = [f"{n}*π/4 ="]
    is_first_term = True
    for i, (c, x) in enumerate(zip(cs, xs)):
        if c == 0:
            continue
        
        sign = ""
        # Add a '+' or '-' sign for terms after the first one.
        if not is_first_term:
            sign = " + " if c > 0 else " - "
        # For a negative first term, prepend a '-'
        elif c < 0:
            sign = "-"
        
        abs_c = abs(c)
        
        # Format the term, omitting the coefficient if it's 1.
        if abs_c == 1:
            term = f"arctan(1/{x})"
        else:
            term = f"{abs_c}*arctan(1/{x})"
        
        equation_parts.append(f"{sign}{term}")
        is_first_term = False
    
    # Print the full equation string.
    print("".join(equation_parts))
    
    # Print the coefficients individually as requested.
    print("\nThe solution for the constants is:")
    print(f"n = {n}")
    for i, c in enumerate(cs, 1):
        print(f"c{i} = {c}")

    # For the final answer format
    final_coeffs = [n] + cs
    final_answer = ",".join(map(str, final_coeffs))
    print("\nAnswer in the required format n,c1,c2,c3,c4,c5,c6:")
    print(final_answer)
    
if __name__ == "__main__":
    solve_machin_like_formula()