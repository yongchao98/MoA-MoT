import sympy

def solve_machin_like_formula():
    """
    Finds the integer coefficients for the given Machin-like formula for pi/4.
    """
    # Set the precision for the calculations. High precision is crucial for PSLQ to work.
    sympy.mpmath.mp.dps = 100

    # The denominators of the arctan arguments.
    y = [122, 239, 682, 1252, 2855, 12943]

    # Create the list of values for which to find an integer relation.
    # The values are arctan(1/y_i) and pi/4.
    values_to_relate = [sympy.atan(sympy.S(1)/i) for i in y] + [sympy.pi/4]

    # Use the PSLQ algorithm implementation in sympy.mpmath.identify to find the coefficients.
    # It finds an integer vector `coeffs` such that `coeffs . values_to_relate` is nearly zero.
    try:
        coeffs = sympy.mpmath.identify(values_to_relate, tol=1e-90)
    except Exception as e:
        print(f"Could not find an integer relation. The error was: {e}")
        return

    # The problem asks for the smallest *positive* n.
    # The relation is sum(c_i * arctan) + c_pi * (pi/4) = 0.
    # We want sum(c_i * arctan) = n * pi/4.
    # So, n = -c_pi.
    # If n turns out to be negative, we can multiply all coefficients by -1
    # to make n positive.
    if coeffs[-1] > 0:
        coeffs = [-c for c in coeffs]

    c = list(coeffs[:-1])
    n = -coeffs[-1]

    # Construct and print the equation string.
    equation_parts = []
    for i in range(len(y)):
        coeff = c[i]
        if coeff != 0:
            # Use '+' for positive coeffs (except the first one) and '-' for negative.
            sign = ""
            if len(equation_parts) > 0:
                sign = " + " if coeff > 0 else " - "
            elif coeff < 0:
                 sign = "-"
            
            # Format the coefficient, omitting it if it is 1 or -1.
            val = abs(coeff)
            coeff_str = str(val) if val != 1 else ""
            if val != 1 and len(equation_parts) == 0 and coeff < 0:
                 coeff_str = str(val)
                 
            term = f"{sign}{coeff_str}*arctan(1/{y[i]})"
            # Special case for first term if it's positive
            if len(equation_parts) == 0 and coeff > 0:
                 term = f"{coeff_str}*arctan(1/{y[i]})"

            equation_parts.append(term)

    # Join the parts of the sum and add the right hand side.
    lhs = " ".join(equation_parts)
    rhs = f"{n}*pi/4"
    print(f"The equation is:")
    print(f"{lhs} = {rhs}")
    
    # Print the coefficients for the final answer format
    print("\nThe coefficients are:")
    print(f"n = {n}")
    for i in range(len(c)):
        print(f"c_{i+1} = {c[i]}")

    # Prepare final answer string
    final_answer = f"<<<{n},{','.join(map(str, c))}>>>"
    print("\nFinal answer in required format:")
    print(final_answer)


solve_machin_like_formula()