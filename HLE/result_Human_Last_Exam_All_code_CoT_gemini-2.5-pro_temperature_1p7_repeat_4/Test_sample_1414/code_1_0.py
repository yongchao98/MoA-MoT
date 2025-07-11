import mpmath

def solve_machin_like_formula():
    """
    Finds the integer coefficients for the given Machin-like formula for pi
    using the PSLQ algorithm.
    """
    # Set the precision for the calculations. High precision is necessary for PSLQ to work.
    mpmath.mp.dps = 200

    # Define the values for which we want to find an integer relation.
    # The first value corresponds to n * pi/4.
    # The other values correspond to the arctan terms.
    args = [122, 239, 682, 1252, 2855, 12943]
    values = [mpmath.pi / 4] + [mpmath.atan(mpmath.mpf(1)/x) for x in args]

    # Use the PSLQ algorithm to find the list of integer coefficients.
    coeffs = mpmath.pslq(values)

    # The first coefficient corresponds to 'n'.
    n = coeffs[0]
    
    # The subsequent coefficients correspond to '-c_i'.
    c = [-k for k in coeffs[1:]]

    # The problem asks for the smallest positive 'n'.
    # If n is negative, we can multiply the entire relation by -1.
    if n < 0:
        n = -n
        c = [-val for val in c]
    
    # Unpack the c_i coefficients
    c1, c2, c3, c4, c5, c6 = c

    # Display the final equation
    print("The discovered relation is:")
    equation_parts = [f"{n} * pi/4 ="]
    for i, (val, coeff) in enumerate(zip(args, c)):
        if i > 0:
            equation_parts.append(f" + {coeff}*arctan(1/{val})")
        else:
            equation_parts.append(f" {coeff}*arctan(1/{val})")
    print("".join(equation_parts))
    print("\nWhich is:")
    print(f"{n} * pi/4 = {c1}*arctan(1/{args[0]}) + {c2}*arctan(1/{args[1]}) + {c3}*arctan(1/{args[2]}) + {c4}*arctan(1/{args[3]}) + {c5}*arctan(1/{args[4]}) + {c6}*arctan(1/{args[5]})")

    # Print the coefficients individually
    print("\nThe integer coefficients are:")
    print(f"n  = {n}")
    print(f"c1 = {c1}")
    print(f"c2 = {c2}")
    print(f"c3 = {c3}")
    print(f"c4 = {c4}")
    print(f"c5 = {c5}")
    print(f"c6 = {c6}")

if __name__ == '__main__':
    solve_machin_like_formula()
