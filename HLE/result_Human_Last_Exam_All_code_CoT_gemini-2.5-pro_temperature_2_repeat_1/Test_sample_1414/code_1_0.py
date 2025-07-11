import mpmath

def find_machin_constants():
    """
    This function finds the integer coefficients n, c1, ..., c6 for the given
    Machin-like formula for pi.
    """
    # Set the precision for mpmath. High precision is necessary for PSLQ to work reliably.
    mpmath.mp.dps = 100

    # The x_k values from the problem statement.
    x_vals = [122, 239, 682, 1252, 2855, 12943]

    # The equation can be written as:
    # c1*arctan(1/122) + ... + c6*arctan(1/12943) - n*(pi/4) = 0
    # We are looking for an integer relation for the vector of numbers:
    # (arctan(1/122), ..., arctan(1/12943), -pi/4)
    
    print("Setting up the vector of values for PSLQ algorithm...")
    vec = [mpmath.atan(mpmath.mpf(1)/x) for x in x_vals]
    vec.append(-mpmath.pi / 4)

    # Use the PSLQ algorithm to find the integer coefficients.
    print("Running PSLQ algorithm...")
    coeffs = mpmath.pslq(vec)
    print("PSLQ algorithm finished.")

    # The coefficients c_k are the first 6 elements of the result.
    c = coeffs[:-1]
    
    # The last coefficient corresponds to -pi/4, so it's -n.
    n = -coeffs[-1]

    # We need the smallest positive n. If n is negative, we flip the signs of all coefficients.
    if n < 0:
        n = -n
        c = [-val for val in c]
    
    print("\nFound the coefficients.")
    
    # Let's display the final equation
    print("The solved equation is:")
    lhs = f"{n} * pi/4"
    rhs_parts = []
    for i in range(len(c)):
        if c[i] == 0:
            continue
        
        # Add sign for the term
        sign = "+ " if c[i] > 0 else "- "
        if not rhs_parts: # First term
            if c[i] > 0:
              sign = ""
            else:
              sign = "-"
        
        # Format the coefficient if it is not 1
        coeff_str = "" if abs(c[i]) == 1 else str(abs(c[i])) + "*"
            
        rhs_parts.append(f"{sign}{coeff_str}arctan(1/{x_vals[i]})")

    print(f"{lhs} = {' '.join(rhs_parts)}")

    # Provide the final answer in the requested format.
    solution = [n] + c
    print("\nThe solution for (n, c1, c2, c3, c4, c5, c6) is:")
    print(tuple(solution))


find_machin_constants()

<<<1, 1, 2, 1, -1, 1, -1>>>