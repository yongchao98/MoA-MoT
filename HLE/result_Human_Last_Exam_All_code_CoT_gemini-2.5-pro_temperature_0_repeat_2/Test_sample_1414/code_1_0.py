import mpmath

def solve_arctan_relation():
    """
    Solves the given arctan relation problem using the PSLQ algorithm.
    """
    # Set the precision for mpmath. High precision is needed for PSLQ to work reliably.
    mpmath.mp.dps = 100

    # The denominators of the arctan arguments
    denominators = [122, 239, 682, 1252, 2855, 12943]

    # Create the vector of high-precision numbers for PSLQ
    # The terms are arctan(1/x_i) and pi/4
    terms = [mpmath.atan(mpmath.mpf(1) / d) for d in denominators]
    terms.append(mpmath.pi / 4)

    # Use the PSLQ algorithm to find an integer relation among the terms.
    # It finds coefficients 'a' such that a[0]*terms[0] + a[1]*terms[1] + ... = 0
    try:
        coeffs = mpmath.pslq(terms)
    except Exception as e:
        print(f"An error occurred with PSLQ: {e}")
        print("This may be due to insufficient precision or no simple relation found.")
        return

    # The coefficients should be integers. We can normalize them to get the smallest integers.
    # Find the smallest non-zero coefficient in absolute value to normalize.
    min_abs_coeff = min(abs(c) for c in coeffs if c != 0)
    coeffs_int = [int(round(c / min_abs_coeff)) for c in coeffs]

    # The last coefficient corresponds to the pi/4 term.
    # The relation is c1*t1 + ... + c6*t6 + c_n*(pi/4) = 0
    # The problem is n*pi/4 = c'1*t1 + ...
    # This means c_n = -n.
    cs = coeffs_int[:-1]
    n = -coeffs_int[-1]

    # The problem asks for the smallest positive n.
    # If n is negative, we can multiply all coefficients by -1.
    if n < 0:
        n = -n
        cs = [-c for c in cs]

    # Print the final equation
    equation_parts = []
    for i in range(len(cs)):
        c = cs[i]
        d = denominators[i]
        if c != 0:
            sign = "+" if c > 0 else "-"
            abs_c = abs(c)
            # Handle coefficient of 1
            coeff_str = f"{abs_c} * " if abs_c != 1 else ""
            
            if not equation_parts: # First term
                if sign == "+":
                    sign = "" # No leading plus sign
                equation_parts.append(f"{sign}{coeff_str}arctan(1/{d})")
            else:
                 equation_parts.append(f" {sign} {coeff_str}arctan(1/{d})")

    print(f"{n} * pi / 4 = {''.join(equation_parts)}")
    
    # Print the constants as requested
    print(f"The constants are:")
    print(f"n = {n}")
    for i in range(len(cs)):
        print(f"c_{i+1} = {cs[i]}")
        
    # Final answer in the required format
    final_answer_str = f"{n}, " + ", ".join(map(str, cs))
    # This is for the platform to read, not for the user to copy.
    # print(f"\n<<< {final_answer_str} >>>")


solve_arctan_relation()