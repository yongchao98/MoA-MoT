import mpmath

def find_pi_formula_coeffs():
    """
    Solves for the integer coefficients n, c1, ..., c6 in the given
    Machin-like formula for pi using the PSLQ algorithm.
    """
    # Set the precision for the calculation. PSLQ needs high precision.
    mpmath.mp.dps = 200

    # The denominators of the arctan arguments
    denominators = [122, 239, 682, 1252, 2855, 12943]

    # The vector of real numbers for which we seek an integer relation.
    # The relation is of the form:
    # a_0*x_0 + a_1*x_1 + ... = 0
    # Our equation maps to this with x_0=pi/4, x_1=arctan(1/122), etc.
    # and a_0=n, a_1=-c1, a_2=-c2, etc.
    vals = [mpmath.pi / 4] + [mpmath.atan(mpmath.mpf(1) / d) for d in denominators]

    # Use the PSLQ algorithm to find the integer coefficients.
    coeffs = mpmath.pslq(vals)

    # Ensure n (the first coefficient) is positive, as requested.
    # If it's negative, we can multiply the entire relation by -1.
    if coeffs and coeffs[0] < 0:
        coeffs = [-c for c in coeffs]

    # Extract the coefficients n, c1, ..., c6.
    # coeffs corresponds to [n, -c1, -c2, -c3, -c4, -c5, -c6]
    n = coeffs[0]
    c_coeffs = [-val for val in coeffs[1:]]

    # Print the discovered equation with the found coefficients.
    print("The discovered relation is:")
    
    equation_str = f"{n} * pi/4 = "
    first_term = True
    for i in range(len(c_coeffs)):
        coeff = c_coeffs[i]
        denom = denominators[i]
        
        # This part handles pretty printing of the formula
        if coeff == 0:
            continue
            
        sign = " + " if coeff > 0 else " - "
        coeff_abs = abs(coeff)
        
        if first_term:
            sign = "" if coeff > 0 else "-"
            
        term_str = ""
        if coeff_abs == 1:
            term_str = f"arctan(1/{denom})"
        else:
            term_str = f"{coeff_abs}*arctan(1/{denom})"
        
        equation_str += f"{sign}{term_str}"
        first_term = False
        
    print(equation_str)
    print("\n" + "="*50 + "\n")

    # Print the final solution in the requested format.
    print("The final coefficients are:")
    result_list = [n] + c_coeffs
    print(','.join(map(str, result_list)))

# Run the function to get the answer.
find_pi_formula_coeffs()