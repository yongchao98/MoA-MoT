import mpmath

def solve_arctan_relation():
    """
    Finds the integer coefficients for the given Machin-like formula
    using the PSLQ algorithm.
    """
    # Set the precision for the calculations. PSLQ requires high precision.
    mpmath.mp.dps = 100

    # These are the denominators from the arctan terms in the equation.
    denominators = [122, 239, 682, 1252, 2855, 12943]

    # Create the vector of real numbers for which to find an integer relation.
    # The first term is pi/4, followed by the arctan values.
    terms = [mpmath.pi / 4] + [mpmath.atan(mpmath.mpf(1)/d) for d in denominators]

    # Use the PSLQ algorithm to find the integer relation vector.
    # The relation is of the form: a[0]*x[0] + a[1]*x[1] + ... = 0
    coeffs = mpmath.pslq(terms)

    # From the problem, the relation is n*(pi/4) - c1*arctan(...) - ... = 0
    # So, n = coeffs[0] and c_i = -coeffs[i+1].
    n = coeffs[0]
    c = [-k for k in coeffs[1:]]

    # The problem asks for the smallest positive n.
    # If PSLQ returns a relation with a negative n, we just flip all the signs.
    if n < 0:
        n = -n
        c = [-val for val in c]

    # Now, we format and print the final equation.
    print("The discovered equation is:")
    
    equation_parts = []
    for i, coeff in enumerate(c):
        if coeff == 0:
            continue
        
        # Determine the sign for the term
        sign = " + " if coeff > 0 else " - "
        
        # Format the coefficient (don't show 1)
        abs_coeff = abs(coeff)
        if abs_coeff == 1:
            coeff_str = ""
        else:
            coeff_str = f"{abs_coeff}*"
        
        term = f"{coeff_str}arctan(1/{denominators[i]})"
        equation_parts.append((sign, term))
    
    # Join all parts into the right-hand side of the equation
    rhs = "".join([part[1] if i == 0 and equation_parts[i][0] == " + " else part[0] + part[1] for i, part in enumerate(equation_parts)])
    
    # Correct the starting sign if it's "+ "
    if rhs.startswith(" + "):
        rhs = rhs[3:]

    print(f"{n} * pi/4 = {rhs}")

    # Return the final answer string for the submission format
    return f"{n},{','.join(map(str, c))}"

if __name__ == '__main__':
    # Run the solver and print the result.
    # The problem asks for the final answer in a specific format,
    # but printing the coefficients first can be helpful.
    final_coeffs = solve_arctan_relation().split(',')
    n, c = final_coeffs[0], final_coeffs[1:]
    
    print("\nThe coefficients are:")
    print(f"n = {n}")
    for i in range(len(c)):
        print(f"c_{i+1} = {c[i]}")
