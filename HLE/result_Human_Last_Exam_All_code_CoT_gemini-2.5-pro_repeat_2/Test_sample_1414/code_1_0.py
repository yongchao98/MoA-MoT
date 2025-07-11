import mpmath

def solve_machin_like_formula():
    """
    Solves the given Machin-like formula for PI using the PSLQ algorithm.
    The formula is:
    n*pi/4 = c1*atan(1/122) + c2*atan(1/239) + c3*atan(1/682) + 
             c4*atan(1/1252) + c5*atan(1/2855) + c6*atan(1/12943)
    """
    
    # The mpmath library is required. You can install it with:
    # pip install mpmath
    
    # Set the precision for the calculations.
    mpmath.mp.dps = 100
    
    # Define the denominators of the arctan terms
    denominators = [122, 239, 682, 1252, 2855, 12943]
    
    # Create the vector of high-precision numbers for PSLQ.
    # The equation is n*(pi/4) - c1*atan(1/x1) - ... = 0
    # The input vector for pslq will be [pi/4, atan(1/x1), atan(1/x2), ...]
    pi_4 = mpmath.pi / 4
    arctan_terms = [mpmath.atan(mpmath.mpf(1) / d) for d in denominators]
    input_vector = [pi_4] + arctan_terms
    
    # Find the integer relation using PSLQ.
    # pslq finds a vector r such that r[0]*v[0] + r[1]*v[1] + ... = 0
    relation = mpmath.pslq(input_vector)
    
    # Extract n and c_k from the relation vector.
    # relation[0] corresponds to n.
    # relation[k+1] corresponds to -c_k.
    n = relation[0]
    c = [-k for k in relation[1:]]
    
    # The problem asks for the smallest positive n.
    # If n from PSLQ is negative, we negate the whole solution vector.
    if n < 0:
        n = -n
        c = [-k for k in c]

    # Print the full equation with the found coefficients.
    print("The final equation is:")
    
    equation_str = f"{n} * pi/4 = "
    for i in range(len(c)):
        term_val = c[i]
        den_val = denominators[i]
        
        # Add sign and coefficient
        if i > 0:
            if term_val >= 0:
                equation_str += "+ "
            else:
                equation_str += "- "
        elif term_val < 0:
             equation_str += "- "
        
        # Add coefficient if not 1
        abs_coeff = abs(term_val)
        if abs_coeff != 1:
            equation_str += f"{abs_coeff}*"

        equation_str += f"arctan(1/{den_val}) "

    print(equation_str)
    print("\n" + "="*20 + "\n")

    # Print the final answer in the requested format
    print("The constants [n, c1, c2, c3, c4, c5, c6] are:")
    result = [n] + c
    print(','.join(map(str, result)))

solve_machin_like_formula()
