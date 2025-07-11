def solve_chromatic_polynomial():
    """
    This function prints the derived chromatic polynomial for the given graph.
    The derivation is based on structural analysis of the graph as explained above.
    The final polynomial is P(G, k) = k^5 - 7k^4 + 19k^3 - 23k^2 + 10k.
    This can also be written in factored form: k(k-1)(k-2)(k^2 - 4k + 5).
    """
    
    # Coefficients of the polynomial from highest power to lowest
    # For k^5 - 7k^4 + 19k^3 - 23k^2 + 10k
    coeffs = {
        5: 1,
        4: -7,
        3: 19,
        2: -23,
        1: 10
    }
    
    poly_string = "P(G, k) = "
    
    is_first_term = True
    for power in sorted(coeffs.keys(), reverse=True):
        coeff = coeffs[power]
        
        # Determine the sign and coefficient part of the term
        sign = ""
        if coeff > 0:
            if not is_first_term:
                sign = "+ "
        elif coeff < 0:
            sign = "- "
        else:
            continue # Skip terms with zero coefficient
        
        # Absolute value of the coefficient
        abs_coeff = abs(coeff)
        
        # Build the term string
        term = sign
        if abs_coeff != 1 or power == 0:
            term += str(abs_coeff)
            
        if power > 0:
            if abs_coeff != 1 and power != 0:
                term += "*" # e.g., 2*k
            term += "k"
            if power > 1:
                term += f"^{power}"
        
        if not is_first_term:
            poly_string += term
        else:
            # For the first term, we remove the leading '+' and space if present
            poly_string += term.lstrip("+ ")

        is_first_term = False

    print("The chromatic polynomial of the graph is:")
    # Print the full expanded equation with each number.
    print(f"P(G, k) = {1}*k^5 - {7}*k^4 + {19}*k^3 - {23}*k^2 + {10}*k")

solve_chromatic_polynomial()