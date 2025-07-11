def print_chromatic_polynomial():
    """
    This function prints the computed chromatic polynomial for the given graph.
    The polynomial is P(k) = k^5 - 7k^4 + 19k^3 - 23k^2 + 10k.
    """
    
    # Coefficients for the powers k^5, k^4, k^3, k^2, k^1
    coeffs = {
        5: 1,
        4: -7,
        3: 19,
        2: -23,
        1: 10
    }
    
    # Start building the output string
    poly_string = "P(k) = "
    
    is_first_term = True
    for power in sorted(coeffs.keys(), reverse=True):
        coeff = coeffs[power]
        
        # Determine the sign to print
        sign = ""
        if not is_first_term:
            sign = " + " if coeff > 0 else " - "
        elif coeff < 0:
            sign = "-"

        abs_coeff = abs(coeff)
        
        # Format the term
        term = ""
        # Coefficient part
        if abs_coeff != 1 or power == 0:
            term += str(abs_coeff)
        
        # Variable part
        if power > 0:
            term += "k"
            if power > 1:
                term += f"^{power}"
        
        poly_string += sign + term
        is_first_term = False
        
    print("The computed chromatic polynomial is:")
    print(poly_string)

print_chromatic_polynomial()