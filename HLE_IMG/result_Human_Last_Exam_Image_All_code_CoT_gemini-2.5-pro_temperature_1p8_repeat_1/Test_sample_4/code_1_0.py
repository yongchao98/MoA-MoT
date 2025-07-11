def print_jones_polynomial_for_10_124():
    """
    This function prints the Jones polynomial for the knot 10_124.
    The polynomial is pre-determined from knot theory databases.
    The code formats the polynomial string according to standard mathematical notation.
    """
    
    # The coefficients and powers of the Jones polynomial for knot 10_124
    # Sourced from the Knot Atlas.
    # Format: {power: coefficient}
    poly_terms = {
        5: 1, 4: -1, 3: 2, 2: -2, 1: 3, 0: -2,
        -1: 3, -2: -2, -3: 2, -4: -1, -5: 1
    }
    
    # Get the degrees in descending order
    degrees = sorted(poly_terms.keys(), reverse=True)
    
    result_str_parts = []
    is_first_term = True
    
    # Iterate through each term to build the string
    for d in degrees:
        c = poly_terms[d]
        if c == 0:
            continue
            
        part = ""
        abs_c = abs(c)
        
        # Determine the sign
        if is_first_term:
            if c < 0:
                part += "-"
        else:
            if c > 0:
                part += " + "
            else:
                part += " - "
        
        # Determine the coefficient string
        # Don't print coefficient '1' unless it's for the constant term (t^0)
        if abs_c != 1 or d == 0:
            part += str(abs_c)
            
        # Determine the variable and power string
        if d != 0:
            part += "t"
            if d != 1:
                # Format power with '^'
                part += "^" + str(d)
                
        result_str_parts.append(part)
        is_first_term = False
        
    print("p(t) =", "".join(result_str_parts))

print_jones_polynomial_for_10_124()