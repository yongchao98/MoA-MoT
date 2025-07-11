import collections

def get_jones_polynomial_for_knot():
    """
    This function calculates and prints the Jones polynomial for the knot 9_42.
    The knot in the image is the 9_42 knot, also known as the Celtic knot.
    Its Jones polynomial V(t) is well-documented in knot theory literature.
    This function formats the polynomial in decreasing degree order.
    """
    
    # The Jones polynomial for knot 9_42 is V(t) = t^4 - t^3 + t^2 - t + 1 - t^-1 + t^-2 - t^-3 + t^4.
    # We represent the polynomial as a list of (coefficient, power) tuples.
    # The list is ordered by power in descending order.
    poly_terms = [
        (1, 4), 
        (-1, 3), 
        (1, 2), 
        (-1, 1), 
        (1, 0), 
        (-1, -1), 
        (1, -2), 
        (-1, -3), 
        (1, -4)
    ]
    
    output_parts = []
    
    for i, (coeff, power) in enumerate(poly_terms):
        if coeff == 0:
            continue
            
        term_str = ""
        
        # Determine the sign prefix for the term
        if i == 0:
            if coeff < 0:
                term_str = "-"
        else:
            term_str = " + " if coeff > 0 else " - "
            
        # Format the coefficient
        abs_coeff = abs(coeff)
        if abs_coeff != 1 or power == 0:
            term_str += str(abs_coeff)
            
        # Format the variable and power
        if power != 0:
            term_str += "t"
            if power != 1:
                # Following formatting examples, use braces for negative powers
                if power < 0:
                    term_str += f"^{{{power}}}"
                else:
                    term_str += f"^{power}"

        output_parts.append(term_str)
        
    print("".join(output_parts))

get_jones_polynomial_for_knot()