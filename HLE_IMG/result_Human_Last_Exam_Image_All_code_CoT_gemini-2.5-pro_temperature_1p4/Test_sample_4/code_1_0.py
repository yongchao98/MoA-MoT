def get_jones_polynomial():
    """
    This function calculates and prints the Jones polynomial for the given knot (10_139 mirror).
    The polynomial is constructed from its coefficients and powers and formatted
    in decreasing degree order.
    """
    # The terms of the Jones polynomial for the knot 10_139 (mirror)
    # stored as a list of (coefficient, power) tuples in decreasing power order.
    polynomial_terms = [
        (-1, 2),
        (1, 1),
        (-1, 0),
        (2, -1),
        (-2, -2),
        (3, -3),
        (-2, -4),
        (2, -5),
        (-1, -6),
        (1, -7)
    ]

    output_parts = []
    
    # Format the first term
    coeff, power = polynomial_terms[0]
    if coeff == 1 and power != 0:
        term_str = ""
    elif coeff == -1 and power != 0:
        term_str = "-"
    else:
        term_str = str(coeff)

    if power == 1:
        term_str += "t"
    elif power != 0:
        term_str += f"t^{{{power}}}"
    
    output_parts.append(term_str)

    # Format the remaining terms
    for coeff, power in polynomial_terms[1:]:
        # Determine the sign part
        if coeff > 0:
            sign = " + "
        else:
            sign = " - "
        
        abs_coeff = abs(coeff)
        
        # Determine the coefficient part
        if abs_coeff == 1 and power != 0:
            coeff_str = ""
        else:
            coeff_str = str(abs_coeff)
        
        # Determine the variable part
        if power == 0:
            var_str = ""
        elif power == 1:
            var_str = "t"
        else:
            var_str = f"t^{{{power}}}"
            
        output_parts.append(f"{sign}{coeff_str}{var_str}")

    print("".join(output_parts))

get_jones_polynomial()