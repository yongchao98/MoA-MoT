import collections

def get_jones_polynomial_for_9_42():
    """
    This function returns the Jones polynomial for the 9_42 knot.
    The polynomial is pre-determined from knot theory databases.
    The function formats it as a string in decreasing degree order.
    """
    
    # The Jones polynomial for the 9_42 knot is:
    # t^4 - t^3 + 2t^2 - 2t + 3 - 2t^-1 + 2t^-2 - t^-3 + t^-4
    # We represent this as a list of (coefficient, power) tuples.
    # The list is already in decreasing order of power.
    poly_terms = [
        (1, 4),
        (-1, 3),
        (2, 2),
        (-2, 1),
        (3, 0),
        (-2, -1),
        (2, -2),
        (-1, -3),
        (1, -4),
    ]

    poly_parts = []
    is_first_term = True

    for coeff, power in poly_terms:
        # Build the string for the current term
        term_str = ""

        # Handle the sign
        if is_first_term:
            if coeff < 0:
                term_str += "-"
        else:
            if coeff < 0:
                term_str += " - "
            else:
                term_str += " + "
        
        # Handle the coefficient
        abs_coeff = abs(coeff)
        if power == 0:
            term_str += str(abs_coeff)
        elif abs_coeff != 1:
            term_str += str(abs_coeff)

        # Handle the variable 't' and its power
        if power != 0:
            if power == 1:
                term_str += "t"
            else:
                term_str += f"t^{power}"

        poly_parts.append(term_str)
        is_first_term = False

    final_polynomial = "".join(poly_parts)
    print(final_polynomial)

get_jones_polynomial_for_9_42()