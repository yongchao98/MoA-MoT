def format_jones_polynomial():
    """
    This function formats and prints the Jones polynomial for the 9_42 knot.
    The terms are represented as a list of (coefficient, power) tuples,
    ordered by decreasing power of t.
    """
    # Jones polynomial for the 9_42 knot:
    # t^2 - t + 1 - t^-1 + t^-2 - t^-3 + t^-4
    terms = [
        (1, 2),
        (-1, 1),
        (1, 0),
        (-1, -1),
        (1, -2),
        (-1, -3),
        (1, -4)
    ]

    poly_parts = []
    is_first_term = True

    for coeff, power in terms:
        # Build the string for the current term
        part = ""

        # Determine the sign and absolute coefficient
        sign = " + " if coeff > 0 else " - "
        abs_coeff = abs(coeff)

        # For the very first term, handle the sign separately
        if is_first_term:
            if sign == " - ":
                part += "-"
            is_first_term = False
        else:
            part += sign

        # Format the term based on coefficient and power
        if power == 0:
            part += str(abs_coeff)
        else:
            if abs_coeff != 1:
                part += str(abs_coeff)

            part += "t"

            if power != 1:
                part += f"^{power}"
    
        poly_parts.append(part)

    # Join all the parts into the final polynomial string
    # The first part is added directly, subsequent parts are appended.
    # We need to clean up spaces for the first term.
    final_string = "".join(poly_parts).lstrip(" + ")
    
    print(final_string)

format_jones_polynomial()