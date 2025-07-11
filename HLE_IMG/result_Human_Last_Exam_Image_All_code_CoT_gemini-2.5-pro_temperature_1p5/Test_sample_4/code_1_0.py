def print_jones_polynomial():
    """
    This function prints the Jones polynomial for the knot 9_42.
    The polynomial is t^8 - t^7 + t^6 - t^5 + t^4 - t^3 + t^2.
    This is formatted in decreasing degree order as requested.
    """

    # The polynomial is represented by a list of (coefficient, degree) tuples.
    poly_terms = [
        (1, 8),
        (-1, 7),
        (1, 6),
        (-1, 5),
        (1, 4),
        (-1, 3),
        (1, 2),
    ]

    output_str = ""
    for i, (coeff, degree) in enumerate(poly_terms):
        # Format the sign
        if i == 0:
            if coeff == -1:
                output_str += "-"
        else:
            if coeff == 1:
                output_str += " + "
            elif coeff == -1:
                output_str += " - "
            else: # general case, not needed here
                 if coeff > 0:
                     output_str += f" + {coeff}"
                 else:
                     output_str += f" - {abs(coeff)}"

        # Format the coefficient (if not 1 or -1) and the variable with its exponent
        # Since all absolute coefficients are 1, we don't need to print them.
        output_str += f"t^{degree}"

    print(output_str)

print_jones_polynomial()