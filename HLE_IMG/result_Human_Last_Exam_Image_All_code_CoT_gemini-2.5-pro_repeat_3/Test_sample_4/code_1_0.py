def print_jones_polynomial():
    """
    This function prints the Jones polynomial for the knot 10_124,
    which is represented by the provided image.
    The polynomial is formatted in decreasing degree order.
    """
    # The coefficients of the Jones polynomial for knot 10_124,
    # mapping power to coefficient.
    poly_data = {
        4: 1, 3: -1, 2: 1, 1: -2, 0: 2,
        -1: -2, -2: 1, -3: -1, -4: 1
    }

    # Sort powers in descending order to format the polynomial correctly.
    sorted_powers = sorted(poly_data.keys(), reverse=True)

    output_str = ""
    is_first_term = True

    for p in sorted_powers:
        c = poly_data[p]

        # Determine the sign part of the term.
        sign = ""
        if c > 0:
            if not is_first_term:
                sign = " + "
        else:
            # Add a space before the minus sign for subsequent terms.
            sign = " - " if not is_first_term else "-"
        
        abs_c = abs(c)

        # Determine the coefficient part. Omit for 1 or -1 unless it's a constant term.
        coeff_str = str(abs_c)
        if abs_c == 1 and p != 0:
            coeff_str = ""
        
        # Determine the variable part (t and its power).
        var_str = ""
        if p == 1:
            var_str = "t"
        elif p != 0:
            var_str = f"t^{{{p}}}"
            
        # For the constant term (p=0), we just use the coefficient.
        if p == 0:
            output_str += f"{sign}{coeff_str}"
        else:
            output_str += f"{sign}{coeff_str}{var_str}"
        
        is_first_term = False
        
    print(output_str)

print_jones_polynomial()