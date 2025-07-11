def get_jones_polynomial_9_42():
    """
    This function returns the Jones polynomial for the 9_42 knot.
    The polynomial is represented as a list of tuples,
    where each tuple is (coefficient, exponent).
    The list is ordered by decreasing degree of the exponent.
    """
    # Jones polynomial for knot 9_42 from the Knot Atlas
    # V(t) = t^-1 - t^-2 + 2t^-3 - 2t^-4 + 2t^-5 - 2t^-6 + 2t^-7 - t^-8 + t^-9 - t^-10 + t^-11
    terms = [
        (1, -1), (-1, -2), (2, -3), (-2, -4), (2, -5), (-2, -6),
        (2, -7), (-1, -8), (1, -9), (-1, -10), (1, -11)
    ]
    return terms

def format_polynomial(terms):
    """
    Formats a list of (coefficient, exponent) tuples into a human-readable string.
    """
    poly_string = ""
    for i, (coeff, exp) in enumerate(terms):
        # Format the sign
        sign = ""
        if i == 0:
            if coeff < 0:
                sign = "-"
        else:
            if coeff > 0:
                sign = " + "
            else:
                sign = " - "

        # Format the coefficient
        abs_coeff = abs(coeff)
        if abs_coeff == 1 and exp != 0:
            coeff_str = ""
        else:
            coeff_str = str(abs_coeff)

        # Format the variable and exponent part
        if exp == 0:
            var_str = ""
            # If the variable part is gone, make sure we print the coefficient if it's 1
            if coeff_str == "":
                coeff_str = "1"
        elif exp == 1:
            var_str = "t"
        else:
            # Use curly braces for exponents as per the requested format
            var_str = f"t^{{{exp}}}"
            
        poly_string += f"{sign}{coeff_str}{var_str}"

    return poly_string

# Get the terms for the Jones polynomial of knot 9_42
polynomial_terms = get_jones_polynomial_9_42()

# Format the polynomial into a string
final_equation = format_polynomial(polynomial_terms)

# Print the final equation
print(final_equation)
