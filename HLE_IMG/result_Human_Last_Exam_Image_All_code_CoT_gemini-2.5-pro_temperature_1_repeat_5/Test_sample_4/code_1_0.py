import sys

def format_polynomial():
    """
    This function generates and prints the Jones polynomial for the 9_42 knot.
    """
    # The Jones polynomial for knot 9_42 is represented by a list of (coefficient, power) tuples.
    # V(t) = t^4 - t^3 + t^2 - t + 1 - t^-1 + t^-2 - t^-3 + t^-4
    terms = [
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

    poly_str = ""

    for i, (coeff, power) in enumerate(terms):
        # Determine the sign prefix for the term
        sign = ""
        if i > 0:
            if coeff > 0:
                sign = " + "
            else:
                sign = " - "
        elif coeff < 0:
            sign = "-"

        # Get the absolute value of the coefficient
        abs_coeff = abs(coeff)

        # Format the coefficient part of the term
        coeff_str = ""
        if abs_coeff != 1 or power == 0:
            coeff_str = str(abs_coeff)

        # Format the variable part of the term
        var_str = ""
        if power != 0:
            if power == 1:
                var_str = "t"
            else:
                var_str = f"t^{{{power}}}" # Using curly braces for clarity

        # Combine the parts of the term
        if power == 0:
            term_str = coeff_str
        elif coeff_str and var_str:
            term_str = f"{coeff_str}{var_str}"
        else:
            term_str = var_str
            
        poly_str += sign + term_str

    # Print the final formatted polynomial string
    # Replace curly braces with standard characters for the final output
    print(poly_str.replace("{", "").replace("}", ""))

format_polynomial()