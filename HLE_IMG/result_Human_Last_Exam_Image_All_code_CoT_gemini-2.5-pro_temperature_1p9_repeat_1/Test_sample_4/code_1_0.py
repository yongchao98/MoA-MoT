def get_jones_polynomial_9_42():
    """
    This function calculates and prints the Jones polynomial for the 9_42 knot.
    The polynomial terms (coefficient, power) are hardcoded based on known values.
    The function then constructs and prints the polynomial string in decreasing degree order.
    """
    # Terms of the Jones polynomial for the 9_42 knot
    # Each tuple is in the format (coefficient, power)
    terms = [
        (1, 2),
        (-1, 1),
        (1, 0),
        (-1, -1),
        (1, -2),
        (-1, -3),
        (1, -4)
    ]

    polynomial_parts = []

    for i, (coeff, power) in enumerate(terms):
        # Format the sign
        sign = ""
        if i == 0:
            if coeff < 0:
                sign = "-"
        else:
            sign = " + " if coeff > 0 else " - "
        
        # Format the coefficient
        abs_coeff = abs(coeff)
        coeff_str = ""
        # Only show coefficient if it's not 1, or if it's a constant term (power=0)
        if abs_coeff != 1 or power == 0:
            coeff_str = str(abs_coeff)
            
        # Format the variable 't' and its power
        var_str = ""
        if power != 0:
            var_str = "t"
            if power != 1:
                # Add exponent using '^' for clarity
                var_str += f"^{power}"

        polynomial_parts.append(f"{sign}{coeff_str}{var_str}")
        
    print("".join(polynomial_parts))

if __name__ == '__main__':
    get_jones_polynomial_9_42()
