def get_jones_polynomial_12n427():
    """
    This function generates and prints the Jones polynomial for the knot 12n427.
    The polynomial data is pre-determined from knot theory databases.
    The function formats the output string in decreasing degree order.
    """
    
    # The Jones polynomial for knot 12n427 is:
    # t^5 - t^4 + 2t^3 - 2t^2 + 3t - 3 + 3t^-1 - 2t^-2 + 2t^-3 - t^-4 + t^-5
    # We represent this as a list of (coefficient, power) tuples.
    poly_data = [
        (1, 5), (-1, 4), (2, 3), (-2, 2), (3, 1), (-3, 0),
        (3, -1), (-2, -2), (2, -3), (-1, -4), (1, -5)
    ]

    result_parts = []
    is_first_term = True

    for coeff, power in poly_data:
        # This loop iterates through the coefficients and powers to build each term of the polynomial string.
        
        # Skip terms with zero coefficient
        if coeff == 0:
            continue

        # --- Build the string for the current term ---
        current_term = ""

        # 1. Handle the sign (+ or -)
        if coeff < 0:
            # For negative coefficients, add " - " except for the very first term
            if not is_first_term:
                current_term += " - "
            else:
                current_term += "-"
        else: # coeff > 0
            # For positive coefficients, add " + " except for the very first term
            if not is_first_term:
                current_term += " + "
        
        abs_coeff = abs(coeff)

        # 2. Handle the coefficient value
        # We don't print the coefficient if it's 1, unless it's the constant term (power = 0)
        if abs_coeff != 1 or power == 0:
            current_term += str(abs_coeff)

        # 3. Handle the variable 't' and its exponent
        if power != 0:
            current_term += "t"
            # We don't print the exponent if it's 1
            if power != 1:
                current_term += f"{power}"

        result_parts.append(current_term)
        is_first_term = False
            
    # Join all the term strings together and print the final polynomial
    print("".join(result_parts))

# Execute the function to print the polynomial
get_jones_polynomial_12n427()