import sympy

def get_chromatic_polynomial():
    """
    Computes and formats the chromatic polynomial for the given graph.
    The derivation results in the factored form: k(k-1)(k-2)(k^2 - 4k + 5)
    This function expands it and prints it in a specified format.
    """
    k = sympy.Symbol('k')

    # The chromatic polynomial in factored form, from the deletion-contraction derivation
    poly_factored = k * (k - 1) * (k - 2) * (k**2 - 4*k + 5)

    # Expand the polynomial
    poly_expanded = sympy.expand(poly_factored)

    # Get the polynomial as a dictionary of {power: coefficient}
    poly_dict = poly_expanded.as_poly(k).as_dict()

    # Format the polynomial string to show all numbers
    terms = []
    # Sort terms by power in descending order
    sorted_coeffs = sorted(poly_dict.items(), key=lambda item: item[0][0], reverse=True)

    for i, (power_tuple, coeff) in enumerate(sorted_coeffs):
        power = power_tuple[0]
        if coeff == 0:
            continue

        # Determine the sign for the term
        sign = ""
        if i == 0:
            if coeff < 0:
                sign = "-"
        else:
            if coeff > 0:
                sign = " + "
            else:
                sign = " - "

        # Get the absolute value of the coefficient
        abs_coeff = str(abs(coeff))

        # Format the variable part (k with its power)
        if power == 0:
            var_str = ""
        elif power == 1:
            var_str = "k"
        else:
            var_str = f"k^{power}"
        
        # Combine coefficient and variable, adding a '*' if both exist
        if var_str:
            term = f"{abs_coeff}*{var_str}"
        else:
            term = abs_coeff

        terms.append(f"{sign}{term}")
    
    final_polynomial_str = "".join(terms)
    
    print(f"The chromatic polynomial is P(k) = {final_polynomial_str}")

get_chromatic_polynomial()