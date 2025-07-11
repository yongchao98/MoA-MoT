import sympy

def get_knot_jones_polynomial():
    """
    This function calculates and prints the Jones polynomial for the square knot.
    """
    # Define the symbolic variable for the polynomial
    t = sympy.Symbol('t')

    # The Jones polynomial for the right-handed trefoil knot (3_1)
    # according to the standard convention.
    # V(3_1) = t^-1 + t^-3 - t^-4
    V_3_1 = t**(-1) + t**(-3) - t**(-4)

    # The knot in the image is the square knot (3_1 # -3_1).
    # Its Jones polynomial is the product of the polynomials for the
    # right-handed and left-handed trefoils. The polynomial for the
    # left-handed trefoil is obtained by substituting t with 1/t.
    V_neg_3_1 = V_3_1.subs(t, 1/t)

    # Calculate the product of the two polynomials
    V_square_knot = sympy.expand(V_3_1 * V_neg_3_1)

    # Extract the terms of the resulting polynomial to sort and format them.
    # A Poly object is useful for handling the terms.
    # The domain ZZ[t, 1/t] is for Laurent polynomials with integer coefficients.
    poly = sympy.Poly(V_square_knot, t, 1/t)
    
    # Get terms as a list of [((power_of_t, power_of_1/t), coefficient), ...]
    terms = poly.terms()

    # Sort terms by degree in descending order. The power is p_t - p_1/t.
    sorted_terms = sorted(terms, key=lambda item: item[0][0] - item[0][1], reverse=True)

    # Build the final string representation
    result_string = ""
    for i, term in enumerate(sorted_terms):
        # The power is the difference between the power of t and the power of 1/t
        power = term[0][0] - term[0][1]
        coeff = term[1]

        # Determine the sign to prepend (+, -, or nothing for the first term)
        if i == 0:
            sign = "-" if coeff < 0 else ""
        else:
            sign = " - " if coeff < 0 else " + "
        
        abs_coeff = abs(coeff)

        # Format the coefficient part of the term
        # Hide coefficient if it's 1, unless it's a constant term (power=0)
        if abs_coeff == 1 and power != 0:
            coeff_str = ""
        else:
            coeff_str = str(abs_coeff)

        # Format the variable part of the term (t with its power)
        if power == 0:
            power_str = "" # Constant term
        elif power == 1:
            power_str = "t"
        else:
            power_str = f"t^{power}"
            
        result_string += f"{sign}{coeff_str}{power_str}"

    print(result_string)

get_knot_jones_polynomial()