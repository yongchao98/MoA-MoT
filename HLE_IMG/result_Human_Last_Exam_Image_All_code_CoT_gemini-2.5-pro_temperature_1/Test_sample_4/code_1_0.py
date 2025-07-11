import numpy as np

def get_jones_polynomial_10_132():
    """
    Calculates and prints the Jones polynomial for the knot 10_132.

    The knot in the image is the 10_132 knot.
    Its Jones polynomial V(t) can be calculated from the formula:
    V(t) = t^(-7) * (-1 + 2t - 2t^2 + 3t^3 - 2t^4 + 2t^5 - t^6 + t^7)^2
    This script expands this formula and prints the result in the required format.
    """
    
    # Coefficients of the polynomial P(t) = -1 + 2t - 2t^2 + ... + t^7
    # The i-th element is the coefficient of t^i
    p_coeffs = [-1, 2, -2, 3, -2, 2, -1, 1]

    # Square P(t) by convolving its coefficient list with itself
    p_squared_coeffs = np.convolve(p_coeffs, p_coeffs)

    # Store the terms of the final polynomial as a list of dictionaries
    final_terms = []
    # The leading t^(-7) shifts the power of each term t^i to t^(i-7)
    for i, coeff in enumerate(p_squared_coeffs):
        if coeff != 0:
            power = i - 7
            final_terms.append({'coeff': int(coeff), 'power': power})

    # Sort terms in descending order of power
    final_terms.sort(key=lambda x: x['power'], reverse=True)

    # Helper function to format a single term (e.g., -3t^3)
    def format_term(coeff, power, is_first_term):
        # Determine the sign to prepend
        if is_first_term:
            sign = "-" if coeff < 0 else ""
        else:
            sign = " - " if coeff < 0 else " + "
        
        abs_coeff = abs(coeff)
        
        # Format the coefficient string (hide 1 for non-constant terms)
        if abs_coeff == 1 and power != 0:
            coeff_str = ""
        else:
            coeff_str = str(abs_coeff)
            
        # Format the variable part with its power
        if power == 0:
            var_str = ""
        elif power == 1:
            var_str = "t"
        else:
            var_str = f"t^{power}"
            
        return f"{sign}{coeff_str}{var_str}"

    # Build the final string representation of the polynomial
    output_parts = []
    for i, term in enumerate(final_terms):
        output_parts.append(format_term(term['coeff'], term['power'], i == 0))
    
    # Print the complete polynomial
    print("".join(output_parts))

get_jones_polynomial_10_132()