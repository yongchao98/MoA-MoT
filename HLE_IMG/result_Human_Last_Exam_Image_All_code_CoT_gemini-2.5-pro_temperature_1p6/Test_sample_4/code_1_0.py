import sys
import io

# Set stdout to use UTF-8 encoding
sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding='utf-8')

def get_jones_polynomial():
    """
    This function prints the Jones polynomial for the knot in the image.
    The knot is the mirror image of the standard 9_42 knot.
    The polynomial is derived by taking the known polynomial for 9_42
    and substituting t with t^-1, then arranging in decreasing degree order.
    """
    
    # The final polynomial's coefficients, with powers as keys.
    # p(t) = t^8 - t^7 + 2t^6 - 2t^5 + 3t^4 - 3t^3 + 3t^2 - 2t + 3 - 2t^-1 + t^-2 - t^-3 + t^-4
    poly_coeffs = {
        8: 1, 7: -1, 6: 2, 5: -2, 4: 3, 3: -3, 2: 3, 1: -2, 0: 3,
        -1: -2, -2: 1, -3: -1, -4: 1
    }

    # Sort powers in descending order to format the polynomial correctly.
    sorted_powers = sorted(poly_coeffs.keys(), reverse=True)

    output_parts = []
    
    # Handle the first term separately to avoid a leading " + ".
    first_power = sorted_powers[0]
    first_coeff = poly_coeffs[first_power]
    
    term_str = ""
    # Add negative sign if coefficient is negative
    if first_coeff < 0:
        term_str += "-"
    
    # Add coefficient value if it's not 1 (or -1), or if it's the constant term.
    if abs(first_coeff) != 1 or first_power == 0:
        term_str += str(abs(first_coeff))

    # Add the variable 't' and its exponent.
    if first_power == 1:
        term_str += "t"
    elif first_power != 0:
        term_str += f"t^{first_power}"
    output_parts.append(term_str)

    # Process the remaining terms.
    for power in sorted_powers[1:]:
        coeff = poly_coeffs[power]
        
        # Add the sign (plus or minus) with spaces.
        term_str = " + " if coeff > 0 else " - "
        
        # Add the coefficient.
        if abs(coeff) != 1 or power == 0:
            term_str += str(abs(coeff))
        
        # Add the variable and exponent.
        if power == 1:
            term_str += "t"
        elif power != 0:
            term_str += f"t^{power}"
            
        output_parts.append(term_str)

    print("".join(output_parts))

get_jones_polynomial()