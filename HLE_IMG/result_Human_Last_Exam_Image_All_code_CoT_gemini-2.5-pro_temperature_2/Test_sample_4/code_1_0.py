import sys

# The knot in the image is the mirror of the standard 9_42 knot.
# Its Jones polynomial V(t) is derived from the polynomial of the standard 9_42 knot
# by replacing t with t^-1.
# The standard polynomial is V_pos(t) = t^8 - t^7 + t^6 - t^5 + t^4 - t^3 + t^2 - t + 1.
# So, for the given knot, V_neg(t) = (t^-1)^8 - (t^-1)^7 + ... - t^-1 + 1
# V_neg(t) = t^-8 - t^-7 + t^-6 - t^-5 + t^-4 - t^-3 + t^-2 - t^-1 + 1.
#
# Arranging in decreasing degree order:
# 1 - t^-1 + t^-2 - t^-3 + t^-4 - t^-5 + t^-6 - t^-7 + t^-8
#
# The coefficients for powers 0, -1, -2, ..., -8 are [1, -1, 1, -1, 1, -1, 1, -1, 1].

def get_jones_polynomial_string():
    """
    This function generates the string for the Jones polynomial of the knot 9_42 (mirror).
    """
    coeffs = [1, -1, 1, -1, 1, -1, 1, -1, 1]
    powers = [0, -1, -2, -3, -4, -5, -6, -7, -8]

    # Use a file-like object to build the string to handle printing logic easily
    from io import StringIO
    f = StringIO()

    for i in range(len(coeffs)):
        coeff = coeffs[i]
        power = powers[i]

        # Print sign for terms after the first one
        if i > 0:
            if coeff > 0:
                f.write(" + ")
            else:
                f.write(" - ")
        # Handle the very first term's sign
        elif coeff < 0:
            f.write("-")

        # Print the absolute value of the coefficient
        f.write(str(abs(coeff)))

        # Print the variable and power
        # No variable for constant term
        if power == 0:
            continue
        
        # Standard format for other powers
        if sys.version_info >= (3, 8):
             # This uses f-string debugging specifiers which is a clean way, but let's be more compatible
             f.write(f"t^({power})")
        else: # more compatible way
             f.write("t^(" + str(power) + ")")


    return f.getvalue()

if __name__ == '__main__':
    # Print the full polynomial equation string
    # We simplified the expression according to the example formatting, e.g., not printing "1*"
    # but the prompt insisted on outputting *each number*. 
    # To resolve this, let's go with a format that is explicit but still clean.
    # Format example: -3t^3 + t - t^{-1} + 5t^{-10}
    # This example omits coefficients of 1 and -1.
    
    coeffs = [1, -1, 1, -1, 1, -1, 1, -1, 1]
    powers = [0, -1, -2, -3, -4, -5, -6, -7, -8]

    final_poly = []

    for i, (c, p) in enumerate(zip(coeffs, powers)):
        # Sign
        if i == 0:
            sign = "-" if c < 0 else ""
        else:
            sign = " - " if c < 0 else " + "
        
        # Coefficient string
        coeff_str = ""
        if abs(c) != 1 or p == 0:
             coeff_str = str(abs(c))
        
        # Variable string
        if p == 0:
            var_str = ""
        elif p == 1:
            var_str = "t"
        else:
            var_str = f"t^({p})"

        # Assemble the term
        term = sign + coeff_str + var_str
        final_poly.append(term)
    
    print("".join(final_poly))