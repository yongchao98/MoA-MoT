import numpy as np

def poly_to_string(p, var='x'):
    """Converts a list of polynomial coefficients into a string representation."""
    res = []
    for i, coeff in enumerate(p):
        if coeff == 0:
            continue
        
        # Coefficient part
        if i > 0:
            sign = " + " if coeff > 0 else " - "
            coeff_val = abs(coeff)
        else: # i == 0
            sign = "" if coeff > 0 else "-"
            coeff_val = abs(coeff)

        if coeff_val == 1 and i != 0:
            coeff_str = ""
        else:
            coeff_str = str(coeff_val)

        # Variable part
        if i == 0:
            var_str = ""
        elif i == 1:
            var_str = f"{var}"
        else:
            var_str = f"{var}^{i}"
            
        # Asterisk for multiplication
        if coeff_str and var_str:
            asterisk = "*"
        else:
            asterisk = ""

        term = f"{sign}{coeff_str}{asterisk}{var_str}"
        res.append(term)
        
    # Join terms, handling the leading '+' for positive polynomials
    if res[0].startswith(" + "):
        res[0] = res[0][3:]
    return "".join(res)

def print_poly_equation(p, var='x'):
    """Prints the polynomial in a format showing each number."""
    terms = []
    for i, coeff in enumerate(p):
        if coeff == 0 and len(p) > 1:
            continue
        
        # variable part
        if i == 0:
            var_part = ""
        elif i == 1:
            var_part = f" * {var}"
        else:
            var_part = f" * {var}^{i}"

        terms.append(f"{coeff}{var_part}")
    
    print(" + ".join(terms))


def main():
    """
    Computes the Poincaré polynomial for the given Lie algebra.
    """
    # P_A(x) = 1 + 2x + 2x^2 + 2x^3 + x^4
    p_A = [1, 2, 2, 2, 1]
    
    # P_B(x) = 1 + 2x + 2x^2 + x^3
    p_B = [1, 2, 2, 1]
    
    # P_C(x) = 1 + x
    p_C = [1, 1]

    # Calculate P_g(x) = P_A(x) * P_B(x) / P_C(x)
    # This is equivalent to (P_A(x) / P_C(x)) * P_B(x)
    # P_A(x) / P_C(x) = (1 + 2x + 2x^2 + 2x^3 + x^4) / (1 + x) = 1 + x + x^2 + x^3
    p_A_div_C = np.polydiv(p_A[::-1], p_C[::-1])[0][::-1]
    
    # Now multiply (P_A(x)/P_C(x)) with P_B(x)
    p_g = np.convolve(p_A_div_C, p_B).astype(int).tolist()

    print("The Poincaré polynomial is P(x) = ", end="")
    print_poly_equation(p_g, var='x')

if __name__ == "__main__":
    main()
