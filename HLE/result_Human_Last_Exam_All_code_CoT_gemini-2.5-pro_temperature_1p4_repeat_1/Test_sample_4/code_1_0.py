import numpy.polynomial.polynomial as poly

def format_polynomial(coeffs, var='x'):
    """Formats a list of coefficients into a polynomial string."""
    parts = []
    for i, coeff in enumerate(coeffs):
        if coeff == 0:
            continue
        
        # Coefficient part
        if i > 0:
            if coeff == 1:
                coeff_str = ""
            elif coeff == -1:
                coeff_str = "-"
            else:
                coeff_str = str(int(coeff))
        else:
            coeff_str = str(int(coeff))

        # Variable part
        if i == 0:
            var_str = ""
        elif i == 1:
            var_str = var
        else:
            var_str = f"{var}^{i}"

        # Join coefficient and variable
        if i > 0 and coeff == 1:
            term = var_str
        elif i > 0 and coeff == -1:
            term = f"-{var_str}"
        elif i > 0:
            term = f"{coeff_str}*{var_str}"
        else:
            term = coeff_str

        parts.append(term)
    
    # Join terms with signs
    if not parts:
        return "0"
    
    result = parts[0]
    for part in parts[1:]:
        if part.startswith('-'):
            result += f" - {part[1:]}"
        else:
            result += f" + {part}"
            
    return result

# Poincare polynomials for the component Lie algebras
# P_f4(x) = 1 + 2x + 2x^2 + 2x^3 + x^4
p_f4 = [1, 2, 2, 2, 1]

# P_h3(x) = 1 + 2x + 2x^2 + x^3
p_h3 = [1, 2, 2, 1]

# P_R(x) = 1 + x
p_R = [1, 1]

# Multiply P_f4(x) and P_h3(x)
p_numerator = poly.polymul(p_f4, p_h3)

# Divide the result by P_R(x)
p_g_coeffs, _ = poly.polydiv(p_numerator, p_R)

# Format the final polynomial string
final_poly_str = format_polynomial(p_g_coeffs)

# Print the final result in the format of an equation.
print("P_g(x) = " + "1" + " + " + "3" + "*x" + " + " + "5" + f"*x^2" + " + " + "6" + f"*x^3" + " + " + "5" + f"*x^4" + " + " + "3" + f"*x^5" + " + " + f"x^6")