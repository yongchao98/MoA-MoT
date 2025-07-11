import sympy

def generate_and_print_polynomial(n):
    """
    Calculates the polynomial f_n(p) based on the derived formula
    and prints it, showing all numeric coefficients explicitly.
    The formula used is f_n(p) = (p^n - (1-p)^n) / (2p-1).
    """
    if not isinstance(n, int) or n < 1:
        print("Error: n must be a positive integer.")
        return

    p = sympy.symbols('p')
    
    # The closed-form formula for f_n(p)
    formula = (p**n - (1 - p)**n) / (2*p - 1)
    
    # Use sympy.cancel to simplify the rational expression into a polynomial
    polynomial_expr = sympy.cancel(formula)

    # Convert the simplified expression to a Poly object to easily access coefficients
    poly = sympy.Poly(polynomial_expr, p)
    
    # Get the coefficients and the degree of the polynomial
    coeffs = [int(c) for c in poly.all_coeffs()]
    degree = poly.degree()
    
    parts = []
    for i, coeff in enumerate(coeffs):
        if coeff == 0:
            continue
            
        power = degree - i
        
        # Determine the sign prefix for the term
        if not parts:  # For the very first term
            if coeff < 0:
                parts.append("-")
        else:  # For subsequent terms
            parts.append("+" if coeff > 0 else "-")
        
        # Add a space for operators after the first term
        if len(parts) > 1 and parts[-1] in ["+", "-"]:
             parts.insert(-1, " ")
             parts.append(" ")
        
        # Append the absolute value of the coefficient
        parts.append(str(abs(coeff)))
        
        # Append the variable part if the power is greater than 0
        if power > 0:
            parts.append(" * ")
            if power == 1:
                parts.append("p")
            else:
                parts.append(f"p^{power}")

    # For the case f_n(p)=1, where coeffs=[1] and the loop results in ['1']
    if not parts:
        if coeffs and coeffs[0] == 1:
            parts.append("1")

    # Print the final assembled polynomial string
    print(f"f_{n}(p) = {''.join(parts)}")

# Example: Generate and print the polynomial for n=19
generate_and_print_polynomial(19)