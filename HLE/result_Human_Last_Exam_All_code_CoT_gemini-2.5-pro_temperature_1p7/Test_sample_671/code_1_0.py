import sympy

def find_polynomial_expression(n):
    """
    Computes the polynomial f_n(p) based on the derived formula:
    f_n(p) = (p^n - (1-p)^n) / (2p - 1)

    Args:
        n (int): The index of the polynomial in the sequence.

    Returns:
        sympy.Expr: The symbolic expression for the polynomial f_n(p).
    """
    p = sympy.Symbol('p')
    
    # The formula has a removable singularity at p=1/2.
    # Sympy's `cancel` function correctly simplifies the fraction, 
    # effectively performing the polynomial division and handling the limit.
    formula = (p**n - (1 - p)**n) / (2 * p - 1)
    polynomial = sympy.cancel(formula)
    
    return polynomial

def format_polynomial_for_output(poly_expr, n):
    """
    Formats a sympy polynomial expression into a specific string format
    that resembles the problem description, printing all numeric coefficients.
    
    Args:
        poly_expr (sympy.Expr): The polynomial to format.
        n (int): The index of the polynomial for the label f_n(p).

    Returns:
        str: The formatted string representation of the polynomial.
    """
    p = sympy.Symbol('p')
    poly = sympy.Poly(poly_expr, p)
    
    terms = []
    
    # Get a dictionary of {power: coefficient}
    coeffs_dict = poly.as_dict()
    # Sort terms by power in descending order
    sorted_powers = sorted(coeffs_dict.keys(), reverse=True)
    
    for i, power_tuple in enumerate(sorted_powers):
        power = power_tuple[0]
        coeff = int(coeffs_dict[power_tuple])
        
        if coeff == 0:
            continue
            
        # Determine the sign for the term
        if i == 0:
            sign = "-" if coeff < 0 else ""
        else:
            sign = " - " if coeff < 0 else " + "
            
        abs_coeff = abs(coeff)
        
        # Build the string for the term
        term_str = ""
        # Coefficient part (don't show '1' unless it's a constant term)
        if abs_coeff != 1 or power == 0:
            term_str += str(abs_coeff)
        
        # Variable part
        if power > 0:
            # Add separator for multi-part terms
            if term_str:
                term_str += " \\, "
            
            if power == 1:
                term_str += "p"
            else:
                term_str += f"p^{{{power}}}"
                
        terms.append(sign + term_str)
        
    return f"f_{n}(p) = {''.join(terms)}"

# We will demonstrate the formula by calculating f_18(p).
n_to_calculate = 18
f_n_polynomial = find_polynomial_expression(n_to_calculate)
final_equation_string = format_polynomial_for_output(f_n_polynomial, n_to_calculate)

# Print the final result, which will be the full equation for f_18(p)
# with all its numeric coefficients, formatted as requested.
print(final_equation_string)