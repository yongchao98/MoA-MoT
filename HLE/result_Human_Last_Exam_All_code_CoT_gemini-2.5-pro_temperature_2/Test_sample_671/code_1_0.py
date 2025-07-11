import sympy

def generate_polynomial_string(n):
    """
    Calculates the polynomial f_n(p) based on the derived formula and returns it
    as a formatted string.
    
    The formula is f_n(p) = (p^n - (1-p)^n) / (2p - 1).
    """
    if not isinstance(n, int) or n <= 0:
        return "Error: n must be a positive integer."

    p = sympy.Symbol('p')

    # Handle the first two special cases directly.
    if n <= 2:
        return f"f_{n}(p) = 1"

    # Use sympy to perform symbolic expansion and polynomial division.
    numerator = sympy.expand(p**n - (1-p)**n)
    denominator = 2*p - 1
    
    # sympy.div returns (quotient, remainder)
    fn_poly, remainder = sympy.div(numerator, denominator, domain='ZZ')

    if remainder != 0:
        # This case should not be reached if the formula is correct.
        return f"Error: Calculation for f_{n}(p) resulted in a non-zero remainder."
    
    # Convert the sympy expression to a Poly object for easier coefficient access.
    poly = sympy.Poly(fn_poly, p)
    degree = poly.degree()
    
    result_parts = []
    is_first_term = True
    # Iterate from the highest degree to the constant term.
    for d in range(degree, -1, -1):
        coeff = int(poly.coeff_monomial(p**d))
        
        # Skip terms with a zero coefficient.
        if coeff == 0:
            continue
        
        term_str = ""
        # Determine the sign for the term.
        if is_first_term:
            if coeff < 0:
                term_str += "- "
        else:
            term_str += " + " if coeff > 0 else " - "

        abs_coeff = abs(coeff)
        
        # Append the coefficient's absolute value, unless it's 1 for a non-constant term.
        if abs_coeff != 1 or d == 0:
            term_str += str(abs_coeff)

        # Append the variable part 'p' if the degree is > 0.
        if d > 0:
            # Add a multiplication sign if there's a visible coefficient.
            if abs_coeff != 1:
                term_str += " * "
            
            # Format the power of p.
            if d > 1:
                term_str += f"p^{d}"
            else:
                term_str += "p"
                
        result_parts.append(term_str)
        is_first_term = False
        
    print(f"f_{n}(p) = {''.join(result_parts)}")


if __name__ == '__main__':
    # Demonstrate by calculating f_15(p) as shown in the problem description.
    example_n = 15
    generate_polynomial_string(example_n)