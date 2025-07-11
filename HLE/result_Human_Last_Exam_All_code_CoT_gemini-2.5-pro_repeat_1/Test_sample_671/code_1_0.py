import sympy

def generate_polynomial_string(n):
    """
    Generates the polynomial f_n(p) using the derived formula,
    and returns it as a formatted string.
    """
    p = sympy.Symbol('p')

    # Handle the special case p=1/2 where the denominator is zero.
    # The formula has a removable singularity at p=1/2.
    # We can compute the limit, or let sympy's `cancel` handle it.
    if n <= 0:
        return "n must be a positive integer"
        
    # The derived formula for f_n(p)
    formula = (p**n - (1 - p)**n) / (2*p - 1)

    # Use sympy.cancel to simplify the expression, which correctly handles the division
    poly_expr = sympy.cancel(formula)

    # Convert the expression to a polynomial object to easily access coefficients
    poly_obj = sympy.Poly(poly_expr, p)
    
    # Get the coefficients and degree
    coeffs = poly_obj.all_coeffs()
    degree = poly_obj.degree()

    # Build the string representation as requested
    terms = []
    for i, coeff in enumerate(coeffs):
        power = degree - i
        if coeff == 0:
            continue
        
        # Format coefficient
        term_str = str(coeff)
        
        # Format variable part
        if power > 0:
            if term_str == '1':
                term_str = "" # Omit coefficient if it's 1
            elif term_str == '-1':
                term_str = "-" # Omit coefficient if it's -1
            else:
                 term_str += " * "
            
            if power == 1:
                term_str += "p"
            else:
                term_str += f"p**{power}"
        
        terms.append(term_str)
    
    # Join all terms with ' + ' and clean up formatting
    final_string = " + ".join(terms).replace("+ -", "- ")
    
    return final_string

# --- Main execution ---
# We will calculate f_15(p) to verify against the provided list.
n = 15
fn_poly_str = generate_polynomial_string(n)

# The problem description asks to output each number in the final equation.
# Our generated string represents the polynomial with all its non-zero coefficients.
# To satisfy the prompt more directly, let's also print the equation with coefficients.
p = sympy.Symbol('p')
formula = (p**n - (1 - p)**n) / (2*p - 1)
poly_expr = sympy.cancel(formula)
poly_obj = sympy.Poly(poly_expr, p)
coeffs = poly_obj.all_coeffs()
degree = poly_obj.degree()

equation_parts = []
for i, coeff in enumerate(coeffs):
    power = degree - i
    if power > 1:
        equation_parts.append(f"{coeff} * p**{power}")
    elif power == 1:
        equation_parts.append(f"{coeff} * p")
    elif power == 0:
        equation_parts.append(f"{coeff}")

print(f"f_{n}(p) = " + " + ".join(equation_parts))
