import sympy as sp

def generate_and_print_polynomial(n):
    """
    Generates and prints the polynomial f_n(p) using the derived formula.
    
    The formula is f_n(p) = (p^n - (1-p)^n) / (2p - 1).
    
    Args:
        n (int): The index of the polynomial to generate.
    """
    # Define the symbolic variable 'p'
    p = sp.Symbol('p')

    # The closed-form formula for f_n(p)
    # The division results in a polynomial because (x-y) is a factor of (x^n - y^n)
    f_n_expression = (p**n - (1 - p)**n) / (2 * p - 1)

    # Expand the expression to get the polynomial form
    polynomial = sp.expand(f_n_expression)

    # Convert the symbolic polynomial to a string for custom formatting
    poly_str = str(polynomial)

    # Format the string to be more readable and match the style in the problem description
    # Replace Python's power operator '**' with '^'
    # Add spaces around the multiplication operator '*'
    formatted_poly_str = poly_str.replace('**', '^').replace('*', ' * ')
    
    # Print the final equation with all its numbers (coefficients and powers)
    print(f"f_{n}(p) = {formatted_poly_str}")

# We will now compute f_20(p) as a demonstration of the formula.
generate_and_print_polynomial(20)
