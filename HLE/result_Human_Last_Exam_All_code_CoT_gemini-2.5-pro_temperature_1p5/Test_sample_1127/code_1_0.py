import sympy

def solve_minimal_polynomial():
    """
    This function finds and prints the minimal polynomial for the connective constant.
    The problem leads to the equation x^2 - 2x - 1 = 0.
    """
    x = sympy.Symbol('x')
    
    # The minimal polynomial for the connective constant mu is mu^2 - 2*mu - 1 = 0.
    # We represent this as a polynomial in x.
    poly = x**2 - 2*x - 1
    
    # The coefficients of the polynomial
    coeffs = sympy.Poly(poly, x).all_coeffs()
    
    # We want to print the equation with its coefficients.
    # The format will be c0*x^2 + c1*x + c2 = 0
    
    print("The minimal polynomial is P(x) = 0, where P(x) is:")
    
    equation_parts = []
    
    # The equation is 1*x^2 - 2*x - 1 = 0
    # Let's print each term explicitly.
    print(f"{coeffs[0]}*x^2 + ({coeffs[1]})*x + ({coeffs[2]}) = 0")

solve_minimal_polynomial()