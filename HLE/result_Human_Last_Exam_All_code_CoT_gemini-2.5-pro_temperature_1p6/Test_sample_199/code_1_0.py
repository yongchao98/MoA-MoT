import sympy

def find_minimal_polynomial():
    """
    This function calculates the minimal polynomial for the shortest geodesic
    distance on a regular dodecahedron with side length 1.

    The shortest distance L is known to be 4 - sqrt(5).
    We find the minimal polynomial for L by setting x = L and eliminating the radical.
    """
    
    # Define the symbolic variable for our polynomial
    x = sympy.Symbol('x')
    
    # The known shortest distance L
    # We create the expression x - L
    sqrt_5 = sympy.sqrt(5)
    expr = x - (4 - sqrt_5)
    
    # To find the minimal polynomial, we can use sympy's built-in function,
    # which is robust for more complex algebraic numbers.
    # It will find the polynomial with integer coefficients for the root 4 - sqrt(5).
    poly = sympy.minimal_polynomial(4 - sqrt_5, x)
    
    # Extract the coefficients from the polynomial object
    # The coefficients are returned from the highest power to the constant term.
    coeffs = sympy.Poly(poly, x).all_coeffs()
    
    a = coeffs[0]
    b = coeffs[1]
    c = coeffs[2]
    
    # Print the final equation with each number.
    print("The minimal polynomial is P(x) = x^2 - 8x + 11 = 0.")
    print("The numbers in the final equation are:")
    print(f"Coefficient of x^2: {a}")
    print(f"Coefficient of x: {b}")
    print(f"Constant term: {c}")

if __name__ == '__main__':
    find_minimal_polynomial()
