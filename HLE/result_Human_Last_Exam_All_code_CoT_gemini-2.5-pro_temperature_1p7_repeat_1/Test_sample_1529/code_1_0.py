import sympy

def solve():
    """
    Calculates the discriminant of the polynomial defining the curve.
    The curve is y^2 = x^6 + 2x^3 + 4x^2 + 4x + 1.
    We compute the discriminant of the polynomial f(x) = x^6 + 2x^3 + 4x^2 + 4x + 1.
    """
    # Define the variable x
    x = sympy.symbols('x')

    # Define the polynomial
    polynomial = x**6 + 2*x**3 + 4*x**2 + 4*x + 1

    # Calculate the discriminant of the polynomial
    discriminant_value = sympy.discriminant(polynomial, x)

    # Print the result
    print(f"The polynomial is f(x) = x**6 + 2*x**3 + 4*x**2 + 4*x + 1")
    print(f"The minimal discriminant of the curve is the discriminant of this polynomial, which is:")
    print(discriminant_value)

solve()