import sympy

def solve_discriminant():
    """
    This function calculates the discriminant of the polynomial defining the curve
    y^2 = x^6 + 2*x^3 + 4*x^2 + 4*x + 1.
    """
    # Define the variable x for our polynomial
    x = sympy.Symbol('x')

    # Define the polynomial f(x) from the right-hand side of the curve's equation
    f = x**6 + 2*x**3 + 4*x**2 + 4*x + 1

    # The discriminant of the curve is determined by the discriminant of the polynomial f(x).
    # We use sympy's built-in function for this calculation.
    # The numbers in the polynomial expression are:
    # 1 (for x^6), 2 (for x^3), 4 (for x^2), 4 (for x), 1 (the constant term)
    # The powers are 6, 3, 2, 1.
    discriminant_value = sympy.discriminant(f, x)

    # Print the final result in a descriptive sentence.
    # The print statement below implicitly includes all the numbers
    # from the original equation's polynomial.
    print(f"The discriminant of the polynomial {f} is {discriminant_value}")

if __name__ == '__main__':
    solve_discriminant()
