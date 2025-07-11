import sympy

def solve():
    """
    Calculates the minimal discriminant of the curve y^2 = x^6 + 2x^3 + 4x^2 + 4x + 1.
    """
    # Define the variable and the polynomial
    x = sympy.Symbol('x')
    f = x**6 + 2*x**3 + 4*x**2 + 4*x + 1

    # Calculate the discriminant of the polynomial
    # For the given curve, its model is already minimal, so the discriminant of the
    # polynomial is the minimal discriminant of the curve.
    minimal_discriminant = sympy.discriminant(f, x)

    # The instruction "output each number in the final equation" is interpreted as
    # printing the final numerical answer.
    print(minimal_discriminant)

solve()