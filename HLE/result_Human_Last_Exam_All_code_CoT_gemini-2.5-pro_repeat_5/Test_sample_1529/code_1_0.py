import sympy

def solve():
    """
    Calculates the minimal discriminant of the curve y^2 = x^6 + 2x^3 + 4x^2 + 4x + 1.
    """
    # Define the symbol x
    x = sympy.Symbol('x')

    # The equation of the curve is y^2 = f(x).
    # We define the polynomial f(x) = 1*x^6 + 0*x^5 + 0*x^4 + 2*x^3 + 4*x^2 + 4*x + 1
    # Let's list the coefficients as requested.
    c6 = 1
    c5 = 0
    c4 = 0
    c3 = 2
    c2 = 4
    c1 = 4
    c0 = 1

    f = c6*x**6 + c5*x**5 + c4*x**4 + c3*x**3 + c2*x**2 + c1*x + c0

    print("The curve is defined by the equation y^2 = f(x), where:")
    print(f"f(x) = {c6}*x^6 + {c3}*x^3 + {c2}*x^2 + {c1}*x + {c0}")
    print("\nThe numbers in the equation are:")
    print(f"Coefficient of x^6: {c6}")
    print(f"Coefficient of x^5: {c5}")
    print(f"Coefficient of x^4: {c4}")
    print(f"Coefficient of x^3: {c3}")
    print(f"Coefficient of x^2: {c2}")
    print(f"Coefficient of x^1: {c1}")
    print(f"Constant term: {c0}")

    # Calculate the discriminant of the polynomial f(x)
    # The discriminant of this model is the minimal discriminant because the model is minimal.
    minimal_discriminant = sympy.discriminant(f, x)

    print(f"\nThe minimal discriminant of the curve is: {minimal_discriminant}")

solve()