import sympy

def solve_discriminant():
    """
    Calculates the discriminant of the polynomial from the curve equation.
    The problem asks for the minimal discriminant of the curve defined by
    y^2 = x^6 + 2x^3 + 4x^2 + 4x + 1.
    This is interpreted as finding the discriminant of the polynomial on the
    right-hand side.
    """
    # Define the symbolic variable
    x = sympy.Symbol('x')

    # Define the polynomial f(x)
    # The equation is y^2 = 1*x^6 + 0*x^5 + 0*x^4 + 2*x^3 + 4*x^2 + 4*x + 1
    c6, c5, c4, c3, c2, c1, c0 = 1, 0, 0, 2, 4, 4, 1
    f = c6*x**6 + c5*x**5 + c4*x**4 + c3*x**3 + c2*x**2 + c1*x + c0

    # Per the instructions, we output the numbers in the final equation.
    # We will print the equation of the curve explicitly.
    print(f"The curve is defined by the equation: y^2 = {c6}*x^6 + {c3}*x^3 + {c2}*x^2 + {c1}*x + {c0}")

    # Calculate the discriminant of the polynomial f(x)
    discriminant_value = sympy.discriminant(f, x)

    # Print the result
    print("The minimal discriminant of the curve is:")
    print(discriminant_value)

if __name__ == "__main__":
    solve_discriminant()