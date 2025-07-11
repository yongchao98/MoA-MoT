import sympy
from sympy import sin, cos, pi, piecewise, integrate, Symbol, pprint

def solve():
    """
    This function calculates the smallest possible constant C.
    """
    # Define the symbolic variable
    x = Symbol('x')

    # As reasoned in the plan, the extremal functions are likely:
    # f(x) = sin(x)
    # a(x) is a "bang-bang" function that maximizes the ratio.
    f = sin(x)
    f_prime = cos(x)

    # a(x) = 3 where sin^2(x) > cos^2(x), and a(x) = 1 where sin^2(x) < cos^2(x)
    # The boundaries are where |tan(x)| = 1, i.e., x = pi/4, 3*pi/4, etc.
    a = piecewise(
        (3, (x > pi/4) & (x < 3*pi/4)),
        (3, (x > 5*pi/4) & (x < 7*pi/4)),
        (1, True)
    )

    # First, verify that the condition ∫ a(x)f(x) dx = 0 is met.
    # Note: sin(x) is odd and a(x) is even around x=pi, on the interval [0, 2pi].
    # Their product integrated over [0, 2pi] is zero.
    # Let's compute it to be sure.
    constraint_integral = integrate(a * f, (x, 0, 2*pi))

    # Calculate the numerator: ∫ a(x) f(x)^2 dx
    numerator = integrate(a * f**2, (x, 0, 2*pi))

    # Calculate the denominator: ∫ a(x) f'(x)^2 dx
    denominator = integrate(a * f_prime**2, (x, 0, 2*pi))

    # The constant C is the ratio of the numerator and denominator
    C = numerator / denominator

    print("The extremal functions chosen are:")
    print("f(x) = sin(x)")
    print("a(x) = 3 if |tan(x)| > 1, and 1 if |tan(x)| < 1.")
    print("\nConstraint check ∫ a(x)f(x) dx from 0 to 2π:")
    print(f"Value = {constraint_integral}")
    
    print("\nCalculation of the constant C:")
    print(f"Numerator: ∫ a(x)f(x)² dx = {numerator}")
    print(f"Denominator: ∫ a(x)f'(x)² dx = {denominator}")
    print(f"C = Numerator / Denominator = ({numerator}) / ({denominator})")

    # The final expression for C
    final_C_expression = sympy.simplify(C)
    print("\nSimplified expression for C:")
    pprint(final_C_expression)

    # Numerical value
    print(f"\nNumerical value of C: {final_C_expression.evalf()}")


solve()