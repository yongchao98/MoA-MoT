import sympy as sp
from sympy import sin, cos, tan, sinh, cosh, log, exp, pi

def solve_maclaurin_coefficient():
    """
    This function calculates the 4th Maclaurin series coefficient of the given complex function
    using the sympy library for symbolic computation.
    """
    # Set up the symbolic variable
    x = sp.Symbol('x')

    # Define the function f(x) as given in the problem.
    # The function is a sum of two terms, f(x) = T1(x) + T2(x).

    # First term
    T1 = (9 * x**4) / (16 * sp.E)

    # Second term is a fraction. We define its numerator and denominator separately.

    # Numerator of the second term
    numerator_factor1 = x**4 - sp.Rational(5, 6) * log(x**4 + 1)**2
    numerator_factor2 = exp(tanh(x**3) / 2) - 1
    numerator_factor3 = cos(sin(pi * cosh(x**6))) - 1 / sp.E
    numerator_T2 = 4 * numerator_factor1 * numerator_factor2 * numerator_factor3

    # Denominator of the second term
    denominator_factor1 = tan(x**6) - log(x**8 + 1)
    denominator_factor2 = exp(cos(x**5)**2 + sinh(x**2)) - 1
    denominator_factor3 = cosh(x**3) - 1/cos(x**7)  # Using 1/cos(x) for sec(x)
    denominator_T2 = denominator_factor1 * denominator_factor2 * denominator_factor3

    # The second term
    T2 = numerator_T2 / denominator_T2

    # The full function
    f = T1 + T2

    # The 4th Maclaurin coefficient is the coefficient of x**4 in the series expansion.
    # We compute the Laurent series up to O(x^5) to ensure we can capture the x^4 term.
    series_expansion = f.series(x, 0, 5)

    # Extract the coefficient of x**4 from the series expansion.
    # .removeO() gets rid of the O(x^n) term.
    # .coeff(x, 4) extracts the coefficient of x**4.
    c4 = series_expansion.removeO().coeff(x, 4)

    # The problem asks to output each number in the final equation.
    # Printing the symbolic expression c4 achieves this.
    print("The symbolic expression for the 4th Maclaurin coefficient (c_4) is:")
    print(c4)

    # For completeness, we also print the numerical value.
    print("\nThe numerical value of the coefficient is:")
    print(c4.evalf())

# Execute the function to find and print the coefficient.
solve_maclaurin_coefficient()