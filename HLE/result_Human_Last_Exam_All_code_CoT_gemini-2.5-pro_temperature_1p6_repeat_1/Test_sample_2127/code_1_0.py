import sympy
from sympy import symbols, E, pi, log, tanh, cos, sin, cosh, sinh, tan, sec, Integer

def solve_maclaurin_coefficient():
    """
    This function calculates the 4th Maclaurin series coefficient of the given function
    using the sympy library.
    """
    x = symbols('x')

    # To ensure exact fractions, we use sympy.Integer
    frac_9 = Integer(9)
    frac_16 = Integer(16)
    frac_4 = Integer(4)
    frac_5 = Integer(5)
    frac_6 = Integer(6)
    frac_1 = Integer(1)
    frac_2 = Integer(2)

    # First term of the expression
    t1 = (frac_9 * x**4) / (frac_16 * E)

    # Numerator of the second term
    n1 = x**4 - (frac_5/frac_6) * log(x**4 + 1)**2
    n2 = sympy.exp(tanh(x**3)/frac_2) - 1
    n3 = cos(sin(pi * cosh(x**6))) - 1/E
    numerator = frac_4 * n1 * n2 * n3

    # Denominator of the second term
    d1 = tan(x**6) - log(x**8 + 1)
    d2 = sympy.exp(cos(x**5)**2 + sinh(x**2)) - 1
    d3 = cosh(x**3) - sec(x**7)
    denominator = d1 * d2 * d3
    
    # Second term of the expression
    t2 = numerator / denominator

    # The full function
    f = t1 + t2

    # We need the series expansion up to O(x^5) to get the x^4 coefficient.
    # The third argument '5' in series() means the series is computed up to, 
    # but not including, the term x^5.
    f_series = f.series(x, 0, 5)

    # The series is a Laurent series. We remove the O() term to get the polynomial part.
    f_poly = f_series.removeO()

    # Extract the coefficient of x^4
    coeff_4 = f_poly.coeff(x, 4)

    # For clarity in the output, let's also compute the coefficient for each part.
    c1 = t1.coeff(x, 4)
    
    # It's better to compute the coefficient of t2 from the combined series
    # to avoid recomputing the series expansion.
    # The coefficient of x^4 in t2 is the total coefficient minus the one from t1.
    c2 = coeff_4 - c1
    
    print("Finding the 4th Maclaurin series coefficient of the function.")
    print("Let the function be f(x) = T1(x) + T2(x).")
    print("Coefficient of x^4 from the first term T1(x) = 9*x**4/(16*E):")
    print(f"{c1}")
    print("\nCoefficient of x^4 from the complex second term T2(x):")
    print(f"{c2}")
    print("\nThe total coefficient of x^4 is the sum of the two:")
    print(f"{c1} + ({c2}) = {coeff_4}")
    print("\nEvaluating the expression numerically:")
    print(f"{c1.evalf()} + ({c2.evalf()}) = {coeff_4.evalf()}")

solve_maclaurin_coefficient()