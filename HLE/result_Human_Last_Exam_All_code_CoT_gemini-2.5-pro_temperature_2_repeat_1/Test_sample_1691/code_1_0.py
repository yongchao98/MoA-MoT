import numpy as np

def solve_integral_approximation():
    """
    Calculates the coefficients for the analytical approximation of the integral
    I(epsilon) for small epsilon.
    The approximation has the form C * epsilon^exponent.
    """

    # From the denominator g(x), the dominant term at x=0 is a_n * x^n
    n = 5.0
    a_n = 9.0

    # The exponent of epsilon in the approximation is (1/n) - 1
    exponent = (1 / n) - 1

    # The coefficient C is given by (a_n^(-1/n)) * integral_part
    # where the integral part is (pi/n) / sin(pi/n)
    integral_part = (np.pi / n) / np.sin(np.pi / n)
    coefficient = (a_n**(-1 / n)) * integral_part

    # Print the final analytical formula with the calculated numerical values
    print("The analytical formula for I(epsilon) in the small epsilon regime is:")
    print(f"I(epsilon) ≈ {coefficient} * epsilon^({exponent})")
    # For a slightly cleaner look, we can also present it like this:
    p = (n - 1.0) / n
    print(f"I(epsilon) ≈ {coefficient} / (epsilon^{p})")


solve_integral_approximation()

# The question can be interpreted as asking for the coefficient C in the formula I(epsilon) ≈ C * epsilon^(-4/5).
n = 5.0
a_n = 9.0
coefficient = (a_n**(-1/n)) * (np.pi/n) / np.sin(np.pi/n)
# print(f"<<<{coefficient}>>>")