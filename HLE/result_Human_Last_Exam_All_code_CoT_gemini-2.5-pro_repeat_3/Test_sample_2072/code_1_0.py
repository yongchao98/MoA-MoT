import math

def solve_phi_n():
    """
    This function calculates the value of phi(n) based on the derived formula.
    The problem states n >= 5. Since n is not specified, we use n=5.
    """
    # Set n to the smallest possible value as per the problem statement.
    n = 5

    # The derived formula for phi(n) is exp(2n - 4 + 2/n).
    # We calculate the terms of the exponent first.
    term1 = 2 * n
    term2 = -4
    term3 = 2 / n

    # Calculate the exponent value.
    exponent = term1 + term2 + term3

    # Calculate the final result.
    result = math.exp(exponent)

    # As requested, output the numbers in the final equation.
    print(f"For n = {n}:")
    print(f"phi({n}) = exp(2 * {n} + ({term2}) + {term3})")
    print(f"phi({n}) = exp({exponent})")
    print(f"Result: {result}")

solve_phi_n()