import math
from fractions import Fraction

def factorial(n):
    """Computes the factorial of a non-negative integer."""
    if not isinstance(n, int) or n < 0:
        raise ValueError("Factorial is only defined for non-negative integers")
    return math.factorial(n)

def double_factorial(n):
    """Computes the double factorial of a non-negative integer."""
    if not isinstance(n, int) or n < 0:
        raise ValueError("Double factorial is only defined for non-negative integers")
    if n == 0:
        return 1
    result = 1
    for i in range(n, 0, -2):
        result *= i
    return result

def get_conductance_moment(n):
    """
    Calculates the n-th moment of the dimensionless conductance <g^n> for a class D
    disordered wire at the critical point.
    Formula: <g^n> = (n! * 2^n) / ((2n+1)!!)
    """
    numerator = factorial(n) * (2**n)
    denominator = double_factorial(2 * n + 1)
    return Fraction(numerator, denominator)

# Calculate the average conductance <g>, which is the first moment (n=1)
avg_g = get_conductance_moment(1)

# Calculate the fourth moment of the conductance <g^4> (n=4)
moment4_g = get_conductance_moment(4)

# Calculate the final ratio
ratio = moment4_g / avg_g

# Print the values and the final equation
print(f"The average dimensionless conductance is <g> = {avg_g}")
print(f"The fourth moment of the conductance is <g^4> = {moment4_g}")
print("\nThe ratio between the fourth moment and its average value is calculated as follows:")
print(f"<g^4> / <g> = ({moment4_g}) / ({avg_g}) = {ratio}")
