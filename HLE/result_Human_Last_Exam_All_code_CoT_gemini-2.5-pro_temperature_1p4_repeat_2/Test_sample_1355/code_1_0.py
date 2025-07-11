import math
from fractions import Fraction

def double_factorial(n):
    """Computes the double factorial of an integer n."""
    if n < 0:
        return 0
    result = 1
    for i in range(n, 0, -2):
        result *= i
    return result

def get_conductance_moment(n):
    """
    Calculates the n-th moment of the dimensionless conductance <g^n>
    for a disordered wire at the BDI critical point.
    The formula is <g^n> = (2n - 1)!! / (2n)!!.
    """
    if n < 1:
        raise ValueError("Moment number n must be a positive integer.")
    
    numerator = double_factorial(2 * n - 1)
    denominator = double_factorial(2 * n)
    
    return Fraction(numerator, denominator)

# Calculate the average conductance, which is the 1st moment (n=1)
avg_g = get_conductance_moment(1)

# Calculate the fourth moment of the conductance (n=4)
fourth_moment_g = get_conductance_moment(4)

# Calculate the ratio between the fourth and the first moment
ratio = fourth_moment_g / avg_g

# Print the components of the final equation and the result
print("This script calculates the ratio of the fourth to the first statistical moment of conductance for a disordered Majorana wire at criticality.")
print("-" * 50)
print("The average conductance <g> is the first moment (n=1).")
print(f"Its value is: <g> = {avg_g.numerator} / {avg_g.denominator}")
print("\nThe fourth statistical moment of conductance <g^4> is for n=4.")
print(f"Its value is: <g^4> = {fourth_moment_g.numerator} / {fourth_moment_g.denominator}")
print("\nThe desired ratio is <g^4> / <g>.")
print(f"Ratio = ({fourth_moment_g.numerator}/{fourth_moment_g.denominator}) / ({avg_g.numerator}/{avg_g.denominator}) = {ratio.numerator}/{ratio.denominator}")
print(f"The final numerical value is: {float(ratio)}")
print("-" * 50)
