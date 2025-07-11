import math
from fractions import Fraction

def calculate_moment_fraction(k):
    """
    Calculates the k-th moment <g^k> of the arcsine distribution as an exact fraction.
    The formula is <g^k> = (product_{i=0}^{k-1}(2i+1)) / (2^k * k!).
    """
    if not isinstance(k, int) or k < 0:
        raise ValueError("Moment order k must be a non-negative integer.")
    if k == 0:
        return Fraction(1)
    
    # The numerator is the product of the first k odd numbers, which is (2k-1)!!.
    numerator = 1
    for i in range(k):
        numerator *= (2 * i + 1)
        
    # The denominator is 2^k * k!.
    denominator = (2**k) * math.factorial(k)
    
    return Fraction(numerator, denominator)

# Calculate the first moment (average value) <g> for k=1
avg_g = calculate_moment_fraction(1)

# Calculate the fourth moment <g^4> for k=4
moment_g4 = calculate_moment_fraction(4)

# Calculate the required ratio
ratio = moment_g4 / avg_g

# Output the results, including the numbers in the final equation as requested
print(f"The average conductance is <g> = {avg_g}")
print(f"The fourth statistical moment of conductance is <g^4> = {moment_g4}")
print("\nThe ratio <g^4> / <g> is calculated as follows:")
print(f"{moment_g4} / {avg_g} = {ratio}")
print(f"\nThe numerical value of the ratio is {float(ratio):.6f}")