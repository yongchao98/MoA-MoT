import math
from fractions import Fraction

def get_moment(k):
    """
    Calculates the k-th moment of the dimensionless conductance 
    based on the universal formula for the BDI (LDOS) / D (conductance) class criticality.
    The formula is <g^k> = (2k)! / (4^k * (k!)^2).
    """
    # Use Fraction for exact arithmetic
    numerator = Fraction(math.factorial(2 * k))
    denominator = Fraction(4**k * (math.factorial(k))**2)
    return numerator / denominator

# Calculate the average value (k=1)
avg_g = get_moment(1)

# Calculate the fourth moment (k=4)
g_moment_4 = get_moment(4)

# Calculate the ratio
ratio = g_moment_4 / avg_g

# Print the results in an equation format
print(f"The average value (first moment) is <g> = {avg_g}")
print(f"The fourth statistical moment is <g^4> = {g_moment_4}")
print(f"The ratio is <g^4> / <g> = ({g_moment_4}) / ({avg_g})")
print(f"Final Result: {ratio}")
