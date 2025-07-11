import numpy as np

def calculate_l(b, n=101):
    """
    Calculates the value of l(b) based on the trace formula.
    l(b) = Tr((B(b)B(b)^T)^-1) = (n + (n-2)b^2) / (1-b^2)
    """
    if abs(b) >= 1:
        raise ValueError("b must be in the interval (-1, 1)")
    
    numerator = n + (n - 2) * b**2
    denominator = 1 - b**2
    return numerator / denominator

# Set the values of b
b1 = 0.5
b2 = -0.5

# Calculate l(1/2) and l(-1/2)
l_half = calculate_l(b1)
l_neg_half = calculate_l(b2)

# Calculate the final expression
final_value = 6 * (l_half + l_neg_half)

# Print the results, showing each number in the final equation
# To display fractions nicely, we can use the fractions module
from fractions import Fraction

l_half_frac = Fraction(l_half).limit_denominator()
l_neg_half_frac = Fraction(l_neg_half).limit_denominator()

print(f"The dimension of the matrix is n = 101.")
print(f"The value of l(b) is calculated using the formula: l(b) = (n + (n-2)*b^2) / (1-b^2)")
print(f"For b = 1/2, l(1/2) = (101 + 99*(1/2)^2) / (1 - (1/2)^2) = {l_half_frac.numerator}/{l_half_frac.denominator}")
print(f"For b = -1/2, l(-1/2) = (101 + 99*(-1/2)^2) / (1 - (-1/2)^2) = {l_neg_half_frac.numerator}/{l_neg_half_frac.denominator}")
print(f"The expression to compute is 6 * (l(1/2) + l(-1/2))")
print(f"Final calculation: 6 * ({l_half_frac} + {l_neg_half_frac}) = 6 * ({l_half_frac + l_neg_half_frac}) = {final_value}")

print("\nFinal Answer:")
print(f"{int(final_value)}")