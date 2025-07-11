import numpy as np

# Based on the reasoning outlined above, we proceed with the following assumptions:
# 1. The robotic arm's path is given by y1(x) = 42688 / x^6.
# 2. The integration interval is [2, 6], derived from the given points.

# Define the constant C and the integration bounds
C = 42688
a = 2.0
b = 6.0

# The function to integrate is y1(x) = C * x^(-6).
# The antiderivative (indefinite integral) of y1(x) is F(x) = C * x^(-5) / (-5).
def antiderivative_F(x):
  """Calculates the value of the antiderivative of y1(x)."""
  return -C / (5.0 * x**5)

# Calculate the value of the definite integral F(b) - F(a)
val_at_b = antiderivative_F(b)
val_at_a = antiderivative_F(a)
integral_result = val_at_b - val_at_a

# The problem asks to output each number in the final equation.
print("Calculating the definite integral of y1(x) = 42688 / x^6 from x=2 to x=6.")
print("Antiderivative F(x) = -42688 / (5 * x^5)")
print(f"Value at upper bound F(6) = -42688 / (5 * 6^5) = {val_at_b:.7f}")
print(f"Value at lower bound F(2) = -42688 / (5 * 2^5) = {val_at_a:.7f}")
print("\nFinal Equation:")
# The equation is F(b) - F(a) = result
print(f"({val_at_b:.7f}) - ({val_at_a:.7f}) = {integral_result:.7f}")
