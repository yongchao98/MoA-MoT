import numpy as np
from scipy.integrate import quad

def p1(x):
  """
  Function definition for the interval 0 <= x <= 3.
  """
  return (2 * x**3) / 8

def p2(x):
  """
  Function definition for the interval 3 <= x <= 5.
  """
  # Note: np.exp, np.sin, and np.cos are used for vectorized operations,
  # which is good practice with numerical libraries.
  return (np.exp(x) * (1 + np.sin(x))) / (1 + np.cos(x))

# --- Calculation ---

# Calculate the integral for the first part (from 0 to 3)
# The result from quad is a tuple (value, estimated_error)
integral1, error1 = quad(p1, 0, 3)

# Calculate the integral for the second part (from 3 to 4)
integral2, error2 = quad(p2, 3, 4)

# Calculate the total integral by summing the two parts
total_integral = integral1 + integral2

# --- Output ---

print("The integral from 0 to 4 is calculated by splitting it into two parts:")
print(f"1. The integral of (2*x^3)/8 from x=0 to x=3.")
print(f"2. The integral of (e^x * (1 + sin(x))) / (1 + cos(x)) from x=3 to x=4.\n")

print("The value of the first integral is:")
print(f"Part 1 = {integral1:.4f}")

print("\nThe value of the second integral is:")
print(f"Part 2 = {integral2:.4f}")

print("\nThe final equation combining these parts is:")
# We print each number that contributes to the final sum
print(f"{integral1:.4f} + ({integral2:.4f}) = {total_integral:.4f}")