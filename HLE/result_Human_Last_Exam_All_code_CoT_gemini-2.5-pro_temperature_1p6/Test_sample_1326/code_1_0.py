import numpy as np
from scipy.integrate import quad

def part1(x):
  """Defines the function for 0 <= x <= 3"""
  # This simplifies to x**3 / 4
  return (2 * x**3) / 8.0

def part2(x):
  """Defines the function for 3 < x <= 5"""
  return (np.exp(x) * (1 + np.sin(x))) / (1 + np.cos(x))

# --- Explanation of the method ---
print("To find the integral of the piecewise function from x=0 to x=4, we split the integral at x=3.")
print("The total integral is the sum of two integrals:")
print("1. The integral of (2*x^3)/8 from x=0 to x=3.")
print("2. The integral of (e^x(1+sin(x)))/(1+cos(x)) from x=3 to x=4.")
print("-" * 30)

# --- Calculation of the first integral ---
# The integral is from 0 to 3 for the first part of the function.
integral1, error1 = quad(part1, 0, 3)
print(f"Calculating the first integral...")
print(f"The value of the integral from 0 to 3 is: {integral1}")
print("-" * 30)


# --- Calculation of the second integral ---
# The integral is from 3 to 4 for the second part of the function.
integral2, error2 = quad(part2, 3, 4)
print(f"Calculating the second integral...")
print(f"The value of the integral from 3 to 4 is: {integral2}")
print("-" * 30)

# --- Final Calculation and Output ---
# The total integral is the sum of the two parts.
total_integral = integral1 + integral2
print("The total integral is the sum of the two parts.")
print("Final Equation:")
print(f"{integral1} + ({integral2}) = {total_integral}")

# Output the final numerical answer.
print("\nThe final result of the integral from x=0 to x=4 is:")
print(total_integral)