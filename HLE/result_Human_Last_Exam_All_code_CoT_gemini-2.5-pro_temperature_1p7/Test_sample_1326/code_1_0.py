import numpy as np
from scipy.integrate import quad

# Define the two parts of the function p(x)
def part1(x):
  """Function definition for 0 <= x <= 3"""
  return (2 * x**3) / 8

def part2(x):
  """Function definition for 3 <= x <= 5"""
  # Ensure all math operations use numpy for vectorization and precision
  return (np.exp(x) * (1 + np.sin(x))) / (1 + np.cos(x))

# --- Step 1: Calculate the integral of the first part from 0 to 3 ---
# The quad function returns a tuple (result, error_estimate)
# We are interested in the result, which is the first element
integral1_val, error1 = quad(part1, 0, 3)

# --- Step 2: Calculate the integral of the second part from 3 to 4 ---
integral2_val, error2 = quad(part2, 3, 4)

# --- Step 3: Sum the results to get the total integral ---
total_integral = integral1_val + integral2_val

# --- Display the results ---
print("The integral is split into two parts at x = 3.")
print(f"Integral from 0 to 3: {integral1_val:.8f}")
print(f"Integral from 3 to 4: {integral2_val:.8f}")
print("\nThe total integral is the sum of these two values.")
print(f"Final Equation: {integral1_val:.8f} + {integral2_val:.8f} = {total_integral:.8f}")

# The problem requests the final numerical answer in a specific format at the end.
# print(f"\n<<<final_answer>>>")
# print(f"<<<{total_integral:.8f}>>>")