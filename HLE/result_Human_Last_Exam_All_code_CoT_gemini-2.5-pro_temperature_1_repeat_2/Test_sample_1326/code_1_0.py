import numpy as np
from scipy.integrate import quad

# Define the first part of the function for 0 <= x <= 3
def p1(x):
    return (2 * x**3) / 8

# Define the second part of the function for 3 < x <= 5
def p2(x):
    # Handle the case where cos(x) = -1 to avoid division by zero,
    # although it doesn't occur in the integration range [3, 4].
    if np.isclose(np.cos(x), -1):
        return np.inf
    return (np.exp(x) * (1 + np.sin(x))) / (1 + np.cos(x))

# Set the integration limits
lower_bound_total = 0
mid_point = 3
upper_bound_total = 4

# Calculate the integral for the first part: from 0 to 3
integral1, error1 = quad(p1, lower_bound_total, mid_point)

# Calculate the integral for the second part: from 3 to 4
integral2, error2 = quad(p2, mid_point, upper_bound_total)

# Calculate the total integral
total_integral = integral1 + integral2

# Print the results of each step
print(f"The integral is split into two parts: from {lower_bound_total} to {mid_point} and from {mid_point} to {upper_bound_total}.")
print(f"Integral of the first part (from x={lower_bound_total} to x={mid_point}): {integral1:.4f}")
print(f"Integral of the second part (from x={mid_point} to x={upper_bound_total}): {integral2:.4f}")
print("\nFinal Equation:")
print(f"Total Integral = {integral1:.4f} + ({integral2:.4f})")
print(f"Total Integral = {total_integral:.4f}")

# Final answer in the required format
# print(f"\n<<<{total_integral:.4f}>>>")