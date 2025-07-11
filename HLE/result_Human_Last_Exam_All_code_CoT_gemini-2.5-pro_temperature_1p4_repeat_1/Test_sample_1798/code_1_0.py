import numpy as np
from scipy.integrate import quad

# Step 1: Define the parameters of the circle from the plan
h = 5.5
k = -1.5
r_squared = 22.5

# The equation of the circle is (x - 5.5)^2 + (y - (-1.5))^2 = 22.5

# Define the function y = f(x) for the upper part of the circle
def f(x):
  """
  Represents the equation of the upper semi-circle y=f(x).
  """
  return np.sqrt(r_squared - (x - h)**2) + k

# Step 2: Calculate the normalization constant alpha
# The integral of the PDF alpha * f(x) from 1 to 10 must be 1.
# So, alpha = 1 / (integral of f(x) from 1 to 10)
integral_total, _ = quad(f, 1, 10)
alpha = 1 / integral_total

# Step 3: Calculate the probability P(X < 3)
# P(X < 3) = integral of alpha * f(x) from 1 to 3
integral_partial, _ = quad(f, 1, 3)
probability = alpha * integral_partial

# Print the results
print("The equation of the circle that passes through points A, B, and C is:")
print(f"(x - {h})^2 + (y + {-k})^2 = {r_squared}")
print(f"\nThe value of alpha is: {alpha}")
print(f"The value of P(X < 3) is: {probability}")

# The final answer in the required format will be generated based on these results.