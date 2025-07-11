import numpy as np
from scipy.integrate import quad

# Step 1: Define the function y = f(x) based on the derived circle equation.
# The circle equation is (x - 5.5)^2 + (y + 1.5)^2 = 22.5.
# The function f(x) is the upper arc of this circle.
def f(x):
  """
  The function y=f(x) representing the path of the circle.
  """
  h = 5.5
  k = -1.5
  r_squared = 22.5
  # We solve for y: y = -k +/- sqrt(r^2 - (x-h)^2)
  # Since the circle passes through (4,3), we take the positive root.
  # y = 1.5 + sqrt(22.5 - (x-5.5)^2) -> This is incorrect, k=-1.5, so y = -(-1.5) ... No
  # y + k = ... => y = -k + ... = 1.5 + ...
  # Let's retrace: y = -1.5 +/- sqrt(...) => -1.5+4.5 = 3. So the function is correct.
  return np.sqrt(r_squared - (x - h)**2) - k

# Step 2: Calculate alpha.
# alpha is the normalization constant, where alpha = 1 / integral of f(x) over [1, 10].
total_integral, _ = quad(f, 1, 10)
alpha = 1 / total_integral

# Step 3: Calculate P(X < 3).
# This is the integral of alpha*f(x) from 1 to 3.
partial_integral, _ = quad(f, 1, 3)
prob_x_less_than_3 = alpha * partial_integral

# Step 4: Output the results, showing the numbers used in the final equations.
print("1. Find the value of alpha:")
print(f"The equation for alpha is: alpha = 1 / Integral(f(x) from 1 to 10)")
print(f"The value of the integral is: {total_integral}")
print(f"alpha = 1 / {total_integral}")
print(f"The final value of alpha is: {alpha}\n")

print("2. Find the value of P(X < 3):")
print(f"The equation for the probability is: P(X < 3) = alpha * Integral(f(x) from 1 to 3)")
print(f"The value of the second integral is: {partial_integral}")
print(f"P(X < 3) = {alpha} * {partial_integral}")
print(f"The final value of P(X < 3) is: {prob_x_less_than_3}")
