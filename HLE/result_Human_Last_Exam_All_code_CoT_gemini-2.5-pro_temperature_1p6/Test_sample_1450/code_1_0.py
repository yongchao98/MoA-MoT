import numpy as np
from scipy.special import polylog
from scipy.optimize import root_scalar

def fermi_dirac_integral(n, z):
  """
  Calculates the Fermi-Dirac integral f_n(z) using the polylogarithm function.
  f_n(z) = -Li_n(-z)
  """
  if z < 0:
    # Polylog is not defined for z > 1, which corresponds to -z < -1.
    # Our fugacity z should be positive.
    raise ValueError("Fugacity z must be non-negative.")
  return -polylog(n, -z)

def target_function(z):
  """
  The function for which we want to find the root.
  target_function(z) = f_{3/2}(z) / f_{5/2}(z) - 0.75
  """
  f_3_2 = fermi_dirac_integral(1.5, z)
  f_5_2 = fermi_dirac_integral(2.5, z)
  if f_5_2 == 0:
      # Avoid division by zero, although this won't happen for z > 0
      return np.inf
  return f_3_2 / f_5_2 - 0.75

# Find the root of the target function.
# The fugacity z must be positive. We need to provide a bracket [a, b]
# for the solver such that target_function(a) and target_function(b) have opposite signs.
# For small z, f_n(z) ~ z, so ratio ~ 1, target_function > 0.
# Let's test a larger z, e.g., z=10. target_function(10) is likely negative.
# So a bracket of [0.1, 10] should work.
try:
    solution = root_scalar(target_function, bracket=[0.1, 10], method='brentq')
    z_value = solution.root

    # Verification: output the numbers in the final equation
    f32_sol = fermi_dirac_integral(1.5, z_value)
    f52_sol = fermi_dirac_integral(2.5, z_value)
    ratio = f32_sol / f52_sol

    print("Solving the equation: f_{3/2}(z) / f_{5/2}(z) = 0.75")
    print(f"For the calculated fugacity z = {z_value:.4f}:")
    print(f"f_{3/2}(z) = {f32_sol:.4f}")
    print(f"f_{5/2}(z) = {f52_sol:.4f}")
    print(f"Ratio = {f32_sol:.4f} / {f52_sol:.4f} = {ratio:.4f}\n")

    # Print the final answer to two significant digits
    print(f"The value of the fugacity is: {'{:.2g}'.format(z_value)}")

except (ImportError, ValueError) as e:
    print(f"An error occurred: {e}")
    print("Please ensure you have scipy installed (`pip install scipy`).")
