import numpy as np
from scipy.integrate import quad

def simplified_integrand(x):
  """
  This is the simplified form of the original integrand: arctan(sin(x)) / sin(x).
  SciPy's quad function can handle the removable singularities at x=0 and x=pi.
  """
  # Avoid explicit division by zero if x is exactly 0 or pi, though quad is robust.
  if np.sin(x) == 0:
    return 1.0
  return np.arctan(np.sin(x)) / np.sin(x)

# The limits of the integral I
lower_limit = 0
upper_limit = np.pi

# Perform the numerical integration using quad
# The quad function returns the result of the integral and an estimate of the error
I, error_estimate = quad(simplified_integrand, lower_limit, upper_limit)

# The problem is to determine the value of I.
# The original equation is I = the integral. We print the determined value for I.
print(f"The value of the integral is I = {I}")
