import numpy as np

def h(x):
  """
  Calculates the function h(x) for a given x.
  Note: This function is defined for x > 0.
  """
  if x <= 0:
    return float('nan')
  return 4*x**2 - 6*x + 2 + 2*x*np.log(2*x)

# The problem asks for the function h(x).
# We will print the formula for h(x) as requested.
# The formula is h(x) = 4*x^2 - 6*x + 2 + 2*x*ln(2*x)

print("The function h(x) is given by the formula:")
print("h(x) = 4*x**2 - 6*x + 2 + 2*x*log(2*x)")

# The coefficients and terms in the final equation are:
print("\nBreakdown of the final equation for h(x):")
print("Term 1: 4 * x**2")
print("Term 2: -6 * x")
print("Term 3: 2")
print("Term 4: 2 * x * log(2*x)")
