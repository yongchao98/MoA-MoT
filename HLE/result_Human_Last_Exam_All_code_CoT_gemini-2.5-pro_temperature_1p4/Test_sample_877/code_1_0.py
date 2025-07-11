import math

def h(x):
  """
  Calculates the function h(x) that defines the stability region.
  Note: This function is for demonstration and is valid for x > 0.
  """
  if x <= 0:
    return float('nan')
  # The mathematical function is h(x) = 4x^2 - 6x + 2 + 2x*ln(2x)
  # Here we print the formula as a string, as requested.
  
  # The coefficients and terms of the final equation are:
  c1 = 4
  c2 = -6
  c3 = 2
  c4 = 2
  c5 = 2
  
  # Output each number in the final equation.
  print(f"h(x) = {c1} * x**2 + ({c2}) * x + {c3} + {c4} * x * ln({c5} * x)")

# Let's call the function to print the result.
# The variable in the problem is b(0), but the function is requested in terms of x.
h(0.1)
