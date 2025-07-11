import math

def h(x):
  """
  Calculates the value of the function h(x) which defines the stability boundary.
  The function is defined for x > 0.
  """
  if x <= 0:
    return float('nan') # ln(2x) is not defined for x <= 0
  return 4*x**2 - 6*x + 2 + 2*x*math.log(2*x)

# The problem asks to determine the function h(x).
# The final code will print the formula for h(x), explicitly showing each numerical coefficient.

c_x2 = 4
c_x1 = -6
c_const = 2
c_log_term = 2
c_log_arg = 2

print("The derived function h(x) is:")
print(f"h(x) = ({c_x2})*x**2 + ({c_x1})*x + ({c_const}) + ({c_log_term})*x*ln(({c_log_arg})*x)")
