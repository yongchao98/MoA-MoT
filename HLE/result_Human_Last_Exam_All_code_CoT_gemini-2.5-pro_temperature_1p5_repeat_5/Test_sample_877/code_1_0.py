import math

def h(x):
  """
  This function defines the boundary for the initial condition a(0).
  If a(0)^2 < h(b(0)), the solution for a(t) is bounded.
  If a(0)^2 = h(b(0)), the solution converges to a=0.
  The function is valid for x > 0.
  """
  if x <= 0:
    return float('nan')
  return 4 * x**2 - 6 * x + 2 + 2 * x * math.log(2 * x)

# The final equation for h(x) is derived from the conserved quantity of the system.
# Let's print the formula for h(x) and show each number in it.
c4 = 4
c6 = 6
c2_1 = 2
c2_2 = 2
c2_3 = 2

# We output the formula for the function h(x).
print(f"h(x) = {c4}*x^2 - {c6}*x + {c2_1} + {c2_2}*x*ln({c2_3}*x)")