import math

def h(x):
  """
  Calculates the value of the function h(x) which defines the stability boundary.
  h(x) = 2 - 2*x - 4*x^2
  """
  return 2 - 2*x - 4*x**2

def solve():
  """
  Solves the problem by defining the function h(x) and printing its formula.
  The condition for a(t) to approach 0 is -sqrt(h(b(0))) < a(0) < 0.
  This is derived from analyzing the nullcline of a'(t). The region where a'(t) > 0,
  which pushes a(t) towards 0, is given by a^2 < 2 - 2*b - 4*b^2.
  Thus, h(b) = 2 - 2*b - 4*b^2.
  """

  # The coefficients of the polynomial h(x) = c0 + c1*x + c2*x^2
  c0 = 2
  c1 = -2
  c2 = -4

  print("The function h(x) determines the stability condition for a(0).")
  print("The condition is: -sqrt(h(b(0))) < a(0) < 0 ensures that a(t) -> 0.")
  print("\nThe equation for the function h(x) is:")
  print(f"h(x) = {c0} + ({c1})*x + ({c2})*x^2")
  
  # As an example, let's test for a specific value b(0)=0.25, which satisfies b(0) < 0.5
  b0 = 0.25
  h_val = h(b0)
  a0_limit = -math.sqrt(h_val)
  print(f"\nFor instance, if b(0) = {b0}:")
  print(f"h({b0}) = {h_val}")
  print(f"The condition on a(0) becomes: {a0_limit:.4f} < a(0) < 0")


solve()
