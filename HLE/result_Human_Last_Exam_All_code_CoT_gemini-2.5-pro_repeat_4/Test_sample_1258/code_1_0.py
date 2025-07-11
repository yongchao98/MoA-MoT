import math

def print_demagnetizing_factor_expression():
  """
  Prints the analytical expression for the fluxmetric demagnetizing factor
  for a cylinder.
  """

  # Define the components of the expression as strings
  term1 = "1"
  factor = "(2 * g / pi)"
  elliptic_term_E = "E(k)"
  elliptic_term_F = "F(k)"
  k_term = "(1 - k**2) / 2"

  # Construct the final expression string
  expression = f"N_f = {term1} - {factor} * ({elliptic_term_E} - {k_term} * {elliptic_term_F})"

  # Print the explanation and the final formula
  print("The analytical expression for the fluxmetric demagnetizing factor (N_f) is:")
  print(expression)
  print("\nWhere:")
  print("  g = length-to-diameter ratio")
  print("  pi = The mathematical constant ~3.14159")
  print("  k = The modulus, defined as k**2 = 1 / (1 + g**2 / 4)")
  print("  F(k) = The complete elliptic integral of the first kind")
  print("  E(k) = The complete elliptic integral of the second kind")

print_demagnetizing_factor_expression()