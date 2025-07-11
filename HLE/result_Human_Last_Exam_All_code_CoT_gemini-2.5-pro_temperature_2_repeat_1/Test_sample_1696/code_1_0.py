import math

def final_formula():
  """
  This function prints the closed-form formula for P(n).
  The formula is derived step-by-step as outlined in the explanation.
  """
  n_placeholder = 'n'
  num_2 = 2
  num_1 = 1
  num_4 = 4

  # The formula is P(n) = (2 * pi)^(n * (n+1) / 4) / sqrt(n!)
  # The print statement is formatted to clearly show the components of the equation.
  print("The closed-form formula for P(n) is:")
  print(f"P({n_placeholder}) = ({num_2} * pi)^({n_placeholder} * ({n_placeholder} + {num_1}) / {num_4}) / sqrt({n_placeholder}!)")

final_formula()
