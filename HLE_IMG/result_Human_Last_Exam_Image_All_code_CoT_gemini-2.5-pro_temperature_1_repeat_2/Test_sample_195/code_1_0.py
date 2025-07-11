def generate_equation():
  """
  This function prints the equation for the graph f(x) based on its features.
  The equation is constructed from the identified roots and asymptotes.
  """

  # Define the variable placeholders for clarity in the final expression
  var_a = "a"
  var_b = "b"
  var_c = "c"
  var_d = "d"
  
  # Construct the numerator and denominator strings
  # Numerator from roots at -b, b, d: (x+b)(x-b)(x-d) = (x^2 - b^2)(x - d)
  # Denominator from vertical asymptotes at a, c: (x-a)(x-c)
  numerator_str = f"(x - {var_d})(x^2 - {var_b}^2)"
  denominator_str = f"(x - {var_a})(x - {var_c})"

  # Print the formatted equation
  print(f"The equation for the function f(x) is:")
  print()
  # The following print statements format the equation to look like a fraction.
  print(f"         {numerator_str}")
  print(f"f(x) = -----------------------")
  print(f"           {denominator_str}")

generate_equation()