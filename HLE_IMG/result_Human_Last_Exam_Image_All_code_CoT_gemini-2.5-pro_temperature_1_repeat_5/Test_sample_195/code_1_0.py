def generate_equation():
  """
  This function generates and prints the equation for the function f(x) based on its graphical features.
  The points a, b, c, d are treated as symbolic constants.
  """
  
  # The constants from the x-axis.
  a = "a"
  b = "b"
  c = "c"
  d = "d"

  # Construct the numerator based on the roots and their multiplicities.
  # Roots are at -b, b, and d (double root).
  numerator = f"(x + {b})(x - {b})(x - {d})^2"

  # Construct the denominator based on the vertical asymptotes and their behavior.
  # Asymptotes are at a (double pole) and c (single pole).
  denominator = f"(x - {a})^2(x - {c})"

  # Print the final equation for f(x).
  print(f"f(x) = {numerator} / {denominator}")

generate_equation()