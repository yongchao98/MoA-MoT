def generate_equation():
  """
  This function analyzes the properties of the given graph and constructs the corresponding equation for f(x).
  """

  # Based on the analysis of the graph's features:
  # Numerator:
  # - Zeros at x = -b and x = b, crossing the axis -> factor (x+b)(x-b) = (x^2 - b^2)
  # - Zero at x = d, touching the axis (double root) -> factor (x-d)^2
  numerator = "(x^2 - b^2)(x - d)^2"

  # Denominator:
  # - Vertical asymptote at x = a, with sign change -> factor (x-a)
  # - Vertical asymptote at x = c, without sign change -> factor (x-c)^2
  denominator = "(x - a)(x - c)^2"

  # The combination gives a rational function with a slant asymptote, as deg(num) = deg(den) + 1.
  
  # Printing the final equation for f(x)
  print(f"Based on the graphical features, the equation for the function is:")
  print(f"f(x) = {numerator} / {denominator}")

generate_equation()