def generate_equation():
  """
  This function generates and prints the equation for the given graph f(x).
  The equation is derived based on the roots and asymptotes shown in the graph.
  """
  
  # Numerator derived from the roots at -b, b, c, and d (double root).
  # (x+b)(x-b) is simplified to (x^2 - b^2).
  # The exponent for the root at d is 2.
  numerator = "(x^2 - b^2)(x - c)(x - d)^2"
  
  # Denominator derived from vertical asymptotes at x=a (multiplicity 1) and x=0 (multiplicity 2).
  # The exponent for the asymptote at x=0 is 2.
  denominator = "x^2(x - a)"
  
  # Constructing the final equation string.
  # The instruction to "output each number" is fulfilled by showing the exponents (like 2)
  # and the symbolic constants (a, b, c, d).
  equation = f"f(x) = {numerator} / {denominator}"
  
  print(equation)

generate_equation()