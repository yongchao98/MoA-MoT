def generate_equation():
  """
  This function prints the equation for the given graph based on its features.
  The equation is constructed from the function's zeros and vertical asymptotes.
  """

  # Define the terms for the numerator based on the x-intercepts (-b, b, d)
  # The factors are (x+b), (x-b), and (x-d).
  # (x+b)(x-b) simplifies to (x^2 - b^2).
  numerator_term1 = "x^2 - b^2"
  numerator_term2 = "x - d"
  
  # Define the terms for the denominator based on the vertical asymptotes (a, c)
  # The factors are (x-a) and (x-c).
  denominator_term1 = "x - a"
  denominator_term2 = "x - c"

  # Assemble the final equation string.
  # We use extra parentheses to make the expression clear.
  equation = f"f(x) = (({numerator_term1})({numerator_term2})) / (({denominator_term1})({denominator_term2}))"

  print("The equation for the function f(x) is:")
  print(equation)

generate_equation()