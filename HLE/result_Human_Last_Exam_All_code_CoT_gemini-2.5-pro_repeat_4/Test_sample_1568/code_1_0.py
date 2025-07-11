def get_product_formula():
  """
  Provides the closed-form expression for the infinite product as a string.
  The string represents the formula derived above and includes all the
  relevant numbers (1, 3, 8, 2, 4) as requested. The notation used
  (e.g., Gamma, exp, pi, i) is common in symbolic mathematics.
  """

  # The product starts at n=3. The terms for n=1 and n=2 are in the denominator.
  # Term for n=1: (1 - z**3 / 1**3) simplifies to (1 - z**3)
  # Term for n=2: (1 - z**3 / 2**3) simplifies to (1 - z**3 / 8)
  
  # The formula involves the Gamma function with arguments related to the cube roots of unity.
  # The cube roots of 1 are 1, exp(2*pi*i/3), and exp(4*pi*i/3).
  
  final_formula_string = "1 / ((1 - z**3) * (1 - z**3 / 8) * Gamma(1 - z) * Gamma(1 - z*exp(2*pi*i/3)) * Gamma(1 - z*exp(4*pi*i/3)))"
  
  print(final_formula_string)

get_product_formula()