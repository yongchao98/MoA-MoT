import math

def print_formula():
  """
  This function prints the derived closed-form formula for P(n).
  The derivation involves identifying the given functions with the Gamma function,
  simplifying the expression using Gamma function identities (like the Legendre duplication formula),
  rearranging the product using number-theoretic sum identities, and finally
  evaluating the resulting sum using Gauss's multiplication formula.
  The final formula expresses P(n) in terms of n, pi, and the factorial function.
  """
  
  # The derived closed-form formula for P(n)
  formula = "P(n) = (2 * pi)**(n * (n - 1) / 4) / (n!)**0.5"
  
  # Printing the formula, showing each number involved as requested.
  # The numbers in the equation are 2, 1, 4, and 0.5.
  print("The closed-form formula for P(n) is:")
  print(formula)

print_formula()