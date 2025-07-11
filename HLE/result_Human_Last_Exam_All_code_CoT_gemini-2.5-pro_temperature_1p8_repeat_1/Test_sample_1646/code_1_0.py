import math

def get_function_form():
  """
  This function provides the explicit form of f(z) based on the derivation.
  """
  
  # The functional equation is f(z) = 2^(1 - z) * f(z/2) * f((z+1)/2)
  # with f(1) = sqrt(pi).
  # The derivation shows that f(z) is related to the reciprocal of the Gamma function.
  
  # The final derived form of the function.
  numerator = "sqrt(pi)"
  denominator = "Gamma(z)"
  
  # The equation f(z) = numerator / denominator
  # There are no explicit integer numbers in the final expression,
  # other than the implicit coefficient 1.
  
  final_equation = f"f(z) = {numerator} / {denominator}"
  
  print(final_equation)

get_function_form()