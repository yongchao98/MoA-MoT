import math

def print_rate_equation():
  """
  This function prints the components of the optimal convergence rate equation.
  The rate is of the form: Numerator / T^(Exponent)
  """
  
  # The optimal convergence rate for stochastic convex optimization is Theta(1/sqrt(T)).
  # This can be written as T^(-1/2).
  
  # The numerator in the rate equation (as a function of T)
  numerator = 1
  
  # The exponent of T in the denominator
  exponent_p = 0.5
  
  print("The optimal rate of convergence is of the form: Numerator / T^(p)")
  print(f"The number for the numerator is: {numerator}")
  print(f"The number for the exponent 'p' is: {exponent_p}")

print_rate_equation()
