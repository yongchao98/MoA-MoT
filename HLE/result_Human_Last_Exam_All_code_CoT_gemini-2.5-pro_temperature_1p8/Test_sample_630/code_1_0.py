import math

def display_rate_equation():
  """
  This function displays the derived optimal rate of convergence for the given problem.
  The rate is Theta(log(T) / sqrt(T)), which can be written as log(T) * T^(-0.5).
  This script prints the components of this equation.
  """
  term1 = "log(T)"
  term2 = "T**(-0.5)"
  
  # The final equation describes the dependency of the convergence rate on T.
  # The rate is not a simple power of T due to the logarithmic factor.
  print(f"The derived optimal rate of convergence is proportional to: {term1} * {term2}")
  
  # As requested, outputting the numbers from the final equation.
  # The exponent in the equation is -0.5.
  exponent = -0.5
  print(f"The exponent of T in the equation is: {exponent}")
  
display_rate_equation()