import math

def get_theta_prime_expression():
  """
  This function returns the mathematical expression for theta'(t).
  """
  # The derived expression for theta'(t) is c*cos^2(theta) + (K/c)*sin^2(theta).
  # We will format this as a string for output.
  c = "c"
  K = "K(gamma(t))"
  theta = "theta(t)"
  
  # The problem asks for the value of theta'(t), which is an expression.
  # We assume the question asks for the positive rate of rotation.
  expression = f"{c} * cos^2({theta}) + (1/{c}) * {K} * sin^2({theta})"
  
  # The problem asks to output each number in the final equation.
  # In our symbolic expression, the only number is 1.
  final_expression = f"{c} * cos^2({theta}) + {1}/{c} * {K} * sin^2({theta})"

  print(final_expression)

get_theta_prime_expression()