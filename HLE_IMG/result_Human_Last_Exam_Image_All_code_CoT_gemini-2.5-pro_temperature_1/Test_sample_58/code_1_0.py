import math

def solve_cylinder_height():
  """
  This function prints the formula for the height of the cylinder.
  The formula is derived from the geometric properties of the cylinder when its surface is unrolled.
  """
  
  # The formula for the height 'h' in terms of radius 'r' and angle 'theta' is:
  # h = 2 * r * sqrt(pi * (pi - 2*theta))
  # We will print this formula in a readable format.
  
  # The numbers in the equation are 2 and 2.
  num1 = 2
  num2 = 2
  
  # Using unicode characters for mathematical symbols for better readability.
  radius_symbol = "r"
  pi_symbol = "π"
  theta_symbol = "θ"
  
  print(f"h = {num1} * {radius_symbol} * sqrt({pi_symbol} * ({pi_symbol} - {num2} * {theta_symbol}))")

solve_cylinder_height()