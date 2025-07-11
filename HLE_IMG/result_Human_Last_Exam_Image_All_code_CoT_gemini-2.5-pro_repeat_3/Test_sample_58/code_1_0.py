import math

def print_height_formula():
  """
  This function prints the derived formula for the height of the cylinder.
  """
  # The derived formula is h = pi * r * tan(theta)
  pi_symbol = 'π'
  r_symbol = 'r'
  theta_symbol = 'θ'
  
  # As requested, we output each number/symbol in the final equation.
  # The formula is symbolic, so we print the symbols themselves.
  # The only 'number' is pi.
  print(f"The height 'h' of the cylinder is given by the formula:")
  print(f"h = {pi_symbol} * {r_symbol} * tan({theta_symbol})")
  # Let's break down the components of the formula as requested by the prompt
  # "you still need to output each number in the final equation!"
  # The number in the equation is pi.
  print("\nBreaking down the equation:")
  print(f"The number in the equation is π (pi), which is approximately {math.pi}")
  print(f"h = {math.pi:.4f} * r * tan(θ)")


print_height_formula()