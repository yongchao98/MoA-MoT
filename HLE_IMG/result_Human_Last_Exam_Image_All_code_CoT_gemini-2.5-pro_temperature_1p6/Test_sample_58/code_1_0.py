import sympy

def solve_cylinder_height():
  """
  This function provides the symbolic formula for the height of the cylinder.
  """
  # Define the symbols used in the equation
  h, r, theta = sympy.symbols('h r theta')
  pi = sympy.pi

  # The derived equation is h = pi * r * tan(theta)
  # We will print this relationship in a readable format.
  # The problem asks to output each number, but there are no literal numbers,
  # so we represent the equation symbolically.
  
  # The mathematical constant pi is a 'number' in the equation.
  
  # Final Answer construction
  # h = (pi) * r * tan(theta)
  # The numbers/constants in the equation are essentially '1' for r and '1' for tan, and 'pi'.
  
  print("The height 'h' of the cylinder is given by the formula:")
  # We construct the string representation of the formula.
  # The problem states "output each number in the final equation" which implies printing the component parts.
  # Since it's symbolic, we will represent it as a clear mathematical expression.
  print("h = pi * r * tan(theta)")
  # To satisfy the "each number" constraint in a symbolic context, let's also point them out.
  # The implicit numerical coefficients are 1. The constant is pi.
  # There are no explicit numbers other than the constant pi.
  # A more detailed print might look like:
  # print(f"h = {1} * pi * {1} * r * tan(theta)") 
  # But this is less natural. The primary result is the formula itself.

solve_cylinder_height()