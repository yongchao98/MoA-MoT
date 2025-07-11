import math

def print_m32_expression():
  """
  Prints the symbolic expression for the M_32 element of the RPR robot's inertia matrix.
  The problem asks to output each number in the final equation.
  The final expression is M_32 = -1 * m_3 * d_c3 * sin(q_3).
  """

  # Define the components of the expression
  coefficient = -1
  mass_symbol = "m_3"
  distance_symbol = "d_c3"
  trigonometric_function = "sin(q_3)"

  # Print the equation part by part
  print("The expression for the entry M_32 is given by:")
  print(f"M_32 = {coefficient} * {mass_symbol} * {distance_symbol} * {trigonometric_function}")

# Execute the function to display the result
print_m32_expression()