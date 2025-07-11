import math

def solve_for_sliding_angle():
  """
  This function provides the symbolic expression for the angle at which the rod begins to slide.
  The problem is symbolic, so no numerical values are required.
  """
  
  # The final equation derived from the physics principles is tan(theta) = mu.
  # where:
  # theta is the angle of the rod with the table.
  # mu is the coefficient of static friction.
  
  # The expression for theta is the arctangent of mu.
  final_expression = "theta = arctan(mu)"
  
  print("The final expression for the angle theta at which the rod begins to slide is:")
  print(final_expression)

solve_for_sliding_angle()