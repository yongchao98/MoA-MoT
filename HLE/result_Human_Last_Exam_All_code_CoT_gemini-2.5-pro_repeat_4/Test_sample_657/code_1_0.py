def solve_controller():
  """
  This function formats and prints the derived transfer function for the set of
  all proper stabilizing controllers H_2(s).
  """
  
  # The numerator is represented as a combination of polynomials and the parameter K(s)
  # Numerator = 4*(s**2 + 2*s + 1) + (s**2 - 1)*K(s)
  # Numerator = 4*s**2 + 8*s + 4 + (s**2 - 1)*K(s)
  num_str = "4*s^2 + 8*s + 4 + (s^2 - 1)*K(s)"
  
  # The denominator is represented as a combination of polynomials and the parameter K(s)
  # Denominator = s**2 - 1 - s*K(s)
  den_str = "s^2 - 1 - s*K(s)"
  
  # Print the final transfer function H_2(s)
  print("The set of all proper stabilizing controllers H_2(s) is given by:")
  print("")
  print(f"H_2(s) = ({num_str}) / ({den_str})")
  print("")
  print("where K(s) is any stable and proper rational function.")

solve_controller()